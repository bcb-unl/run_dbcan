import os
from pathlib import Path

import pytest


def _make_minimal_files(tmp_path: Path):
    # Minimal FASTA with 5 proteins
    faa = tmp_path / "input.faa"
    faa.write_text(
        "\n".join(
            [
                ">p1\nAAAA\n",
                ">p2\nAAAA\n",
                ">p3\nAAAA\n",
                ">p4\nAAAA\n",
                ">p5\nAAAA\n",
            ]
        )
    )
    # Minimal placeholder HMM file (content not parsed because we monkeypatch HMMFile)
    hmm = tmp_path / "dbcan.hmm"
    hmm.write_text("# dummy hmm\n")
    return faa, hmm


class _FakeSeqBlock(list):
    """Acts like a DigitalSequenceBlock for our purposes (truthy + len())."""


class _FakeSequenceFile:
    def __init__(self, path, digital=True, alphabet=None):
        self._i = 0
        self._n = 5

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read_block(self, sequences=None, residues=None):
        if sequences is None:
            sequences = self._n
        if self._i >= self._n:
            return _FakeSeqBlock()
        j = min(self._n, self._i + int(sequences))
        block = _FakeSeqBlock(range(self._i, j))
        self._i = j
        return block


class _FakeHMMFile:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def __iter__(self):
        # Iterable placeholder; we never yield real HMMs because fake hmmsearch ignores it.
        return iter([])


def test_batching_uses_read_block_sequences(monkeypatch, tmp_path):
    # Import inside test to avoid importing unrelated heavy modules at collection time
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig

    faa, hmm = _make_minimal_files(tmp_path)
    # Make the HMM file appear "large" (>1GB) so large_mode disables preload_hmms,
    # forcing the per-batch streaming HMMFile path.
    os.truncate(hmm, 1100 * 1024 * 1024)

    # Ensure processor finds the expected file names
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    # default expected by config/constants in repo
    (out_dir / "uniInput.faa").write_text(faa.read_text())

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "dbcan.hmm").write_text(hmm.read_text())

    calls = {"read_block": [], "hmmsearch": 0, "hmmfile_open": 0}

    def fake_hmmsearch(hmms_or_hmmfile, targets, cpus=1, domE=None):
        calls["hmmsearch"] += 1
        # Yield no hits
        if False:
            yield None
        return iter(())

    def fake_seqfile_ctor(path, digital=True, alphabet=None):
        return _FakeSequenceFile(path, digital=digital, alphabet=alphabet)

    def fake_hmmfile_ctor(path):
        calls["hmmfile_open"] += 1
        return _FakeHMMFile(path)

    monkeypatch.setattr("pyhmmer.hmmsearch", fake_hmmsearch)
    monkeypatch.setattr("pyhmmer.easel.SequenceFile", fake_seqfile_ctor)
    monkeypatch.setattr("pyhmmer.plan7.HMMFile", fake_hmmfile_ctor)

    # Track read_block sequences argument by wrapping method
    orig_read_block = _FakeSequenceFile.read_block

    def tracked_read_block(self, sequences=None, residues=None):
        calls["read_block"].append(sequences)
        return orig_read_block(self, sequences=sequences, residues=residues)

    monkeypatch.setattr(_FakeSequenceFile, "read_block", tracked_read_block, raising=True)

    cfg = PyHMMERDBCANConfig(
        db_dir=str(db_dir),
        output_dir=str(out_dir),
        threads=1,
        hmm_file="dbcan.hmm",
        batch_size=2,  # force batching
        large_mode=True,  # force safe streaming path
        enable_memory_monitoring=False,
    )
    p = PyHMMERDBCANProcessor(cfg)
    p.hmmsearch()

    # With 5 seqs and batch_size=2 => expected blocks: 2,2,1 then empty
    assert calls["read_block"][:3] == [2, 2, 2]
    assert calls["hmmsearch"] >= 3
    # In streaming HMM mode, HMMFile should be opened per batch
    assert calls["hmmfile_open"] >= 3


def test_retry_on_memoryerror(monkeypatch, tmp_path):
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig

    faa, hmm = _make_minimal_files(tmp_path)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "uniInput.faa").write_text(faa.read_text())

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "dbcan.hmm").write_text(hmm.read_text())

    calls = {"hmmsearch": 0}

    def flaky_hmmsearch(hmms_or_hmmfile, targets, cpus=1, domE=None):
        calls["hmmsearch"] += 1
        if calls["hmmsearch"] == 1:
            raise MemoryError("simulated OOM")
        if False:
            yield None
        return iter(())

    monkeypatch.setattr("pyhmmer.hmmsearch", flaky_hmmsearch)
    monkeypatch.setattr("pyhmmer.easel.SequenceFile", lambda *a, **k: _FakeSequenceFile(*a, **k))
    monkeypatch.setattr("pyhmmer.plan7.HMMFile", lambda p: _FakeHMMFile(p))

    cfg = PyHMMERDBCANConfig(
        db_dir=str(db_dir),
        output_dir=str(out_dir),
        threads=1,
        hmm_file="dbcan.hmm",
        batch_size=2,
        large_mode=True,
        enable_memory_monitoring=False,
        max_retries=2,
    )
    p = PyHMMERDBCANProcessor(cfg)
    p.hmmsearch()

    assert calls["hmmsearch"] >= 2, "Should retry after MemoryError"


def test_auto_large_mode_triggers_batching(monkeypatch, tmp_path):
    """
    Simulate a "large input" by making the FASTA file size exceed the threshold,
    and ensure the code enters the batching path even if batch_size is None.
    """
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig

    faa, hmm = _make_minimal_files(tmp_path)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    fasta_out = out_dir / "uniInput.faa"
    fasta_out.write_text(faa.read_text())
    # Make input appear > 1MB so it exceeds the low threshold below.
    os.truncate(fasta_out, 2 * 1024 * 1024)

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "dbcan.hmm").write_text(hmm.read_text())

    calls = {"read_block": []}

    def fake_hmmsearch(hmms_or_hmmfile, targets, cpus=1, domE=None):
        if False:
            yield None
        return iter(())

    monkeypatch.setattr("pyhmmer.hmmsearch", fake_hmmsearch)
    monkeypatch.setattr("pyhmmer.easel.SequenceFile", lambda *a, **k: _FakeSequenceFile(*a, **k))
    monkeypatch.setattr("pyhmmer.plan7.HMMFile", lambda p: _FakeHMMFile(p))

    orig_read_block = _FakeSequenceFile.read_block

    def tracked_read_block(self, sequences=None, residues=None):
        calls["read_block"].append(sequences)
        return orig_read_block(self, sequences=sequences, residues=residues)

    monkeypatch.setattr(_FakeSequenceFile, "read_block", tracked_read_block, raising=True)

    cfg = PyHMMERDBCANConfig(
        db_dir=str(db_dir),
        output_dir=str(out_dir),
        threads=1,
        hmm_file="dbcan.hmm",
        large_mode=False,
        large_input_threshold_mb=1,  # very low threshold to trigger auto large_mode
        batch_size=None,
        enable_memory_monitoring=False,
    )
    PyHMMERDBCANProcessor(cfg).hmmsearch()

    assert any(isinstance(x, int) and x > 0 for x in calls["read_block"]), "Expected batched read_block(sequences=N)"


def test_very_large_input_sparse_file_no_oom(monkeypatch, tmp_path):
    """
    Create a sparse >10GB FASTA to ensure the large-input heuristics kick in
    (large_mode auto-enable + batching + no aggressive preloading), without actually
    consuming disk/memory in the test environment.
    """
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig

    faa, hmm = _make_minimal_files(tmp_path)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    fasta_out = out_dir / "uniInput.faa"
    fasta_out.write_text(faa.read_text())
    # Make it appear ~11GB
    os.truncate(fasta_out, 11 * 1024 * 1024 * 1024)

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    hmm_path = db_dir / "dbcan.hmm"
    hmm_path.write_text(hmm.read_text())
    # Make HMM appear >1GB so large_mode disables preload_hmms
    os.truncate(hmm_path, 1100 * 1024 * 1024)

    calls = {"read_block": 0, "hmmfile_open": 0}

    def fake_hmmsearch(hmms_or_hmmfile, targets, cpus=1, domE=None):
        if False:
            yield None
        return iter(())

    monkeypatch.setattr("pyhmmer.hmmsearch", fake_hmmsearch)
    monkeypatch.setattr("pyhmmer.easel.SequenceFile", lambda *a, **k: _FakeSequenceFile(*a, **k))

    def counted_hmmfile(path):
        calls["hmmfile_open"] += 1
        return _FakeHMMFile(path)

    monkeypatch.setattr("pyhmmer.plan7.HMMFile", counted_hmmfile)

    orig_read_block = _FakeSequenceFile.read_block

    def tracked_read_block(self, sequences=None, residues=None):
        calls["read_block"] += 1
        return orig_read_block(self, sequences=sequences, residues=residues)

    monkeypatch.setattr(_FakeSequenceFile, "read_block", tracked_read_block, raising=True)

    cfg = PyHMMERDBCANConfig(
        db_dir=str(db_dir),
        output_dir=str(out_dir),
        threads=1,
        hmm_file="dbcan.hmm",
        batch_size=2,  # keep batches small in test
        large_mode=False,  # should auto-enable due to 11GB input size
        enable_memory_monitoring=False,
    )
    PyHMMERDBCANProcessor(cfg).hmmsearch()

    assert calls["read_block"] > 0
    assert calls["hmmfile_open"] >= 1

