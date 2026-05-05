from pathlib import Path

import pytest

from dbcan.configs.database_config import DBDownloaderConfig
from dbcan.utils.database import DBDownloader, DatabaseDownloadError


class _FakeSession:
    def close(self):
        pass


def test_download_file_raises_when_any_file_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(DBDownloader, "databases", property(lambda self: {"CAZy.dmnd": "https://example/db"}))
    monkeypatch.setattr(DBDownloader, "_prepare_session", lambda self: _FakeSession())

    def fail_download(self, session, url, output_path):
        raise RuntimeError("network down")

    monkeypatch.setattr(DBDownloader, "_download_single_file", fail_download)

    downloader = DBDownloader(DBDownloaderConfig(db_dir=str(tmp_path)))

    with pytest.raises(DatabaseDownloadError) as exc_info:
        downloader.download_file()

    assert "CAZy.dmnd" in exc_info.value.failures
    assert "network down" in exc_info.value.failures["CAZy.dmnd"]


def test_no_overwrite_skips_existing_file(monkeypatch, tmp_path):
    existing = tmp_path / "CAZy.dmnd"
    existing.write_text("already present")

    monkeypatch.setattr(DBDownloader, "databases", property(lambda self: {"CAZy.dmnd": "https://example/db"}))
    monkeypatch.setattr(DBDownloader, "_prepare_session", lambda self: _FakeSession())

    def unexpected_download(self, session, url, output_path):
        raise AssertionError("download should be skipped")

    monkeypatch.setattr(DBDownloader, "_download_single_file", unexpected_download)

    downloader = DBDownloader(DBDownloaderConfig(db_dir=str(tmp_path), no_overwrite=True))
    downloader.download_file()

    assert existing.read_text() == "already present"


def test_extract_tar_file_rejects_invalid_archive(tmp_path):
    tar_path = tmp_path / "dbCAN-PUL.tar.gz"
    tar_path.write_text("not a tar archive")

    downloader = DBDownloader(DBDownloaderConfig(db_dir=str(tmp_path)))

    with pytest.raises(ValueError, match="not a valid tar archive"):
        downloader._extract_tar_file(tar_path, Path(tmp_path))
