.. _signalp6-deeptmhmm:

SignalP 6.0 and DeepTMHMM (optional tools)
==========================================

``run_dbcan`` can add **signal peptide** and **transmembrane topology** columns to ``overview.tsv`` during ``CAZyme_annotation``. These features rely on **SignalP 6.0** and/or **DeepTMHMM**, which are **not** distributed or installed as dependencies of the ``dbcan`` package. You must install them yourself, verify that they run in your environment, and only then enable the corresponding flags.

Why separate installation?
---------------------------

* **Licensing and distribution**: Upstream tools have their own licenses and release models; we integrate them only as optional wrappers.
* **Heavy dependencies**: SignalP 6.0 and DeepTMHMM typically need PyTorch and other large stacks; keeping them optional avoids forcing every user to install them.
* **Environment-specific failures**: GPU/CPU, ``libstdc++``, and PyTorch versions differ across clusters; local testing is essential.

SignalP 6.0
-------------

**Role in run_dbcan**

When you pass ``--run_signalp`` to ``CAZyme_annotation``, the pipeline invokes the ``signalp6`` executable, writes results under the output directory, and fills the **SignalP** column in ``overview.tsv``.

**Install**

Follow the official SignalP 6.0 installation instructions from the upstream project (e.g. `installation instructions <https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md>`_). Ensure the ``signalp6`` command is on your ``PATH`` in the same environment you use to run ``run_dbcan``.

**Recommended checks before using run_dbcan**

1. Run ``signalp6 --help`` (or your install’s equivalent) and confirm the CLI starts without import errors.
2. Run a **small FASTA** through SignalP 6.0 alone and confirm you get the expected output files.

**Common environment issue (Linux)**

If you see errors such as ``CXXABI_1.3.15`` missing when importing **matplotlib** inside SignalP, the dynamic linker may be picking an older system ``libstdc++`` instead of the one from your conda environment. Typical mitigations:

* Activate the conda environment that contains SignalP and run ``run_dbcan`` from that same environment; the pipeline prepends that environment’s ``lib`` directory to ``LD_LIBRARY_PATH`` for the SignalP subprocess when possible.
* If you still invoke ``signalp6`` manually in a bare shell, try::

     export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

  (adjust if you are not using conda).

**run_dbcan options (SignalP)**

* ``--run_signalp`` / ``--no-run_signalp`` — enable or disable SignalP.
* ``--signalp_org`` — ``other`` or ``euk`` (organism class passed to SignalP).

DeepTMHMM
----------

**Role in run_dbcan**

When you pass ``--run_deeptmhmm`` and ``--deeptmhmm_dir``, the pipeline runs a local **DeepTMHMM** ``predict.py`` from that directory (with working directory set there), writes outputs to a ``deeptmhmm_out`` folder under your ``CAZyme_annotation`` output directory, and fills the **DeepTMHMM** column in ``overview.tsv``.

**Install**

Obtain a DeepTMHMM (or compatible) checkout that contains ``predict.py`` and its model/assets, as required by the upstream tool. Point ``--deeptmhmm_dir`` at the **directory that contains** ``predict.py``.

**Recommended checks before using run_dbcan**

1. ``cd`` into that directory and run your chosen Python on ``predict.py --help`` (or the upstream documented test command).
2. Run **one short FASTA** through ``predict.py`` with a fresh output directory and confirm topology files are produced.

**PyTorch 2.6+ and saved checkpoints**

Some DeepTMHMM checkpoints are pickled objects (not “weights only”). If you see ``weights_only`` / ``UnpicklingError`` when loading ``esm_model_args.pt`` or similar, use a DeepTMHMM tree that loads trusted checkpoints with ``weights_only=False`` (or follow PyTorch’s guidance for allowlisted globals). Only load checkpoints you trust.

**run_dbcan options (DeepTMHMM)**

* ``--run_deeptmhmm`` / ``--no-run_deeptmhmm`` — enable or disable DeepTMHMM.
* ``--deeptmhmm_dir`` — path to the DeepTMHMM install (directory containing ``predict.py``).
* ``--deeptmhmm_python`` — Python executable to use (default: ``python``). Use this if DeepTMHMM must run in a specific conda env interpreter.

**Column overwrite**

Use ``--force_topology`` if you want to **replace** existing **SignalP** / **DeepTMHMM** values in ``overview.tsv`` instead of only filling empty cells.

Example (both tools enabled)
-----------------------------

After you have verified SignalP and DeepTMHMM independently, you can run (adjust paths)::

   run_dbcan CAZyme_annotation \
     --input_raw_data proteins.faa \
     --output_dir out \
     --db_dir db \
     --mode protein \
     --run_signalp \
     --run_deeptmhmm \
     --deeptmhmm_dir /path/to/DeepTMHMM

See also :doc:`CAZyme_annotation` for the full annotation workflow.
