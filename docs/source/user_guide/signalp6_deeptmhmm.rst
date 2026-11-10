.. _signalp6-deeptmhmm:

SignalP 6.0 and DeepTMHMM (optional)
====================================

``CAZyme_annotation`` can add a **SignalP** (signal peptide) and a **DeepTMHMM** (transmembrane topology) column to
``overview.tsv``. Both tools are **not** installed with ``dbcan`` (own licences, large PyTorch stacks): install them,
test them on a small FASTA, then turn them on.

.. code-block:: shell

   run_dbcan CAZyme_annotation --mode protein --input_raw_data proteins.faa --output_dir out --db_dir db \
     --run_signalp --signalp_org euk \
     --run_deeptmhmm --deeptmhmm_dir /path/to/DeepTMHMM

* Only proteins with ``#ofTools`` >= 2 are analysed (written to ``topology_input_recommended.faa``); other rows get
  ``-``.
* The flags work alone or together. Outputs: ``signalp6_out/`` and ``deeptmhmm_out/`` (``_1``, ``_2``, ... if it
  exists) in ``--output_dir``.
* The columns are recomputed on every run; ``--force_topology`` is deprecated and ignored.

.. list-table::
   :widths: 36 64
   :header-rows: 1

   * - Option
     - Description (default)
   * - ``--run_signalp``
     - Run ``signalp6`` (fast mode); it must be on ``PATH`` (off)
   * - ``--signalp_org``
     - ``other`` (bacteria/archaea) or ``euk`` (``other``)
   * - ``--run_deeptmhmm``
     - Run DeepTMHMM ``predict.py`` (off); needs ``--deeptmhmm_dir``
   * - ``--deeptmhmm_dir``
     - Folder that contains ``predict.py`` and its model files; used as working directory
   * - ``--deeptmhmm_python``
     - Python that runs ``predict.py``, e.g. of another conda env (``python``)

**Failures.** A requested tool that cannot start (``signalp6`` not on ``PATH``, no ``predict.py`` in
``--deeptmhmm_dir``) stops the command with a usage error (exit status 2) *before* the CAZyme searches. A tool that
fails while running makes the command exit non-zero after the CAZyme results are written; the column of the other
tool is still added. Check the log and the output folders.

Installing and testing
----------------------

**SignalP 6.0**: follow the `SignalP 6.0 installation instructions
<https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md>`_ in the environment you run
``run_dbcan`` from, then check that ``signalp6 --help`` starts and that a small FASTA gives output files.

**DeepTMHMM**: get a checkout with ``predict.py`` and its models, then run ``python predict.py --help`` and one short
FASTA from that folder.

.. note::

   **Common errors.**

   * ``CXXABI_1.3.15`` missing (Linux, inside SignalP's matplotlib): an old system ``libstdc++`` is loaded. Run
     ``run_dbcan`` from the conda environment that holds SignalP (run_dbcan puts its ``lib`` first on
     ``LD_LIBRARY_PATH`` for SignalP); in a bare shell use
     ``export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"``.
   * ``weights_only`` / ``UnpicklingError`` with PyTorch >= 2.6 when DeepTMHMM loads ``esm_model_args.pt``: the
     checkpoint is a pickled object. Use a DeepTMHMM tree that loads trusted checkpoints with
     ``weights_only=False`` (only for checkpoints you trust).
