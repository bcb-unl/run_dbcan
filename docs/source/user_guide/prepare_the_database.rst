.. _prepare-database:

Preparing Databases
===================

Download all databases once into one folder and reuse it for every run:

.. code-block:: shell

   run_dbcan database --db_dir db               # HTTP, moving snapshot "db_current"
   run_dbcan database --db_dir db --aws_s3      # pinned S3 release: reproducible, often faster
   run_dbcan database --db_dir db --no-cgc      # CAZyme annotation databases only

A full download needs about **8 GB** (about 7.4 GB with ``--no-cgc``). The files are described in
:doc:`../getting_started/database_description`.

HTTP snapshot or pinned S3 release
----------------------------------

.. list-table::
   :widths: 22 78
   :header-rows: 1

   * - Source
     - Behaviour
   * - HTTP (default)
     - ``https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/<file>``.
       ``db_current`` is replaced when the dbCAN team publishes an update, so two downloads months apart may differ.
   * - ``--aws_s3``
     - ``https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/<file>``. The folder name is fixed in the
       installed package (``S3_DATABASE_RELEASE`` in ``dbcan/constants/databases_constants.py``), so one ``dbcan``
       version always downloads the same files.

Use ``--aws_s3`` when results must be reproducible across machines or over time.

Options
-------

.. list-table::
   :widths: 32 68
   :header-rows: 1

   * - Option
     - Description
   * - ``--db_dir``
     - Target directory (created if needed). Required.
   * - ``--cgc`` / ``--no-cgc``
     - Download (default) or skip the CGC and substrate databases (TCDB, TF, STP, peptidase, sulfatase, dbCAN-PUL).
   * - ``--aws_s3``
     - Pinned S3 release instead of the HTTP ``db_current`` snapshot.
   * - ``--no-overwrite``
     - Skip every file that already exists, without checking its size.
   * - ``--resume`` / ``--no-resume``
     - Resume an interrupted download from ``<file>.part`` (default) or start again.
   * - ``--timeout`` / ``--retries``
     - Per-request timeout in seconds (default 30) and retries of transient failures (default 3).
   * - ``--verify-ssl`` / ``--no-verify-ssl``
     - Verify TLS certificates (default). ``--no-verify-ssl`` is insecure.

* Without ``--no-overwrite`` an existing file is skipped only when its size matches the remote ``Content-Length``.
* ``dbCAN-PUL.tar.gz`` is extracted into ``dbCAN-PUL/`` and deleted.
* If a file fails, the command exits with an error listing it; run it again to fetch only the missing files.
  On a slow network: ``run_dbcan database --db_dir db --aws_s3 --timeout 60 --retries 5``.

Files and sizes
---------------

Sizes of the S3 release ``db_v5-2-9_5-5-2026`` (the HTTP snapshot may differ slightly).

.. list-table::
   :widths: 34 12 54
   :header-rows: 1

   * - File
     - Size
     - Used by
   * - ``dbCAN-sub.hmm`` (remote ``dbCAN_sub.hmm``)
     - ~5.1 GB
     - ``CAZyme_annotation``
   * - ``CAZy.dmnd``
     - ~2.2 GB
     - ``CAZyme_annotation``
   * - ``dbCAN.hmm``
     - ~130 MB
     - ``CAZyme_annotation``
   * - ``fam-substrate-mapping.tsv``
     - <1 MB
     - ``CAZyme_annotation`` (``Substrate`` column); the 5-column file and the curated web format (comment line,
       ``new_Substrate_high_level`` and ``type`` columns, ``REMOVED`` rows ignored) are both accepted
   * - ``dbCAN-sub-substrate.tsv`` (optional, not downloaded)
     - <1 MB
     - ``CAZyme_annotation``: when present, substrates of subfamilies with EC numbers match the dbCAN-sub web
       server exactly
   * - ``peptidase_db.dmnd``
     - ~330 MB
     - ``gff_process``
   * - ``sulfatlas_db.dmnd``
     - ~100 MB
     - ``gff_process``
   * - ``TCDB.dmnd``, ``TF.hmm``, ``TF.dmnd``, ``STP.hmm``
     - <15 MB each
     - ``gff_process``
   * - ``PUL.dmnd``, ``dbCAN-PUL.xlsx``, ``dbCAN-PUL/``
     - <5 MB
     - ``substrate_prediction``
   * - ``Pfam-A.hmm`` (optional, not downloaded)
     - ~1.6 GB
     - ``Pfam_null_cgc`` only; download it by hand (below)

Recording the database version
------------------------------

.. code-block:: shell

   run_dbcan version
   ls -l db/
   md5sum db/*.dmnd db/*.hmm > db/checksums.md5

Report the S3 release name (``db_v5-2-9_5-5-2026``) if you used ``--aws_s3``, otherwise the download date.

Manual download
---------------

Use the same file names as the downloader; only ``dbCAN_sub.hmm`` is saved under another name (``dbCAN-sub.hmm``,
the name run_dbcan looks for).

.. code-block:: shell

   mkdir -p db && cd db
   BASE="https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current"
   # or: BASE="https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026"
   for f in CAZy.dmnd dbCAN.hmm fam-substrate-mapping.tsv \
            TCDB.dmnd TF.hmm TF.dmnd STP.hmm peptidase_db.dmnd sulfatlas_db.dmnd \
            PUL.dmnd dbCAN-PUL.xlsx dbCAN-PUL.tar.gz; do
     wget -O "$f" "$BASE/$f"
   done
   wget -O dbCAN-sub.hmm "$BASE/dbCAN_sub.hmm"
   tar -xzf dbCAN-PUL.tar.gz

   # optional, only for Pfam_null_cgc
   wget -O Pfam-A.hmm "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/Pfam-A.hmm"

Or copy the whole S3 release (about 12 GB: it also holds ``Pfam-A.hmm`` and the CAZy FASTA ``CAZyDB.<date>.fa``):

.. code-block:: shell

   aws s3 cp s3://dbcan/db_v5-2-9_5-5-2026/ ./db --recursive --no-sign-request
   mv db/dbCAN_sub.hmm db/dbCAN-sub.hmm

Test files are in the `dbCAN2 test data folder <https://bcb.unl.edu/dbCAN2/download/test>`_; the
:doc:`../getting_started/quick_start` uses them.

.. _rebuild-peptidase-db:

Rebuilding ``peptidase_db.dmnd`` (maintainers)
----------------------------------------------

Users download the finished file. For a new MEROPS release the dbCAN team rebuilds it from ``pepunit.lib``
(DIAMOND must be on ``PATH``):

.. code-block:: shell

   wget https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib     # ~440 MB
   python -m dbcan.utils.build_peptidase_db \
       --input pepunit.lib --output peptidase_db.dmnd \
       --report peptidase_db.dropped.tsv --threads 4

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Step
     - What the script does
   * - Identifier
     - Read from the bracket right before ``#<family>#{... unit: ...}``:
       ``>MER0003996 - ubiquitin-specific endopeptidase 39 [misleading] (Mus musculus) [C19.972]#C19#{peptidase unit:
       187-549}`` gives ``C19.972``. The currently distributed database took the *first* bracket and has 322 entries
       named ``misleading``, ``Clostridium``, ..., which run_dbcan never reports.
   * - Record name
     - ``>MER0003996|C19.972``; run_dbcan keeps the part after ``|`` as ``Peptidase ID``.
   * - Kept
     - Peptidase units and, by default, inhibitor units (families ``I*``) as decoys (``--drop-inhibitors`` leaves
       them out).
   * - Left out
     - Identifiers not of the form ``<type><family>.<id>`` (e.g. ``a01.UNa``), records without ``MER`` accession or
       sequence, exact duplicates (first accession kept). Counts are logged; ``--report`` lists every record.
   * - Output
     - Cleaned FASTA (``--fasta``, default ``peptidase_db.faa``) and ``diamond makedb`` (skip with
       ``--no-makedb``). Non-zero exit when no peptidase unit is found or ``makedb`` fails.

MEROPS release of 2023-02-22 (1,227,939 records): 984,975 peptidase units (287 families) and 147,401 inhibitor
units are kept; the database is about 313 MB (about 290 MB with ``--drop-inhibitors``).

.. note::

   **Why keep inhibitor units?** Each protein keeps only its best DIAMOND hit (``--max-target-seqs 1``) and hits to
   inhibitor units are never reported as ``Peptidase``. With the decoys, a protein whose best match is an inhibitor
   cannot take a weaker peptidase hit instead: on *Xylona heveae* one protein whose best hit is ``I51.001``
   (carboxypeptidase Y inhibitor) is labelled ``C18.UPW`` (a viral peptidase family) when they are dropped. They also
   keep DIAMOND in its fast query-indexed mode on small genomes.

Next: :doc:`CAZyme_annotation`.
