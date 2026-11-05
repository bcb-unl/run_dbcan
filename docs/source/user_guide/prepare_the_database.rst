.. _prepare-database:

Preparing Databases
====================

Database files are essential for running dbCAN analysis. This guide explains how to install and prepare these databases.

Automatic Database Installation
----------------------------------

The simplest way to obtain all required databases is using the built-in command in run_dbCAN:

.. code-block:: shell

   run_dbcan database --db_dir db

To download from AWS S3 (faster and more stable in many regions), add the ``--aws_s3`` flag:

.. code-block:: shell

   run_dbcan database --db_dir db --aws_s3

To download only CAZyme-related databases (without CGC-related databases), use ``--no-cgc``:

.. code-block:: shell

   run_dbcan database --db_dir db --no-cgc

This command will automatically download all necessary databases and organize them in the specified directory.

HTTP snapshot vs. pinned S3 release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two download sources are not equivalent:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Source
     - Behaviour
   * - HTTP (default)
     - ``https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/<filename>``.
       ``db_current`` is a **moving snapshot**: it is replaced when the dbCAN team publishes a database update,
       so two downloads made months apart may differ.
   * - AWS S3 (``--aws_s3``)
     - ``https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_11-3-2026/<filename>``.
       A **pinned release** whose folder name is fixed in the installed ``dbcan`` package
       (``S3_DATABASE_RELEASE`` in ``dbcan/constants/databases_constants.py``). The same ``dbcan`` version always
       downloads the same files.

Use ``--aws_s3`` when you need reproducible database contents across machines or over time.

Downloader options
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Option
     - Description
   * - ``--db_dir``
     - Target directory (created if needed). Required.
   * - ``--cgc`` / ``--no-cgc``
     - Download (default) or skip the CGC and substrate-prediction databases (TCDB, TF, STP, peptidase, sulfatase, dbCAN-PUL).
   * - ``--aws_s3``
     - Use the pinned S3 release instead of the HTTP ``db_current`` snapshot.
   * - ``--no-overwrite``
     - Skip every file that already exists in ``--db_dir``, without checking its size.
   * - ``--resume`` / ``--no-resume``
     - Resume an interrupted download from its ``<filename>.part`` file (default), or always start from scratch.
   * - ``--timeout``
     - Per-request timeout in seconds (default: 30).
   * - ``--retries``
     - Retries for transient HTTP/S3 failures (default: 3).
   * - ``--verify-ssl`` / ``--no-verify-ssl``
     - Verify TLS certificates (default). ``--no-verify-ssl`` is insecure and not recommended.

Without ``--no-overwrite``, an existing file is skipped only when its size matches the remote ``Content-Length``;
otherwise it is downloaded again. ``dbCAN-PUL.tar.gz`` is extracted into ``dbCAN-PUL/`` and then deleted.
If any file fails, the command exits with an error that lists the failed files; rerun it to fetch only the missing ones.

For unreliable networks:

.. code-block:: shell

   run_dbcan database --db_dir db --aws_s3 --timeout 60 --retries 5

.. tip::

   We recommend creating a dedicated directory for databases, as they will be reused for multiple analyses.

Files and approximate sizes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sizes below are those of the pinned S3 release ``db_v5-2-9_11-3-2026`` (dbCAN HMM database V15,
CAZy DIAMOND database built from ``CAZyDB.08052026``); the HTTP snapshot may differ slightly.
Plan for roughly **8.5 GB** of disk space for a full download (about 7.8 GB with ``--no-cgc``).

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Local filename
     - Size
     - Used by
   * - ``dbCAN-sub.hmm`` (remote name ``dbCAN_sub.hmm``)
     - ~5.1 GB
     - ``CAZyme_annotation`` (dbCAN-sub)
   * - ``CAZy.dmnd``
     - ~2.6 GB
     - ``CAZyme_annotation`` (DIAMOND)
   * - ``dbCAN.hmm``
     - ~120 MB
     - ``CAZyme_annotation`` (dbCAN HMM)
   * - ``fam-substrate-mapping.tsv``
     - <1 MB
     - ``CAZyme_annotation`` (``Substrate`` column)
   * - ``peptidase_db.dmnd``
     - ~330 MB
     - ``gff_process`` (CGC)
   * - ``sulfatlas_db.dmnd``
     - ~100 MB
     - ``gff_process`` (CGC)
   * - ``TCDB.dmnd``, ``TF.hmm``, ``TF.dmnd``, ``STP.hmm``
     - <15 MB each
     - ``gff_process`` (CGC)
   * - ``PUL.dmnd``, ``dbCAN-PUL.xlsx``, ``dbCAN-PUL/``
     - <5 MB
     - ``substrate_prediction``
   * - ``Pfam-A.hmm``
     - ~1.6 GB
     - ``Pfam_null_cgc`` only. **Not downloaded** by ``run_dbcan database``; download it manually (see below).

.. note::

   ``dbCAN_sub.hmm`` is saved as ``dbCAN-sub.hmm`` (hyphen). run_dbcan looks for the hyphenated name, so keep
   it when downloading manually.

Recording the database version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For reproducible results, record which databases you used together with the software version:

.. code-block:: shell

   run_dbcan version
   ls -l db/              # file sizes and modification dates
   md5sum db/*.dmnd db/*.hmm > db/checksums.md5

If you used ``--aws_s3``, report the S3 release name (``db_v5-2-9_11-3-2026`` for the current package).
If you used the HTTP source, report the download date, because ``db_current`` changes over time.

Manual Database Installation
-----------------------------

If you prefer to download database files manually or face connectivity issues, use the same paths as the built-in
downloader. The HTTP source is:

.. code-block:: text

   https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/<filename>

The pinned S3 source is:

.. code-block:: text

   https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_11-3-2026/<filename>

Step-by-Step Instructions:

1. Create a database directory:

   .. code-block:: shell

      mkdir -p db && cd db

2. Download the required database files:

   .. code-block:: shell

      # CAZyme databases
      wget -O CAZy.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/CAZy.dmnd"
      wget -O dbCAN.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/dbCAN.hmm"
      wget -O dbCAN-sub.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/dbCAN_sub.hmm"
      wget -O fam-substrate-mapping.tsv "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/fam-substrate-mapping.tsv"

      # CGC-related databases
      wget -O TCDB.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/TCDB.dmnd"
      wget -O TF.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/TF.hmm"
      wget -O TF.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/TF.dmnd"
      wget -O STP.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/STP.hmm"
      wget -O PUL.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/PUL.dmnd"
      wget -O dbCAN-PUL.xlsx "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/dbCAN-PUL.xlsx"
      wget -O dbCAN-PUL.tar.gz "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/dbCAN-PUL.tar.gz"
      wget -O peptidase_db.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/peptidase_db.dmnd"
      wget -O sulfatlas_db.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/sulfatlas_db.dmnd"

      tar -xzf dbCAN-PUL.tar.gz

      # Optional: Pfam-A HMMs, only needed for Pfam_null_cgc
      wget -O Pfam-A.hmm "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_11-3-2026/Pfam-A.hmm"

   To download the whole pinned release with the AWS CLI instead:

   .. code-block:: shell

      aws s3 cp s3://dbcan/db_v5-2-9_11-3-2026/ ./db --recursive --no-sign-request
      mv db/dbCAN_sub.hmm db/dbCAN-sub.hmm

   The S3 folder also contains files that the downloader skips (``Pfam-A.hmm`` and the CAZy FASTA
   ``CAZyDB.<date>.fa``), so a recursive copy needs about 12 GB.

3. Verify your downloads:

   .. code-block:: shell

      ls -lh
      # You should see CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm, fam-substrate-mapping.tsv, and other files.

.. note::

   For a complete list of required database files and their descriptions,
   see :doc:`../getting_started/database_description` and the `dbCAN2 download page <https://bcb.unl.edu/dbCAN2/download/>`_.

Testing Your Installation
---------------------------

After installing the databases, you can verify your setup using example files from our `test data repository <https://bcb.unl.edu/dbCAN2/download/test>`_.

.. admonition:: Next Steps

   Once your databases are installed, proceed to the :doc:`CAZyme Annotation <CAZyme_annotation>` section to start analyzing your sequences.
