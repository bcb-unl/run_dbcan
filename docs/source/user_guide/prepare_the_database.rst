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
The default HTTP source points at the moving ``db_current`` snapshot. The ``--aws_s3`` source points at the pinned
``db_v5-2_9-13-2025`` release, which is preferable when you need reproducible database contents across machines.

For unreliable networks, the downloader also supports retry and resume controls:

.. code-block:: shell

   run_dbcan database --db_dir db --aws_s3 --timeout 60 --retries 5 --resume

Use ``--no-overwrite`` to keep existing files untouched. By default, existing files may be checked by remote size and
completed partial downloads are resumed when possible.

.. tip::

   We recommend creating a dedicated directory for databases, as they will be reused for multiple analyses.

Manual Database Installation
-----------------------------

If you prefer to download database files manually or face connectivity issues, use the same paths as the built-in
downloader. The HTTP source is:

.. code-block:: text

   https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_current/<filename>

The pinned S3 source is:

.. code-block:: text

   https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2_9-13-2025/<filename>

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

3. Verify your downloads:

   .. code-block:: shell

      ls -lh
      # You should see CAZy.dmnd, dbCAN.hmm, dbCAN-sub.hmm, fam-substrate-mapping.tsv, and other files.

.. note::

   For a complete list of required database files and their descriptions,
   please visit the `dbCAN2 database documentation <https://bcb.unl.edu/dbCAN2/download/>`_.

Testing Your Installation
---------------------------

After installing the databases, you can verify your setup using example files from our `test data repository <https://bcb.unl.edu/dbCAN2/download/test>`_.

.. admonition:: Next Steps

   Once your databases are installed, proceed to the :doc:`CAZyme Annotation <CAZyme_annotation>` section to start analyzing your sequences.
