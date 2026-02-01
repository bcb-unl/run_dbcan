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

.. tip::

   We recommend creating a dedicated directory for databases, as they will be reused for multiple analyses.

Manual Database Installation
-----------------------------

If you prefer to download database files manually or face connectivity issues, you can obtain all required files directly from the `dbCAN2 website <http://bcb.unl.edu/dbCAN2/download/Databases/>`_.

Step-by-Step Instructions:

1. Create a database directory:

   .. code-block:: shell

      mkdir -p db && cd db

2. Download the required database files:

   .. code-block:: shell

      # Download CAZy database
      wget https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/CAZy.dmnd

      # Download dbCAN HMM database
      wget https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN_sub.hmm

      #and download other files as needed

3. Verify your downloads:

   .. code-block:: shell

      ls -lh
      # You should see CAZy.dmnd, dbCAN_sub.hmm, and other files

.. note::

   For a complete list of required database files and their descriptions,
   please visit the `dbCAN2 database documentation <https://bcb.unl.edu/dbCAN2/download/>`_.

Testing Your Installation
---------------------------

After installing the databases, you can verify your setup using example files from our `test data repository <https://bcb.unl.edu/dbCAN2/download/test>`_.

.. admonition:: Next Steps

   Once your databases are installed, proceed to the :doc:`CAZyme Annotation <CAZyme_annotation>` section to start analyzing your sequences.
