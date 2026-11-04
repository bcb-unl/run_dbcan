Database Description
====================

The built-in downloader uses the ``db_current`` HTTP snapshot from ``pro.unl.edu`` by default. This snapshot moves
when the databases are updated. Use the ``--aws_s3`` flag with the ``run_dbcan database`` command to download the
pinned S3 release instead, which is faster in many regions and reproducible
(see :doc:`Preparing Databases <../user_guide/prepare_the_database>` for file sizes and how to record the database version).
The databases are generally updated annually between July and September.

The databases used by run_dbCAN are described below.

CAZyme Databases
----------------

- DIAMOND database

  - Description: Fast protein alignment database used to identify CAZyme sequences.
  - Filename: ``CAZy.dmnd``

- dbCAN HMM database

  - Description: HMM profiles for CAZyme families; sensitive identification.
  - Filename: ``dbCAN.hmm``

- dbCAN-sub HMM database

  - Description: Subfamily-level HMM profiles for fine-grained CAZyme subfamily identification.
  - Filename: ``dbCAN-sub.hmm`` (published remotely as ``dbCAN_sub.hmm``)


CGC Databases
-------------

- Transporter DIAMOND database (from TCDB)

  - Description: Transporter proteins used to identify transporters in CGCs.
  - Filename: ``TCDB.dmnd``

- Transcription Factor HMM database (fungi, from MycoCosm)

  - Description: HMM profiles for transcription factors used in fungal datasets.
  - Filename: ``TF.hmm``

- Transcription Factor DIAMOND database (prokaryotes, from PRODORIC)

  - Description: TF protein database used to identify TFs in prokaryotic CGCs.
  - Filename: ``TF.dmnd``

- Signal Transduction Protein HMM database

  - Description: HMM profiles for signal transduction proteins.
  - Filename: ``STP.hmm``

- Sulfatase DIAMOND database (from SulfAtlas)

  - Description: Sulfatase protein database used to identify sulfatases.
  - Filename: ``sulfatlas_db.dmnd``

- Peptidase DIAMOND database (from MEROPS)

  - Description: Peptidase protein database used to identify peptidases.
  - Filename: ``peptidase_db.dmnd``

- dbCAN-PUL DIAMOND database

  - Description: Polysaccharide Utilization Loci (PUL) protein database for PUL identification.
  - Filename: ``PUL.dmnd``

- Pfam-A HMM database (optional)

  - Description: Pfam domain HMMs used by ``Pfam_null_cgc`` to annotate null genes inside CGCs.
  - Filename: ``Pfam-A.hmm``
  - Not downloaded by ``run_dbcan database``; see :doc:`../user_guide/prepare_the_database`.


Substrate Prediction Databases
------------------------------

- Substrate mapping table

  - Description: Mapping from CAZyme family/EC to known substrates. Used during ``CAZyme_annotation`` to fill the
    ``Substrate`` column of ``overview.tsv`` and ``dbCANsub_hmm_results.tsv``.
  - Filename: ``fam-substrate-mapping.tsv``

- dbCAN-PUL substrate table

  - Description: Substrate mapping associated with PULs from the dbCAN-PUL database.
  - Directory: ``dbCAN-PUL/``
  - Spreadsheet: ``dbCAN-PUL.xlsx``
