Database Description
====================

All files are available on the `Web server <https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/>`_.
An AWS-hosted dataset is planned and will be released soon.
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
  - Filename: ``dbCAN_sub.hmm``


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


Substrate Prediction Databases
------------------------------

- Substrate mapping table
  - Description: Mapping from CAZyme family/EC to known substrates.
  - Filename: ``fam-substrate-mapping.tsv``

- dbCAN-PUL substrate table
  - Description: Substrate mapping associated with PULs from the dbCAN-PUL database.
  - Directory: ``dbCAN-PUL/``
  - Spreadsheet: ``dbCAN-PUL.xlsx``
