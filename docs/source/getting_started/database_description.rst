Database Description
====================

``run_dbcan database --db_dir db`` downloads the files below (except the two optional ones) into one folder (see
:doc:`../user_guide/prepare_the_database` for sizes, download sources and versions). The databases are usually
updated once a year, between July and September.

.. list-table::
   :widths: 27 43 30
   :header-rows: 1

   * - File
     - Content
     - Used by
   * - ``CAZy.dmnd``
     - CAZy proteins (DIAMOND)
     - ``CAZyme_annotation``
   * - ``dbCAN.hmm``
     - dbCAN family HMMs
     - ``CAZyme_annotation``
   * - ``dbCAN-sub.hmm``
     - dbCAN-sub subfamily HMMs (remote name ``dbCAN_sub.hmm``)
     - ``CAZyme_annotation``
   * - ``fam-substrate-mapping.tsv``
     - CAZy family/EC to substrate table; fills the ``Substrate`` column of ``overview.tsv``
     - ``CAZyme_annotation``
   * - ``dbCAN-sub-substrate.tsv`` (optional)
     - Per-subfamily (Subfamily, EC, Count, Substrate) table of the dbCAN-sub web server; used before
       ``fam-substrate-mapping.tsv`` when present; not downloaded, copy it into ``--db_dir``
     - ``CAZyme_annotation``
   * - ``TCDB.dmnd``
     - Transporters (TCDB)
     - ``gff_process`` (TC)
   * - ``TF.dmnd``
     - Prokaryotic transcription factors (PRODORIC)
     - ``gff_process`` (TF, default)
   * - ``TF.hmm``
     - Fungal transcription factor HMMs (MycoCosm)
     - ``gff_process --fungi``
   * - ``STP.hmm``
     - Signal transduction protein HMMs
     - ``gff_process`` (STP)
   * - ``sulfatlas_db.dmnd``
     - Sulfatases (SulfAtlas)
     - ``gff_process`` (Sulfatase)
   * - ``peptidase_db.dmnd``
     - MEROPS peptidase units (see :ref:`rebuild-peptidase-db`)
     - ``gff_process`` (Peptidase)
   * - ``PUL.dmnd``, ``dbCAN-PUL.xlsx``, ``dbCAN-PUL/``
     - dbCAN-PUL proteins, their substrates and their gene maps (synteny plots)
     - ``substrate_prediction``
   * - ``Pfam-A.hmm`` (optional)
     - Pfam domains; **not** downloaded by ``run_dbcan database``
     - ``Pfam_null_cgc``

``--no-cgc`` downloads only the first four files (CAZyme annotation).
