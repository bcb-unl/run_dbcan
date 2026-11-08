CAZyme annotation comparison
==============================

This section compares the performance and results of run_dbCAN v5 with previous versions, highlighting improvements in accuracy, speed, and output formatting.

Annotation Results Comparison
-------------------------------

We compared the annotation results between v4 and v5 using identical input datasets (E. coli K-12 MG1655 proteome). The comparison shows that core CAZyme predictions **remain consistent** between versions, confirming that the accuracy has been maintained while making significant improvements to code structure and performance.

.. figure:: ../_static/img/overview_v4.png
   :width: 100%
   :alt: Comparison of annotation results between v4 and v5

   Figure 1: Comparison of annotation results:run_dbCAN v4


.. figure:: ../_static/img/overview_v5.png
   :width: 100%
   :alt: Comparison of annotation results between v4 and v5

   Figure 2: Comparison of annotation results: run_dbCAN v5

Key Improvements in Output Format
-----------------------------------

The new version (v5) provides several improvements in the output formatting:

1. **More Precise Domain Boundaries**

   The v5 output now includes precise domain boundary information for both dbCAN-HMM and dbCAN-sub HMM:

   .. code-block:: text

       # v5 format with precise domain boundaries
       NP_414747.1  -|-|-  GH23(101-244)  GH23_e819(102-244)+CBM50_e338(344-384)+CBM50_e338(403-442)  CBM50+GH23  3  GH23_e819|CBM50_e338|CBM50_e338

       # v4 format with limited domain information
       NP_414747.1  -|-|-  GH23(101-244)  GH23_e819+CBM50_e338+CBM50_e338  CBM50+GH23  3

2. **New "Recommend Results" Column**

   The v5 version adds a new column showing the recommended annotation results, making it easier for users to interpret findings.
   Now we follow the rule: `CAZy-sub in dbCAN-HMM > dbCAN-subfam in dbCAN-sub-HMM > dbCAN-fam in dbCAN-HMM` for the final results:

   .. code-block:: text

       Gene ID  EC#  dbCAN_hmm  dbCAN_sub  DIAMOND  #ofTools  Recommend Results
       NP_414632.1  2.4.1.227:11  GT28(185-341)  GT28_e46(185-341)  GT28  3  GT28_e46

   Current versions also append a ``Substrate`` column. All columns are described in :ref:`overview-columns`.

3. **Cleaner DIAMOND Results**

   The v5 version eliminates extraneous file paths from DIAMOND results, providing cleaner output:

   .. code-block:: text

       # v5 format with clean DIAMOND results
       NP_414555.1  -  -  -  GT1  1  -

       # v4 format with file paths in results
       NP_414555.1  -  -  -  Melli1_GeneCatalog_proteins_20150227.aa.fasta+GT1  1

Performance Comparison
------------------------

The new version shows significant performance improvements due to the implementation of pyHMMER and pyrodigal (tested on 40 CPUs):

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Dataset
     - V4 Runtime
     - V5 Runtime
   * - E. coli genome (~4.5 Mb)
     - 32 min 24 sec
     - 5 min 58 sec
   * - Xylona heveae genome (~4.5 Mb)
     - 1 hr 02 min 02 sec
     - 5 min 50 sec

.. hint::
      Since the IO and computing capabilities of different server CPUs are different, this data is for reference only.


