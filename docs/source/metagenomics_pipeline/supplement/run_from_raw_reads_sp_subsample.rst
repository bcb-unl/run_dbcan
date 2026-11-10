Supplementary protocol: sub-sampling (Carter 2023)
==================================================

Samples sequenced to different depths assemble differently. To check that a result does not come from unequal
depth, draw the same number of read pairs (20 million) from each sample and repeat the protocol. Sub-sampling
also lowers CPU time and RAM, which helps with co-assembly of many samples.

The steps are those of :ref:`metagenomics-protocol` (same data, ``SAMPLES="Wet2014 Dry2014"``); only P2 changes.
Add ``bbmap`` (BBTools, for ``reformat.sh``) to the conda environment.

**P2. Trim, then sub-sample** (~36 min). Run :ref:`P2 <mg-p2>`, then replace the trimmed reads by a random sample
of 20 M pairs, so that P3-P15 run unchanged:

.. code-block:: shell

   for s in $SAMPLES; do
     mv ${s}_1_val_1.fq.gz ${s}_1_val_1.full.fq.gz
     mv ${s}_2_val_2.fq.gz ${s}_2_val_2.full.fq.gz
     reformat.sh in=${s}_1_val_1.full.fq.gz in2=${s}_2_val_2.full.fq.gz \
         out=${s}_1_val_1.fq.gz out2=${s}_2_val_2.fq.gz sample=20000000
   done

Then run :ref:`P3-P15 <mg-p3>`. With fewer reads, MEGAHIT (P3) takes ~50 min instead of ~4 h 20 min.

.. note::
   The read mapping (P8) can use either the sub-sampled or the full reads; the full reads give more counts per
   gene. Precomputed results: https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/.
