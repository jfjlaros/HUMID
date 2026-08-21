Tutorial
========

This tutorial assume you have already installed HUMID and MultiQC.


Download FASTQ data
-------------------

First, we need to download some FASTQ files to analyze with HUMID. Here, we
will use data from the Genome in a Bottle project, specifically the son of the
Ashkenazim trio (HG002). See this GIAB_ repository for more information.

::

  wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
  wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz


Optionally, you can verify the md5 checksum to ensure the data was not
corrupted during the download

::

 md5sum *gz

 c2ae5e412fb211974f9a9a46a5392428  2A1_CGATGT_L001_R1_001.fastq.gz
 83826a956fc90c501645391314b2abf3  2A1_CGATGT_L001_R2_001.fastq.gz


Run HUMID
---------
Next, we will analyze the FASTQ files with HUMID. Note that these files do not
contain an UMI, but HUMID will detect this and use the sequences from the reads
themselves to detect duplicates.

The following command will run HUMID and write the deduplicated FASTQ files and
statistics into the GIAB folder. It should not take long to analyze the 4
million read pairs from the FASTQ files.

::

  humid -d GIAB -s 2A1_CGATGT_L001_R1_001.fastq.gz 2A1_CGATGT_L001_R2_001.fastq.gz

  Reading data... done. (0m9s)
  Calculating neighbours using Hamming distance... done. (0m26s)
  Calculating directional clusters... done. (0m3s)
  Writing filtered results... done. (0m31s)
  Calculating count and neighbour stats... done. (0m4s)


Visualize the statistics with MultiQC
-------------------------------------
If you have a recent (>=1.14) version of MultiQC installed you can use it to
visualize the duplication statistics

::

  multiqc GIAB

If you inspect the produced multiqc_report.html file you will see that around
99.3% of the 4 million analyzed reads were identified as unique reads. In the
detailed section under the HUMID heading you can see that HUMID identified
around 26,000 duplicate reads, while almost 4,000 reads were filtered out
because they contained uncertain characters (N) or were too short.

.. _GIAB: https://github.com/genome-in-a-bottle/giab_data_indexes
