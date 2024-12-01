=========================================
HUMID: reference free FastQ deduplication
=========================================

.. image:: https://img.shields.io/github/last-commit/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID/graphs/commit-activity
.. image:: https://github.com/jfjlaros/HUMID/actions/workflows/cpp-library.yml/badge.svg
   :target: https://github.com/jfjlaros/HUMID/actions/workflows/cpp-library.yml
.. image:: https://readthedocs.org/projects/HUMID/badge/?version=latest
   :target: https://HUMID.readthedocs.io/en/latest
.. image:: https://img.shields.io/github/release-date/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID/releases
.. image:: https://img.shields.io/github/release/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID/releases
.. image:: https://img.shields.io/github/languages/code-size/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID
.. image:: https://img.shields.io/github/languages/count/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID
.. image:: https://img.shields.io/github/languages/top/jfjlaros/HUMID.svg
   :target: https://github.com/jfjlaros/HUMID
.. image:: https://img.shields.io/github/license/jfjlaros/HUMID.svg
   :target: https://raw.githubusercontent.com/jfjlaros/HUMID/master/LICENSE.md

----

HUMID is a tool to quickly and easily remove duplicate reads from FastQ files, with or without UMIs.


Installation
============

You can install HUMID from conda

.. code-block:: bash

    conda install -c bioconda humid

If you want to, you can also install HUMID `from source <https://humid.readthedocs.io/en/latest/install.html#from-source>`_.


Usage
=====

Both the input and output of HUMID are plain FastQ files, no alignment
required! This means that you can use HUMID to remove duplicates as a
pre-processing step before starting your analysis. If your project was
sequenced without UMIs, or if the UMIs are present in the headers of the FastQ
reads (as is done by BCL Convert), you can use the following command:

.. code-block:: bash

    humid forward.fastq.gz reverse.fastq.gz


If the UMIs are located in a separate FastQ file use

.. code-block:: bash

    humid forward.fastq.gz reverse.fastq.gz umi.fast.gz


For other use cases, we recommend that you use `fastp
<https://github.com/OpenGene/fastp#unique-molecular-identifier-umi-processing>`_
to move the UMIs to the header of the forward FastQ file before deduplicating
them with HUMID.

Please see the `usage <https://humid.readthedocs.io/en/latest/usage.html>`_
section of the documentation for more details.
