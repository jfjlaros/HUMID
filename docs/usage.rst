Usage
=====

Word length
-----------
By default, the word length is 24 nucleotides, this can be changed with the
``-n`` flag. The number of nucleotides will be equally divided over all input
files. If this is not possible, the remainder is taken from the **last** input
file. For example, if three input files are specified, and ``-n`` is set to 23, 7
nucleotides will be taken from the first and second input file, and 9 from the
last input file.

If the UMI is present in the header, all nucleotides from the UMI will be used,
and the remainder will be divided between the input files as described.

UMIs in separate file
---------------------
As this tool can handle an arbitrary number of input files, simply specify the
file with the UMIs as an input file. Note that with a separate UMI file, there
will be no preference to include all UMI nucleotides.

::

    dedup forward.fastq.gz reverse.fastq.gz umi.fastq.gz


UMIs in the header with underscore
-----------------------------------
If UMIs are include in the header, the file with the UMI in the header must be
the **first** specified input file, and the UMI must be before any spaces.

Below are some examples of headers with the UMI ``AGTA`` after an
underscore, which will all be recognised.

::

    @A31886:289:T5D5W10Y2:2:12686:4678:1110_AGTA
    @1_AGTA
    @1_AGTA with spaces after


This format of specifying UMIs in the read header is used by UMI-Tools.

UMIs in the header with colon
-----------------------------
If UMIs are include in the header, the file with the UMI in the header must be
the **first** specified input file, and the UMI must be before any spaces.

Below are some examples of headers with the UMI ``AGTA`` that will be
recognised.

::

    @A31886:289:T5D5W10Y2:2:12686:4678:1110:AGTA
    @A31886:289:T5D5W10Y2:2:12686:4678:1110:AGTA with spaces after
    @1:::::::AGTA


This format of specifying UMIs in the read header is used by ``BCL Convert``
and ``fastp``.

Other UMI formats
-----------------
If your data does not confirm to any of the supported schemas described above,
we recommend fastp_ to process the UMI to one of the supported formats.


Deduplication without UMI
-------------------------
If a project was sequenced without UMIs, you can still remove duplicates using
this tool, since it will automatically take the requested number of nucleotides
from one (single end) or two (paired end) FastQ files. Note that without random
UMIs to distinguish identical but independent molecules, the number of
duplicates will most likely be an overestimation, similar to using picard
MarkDuplicates.

.. _fastp: https://github.com/OpenGene/fastp#unique-molecular-identifier-umi-processing
