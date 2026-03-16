Output
======

FastQ output
------------
By default, HUMID will write the deduplicated FastQ files in the current
folder, using the `_dedup` suffix in the file name to distinguish them from the
input FastQ files.

By specifying the `-a` flag, HUMID will output the annotated FastQ files using
the `_annotated` suffix in the file name. For each read in the output, the
`cluster_id` will be appended to the end of the read header, using a colon
(`:`) as a separator.

**Special case**: The cluster with id `0` has been reserved for reads that
could not be classified. For example because there were not enough bases
available to create a `word`, or because the word contains one or more N bases.


Statistics
----------
Run HUMID with the `-s` flag to generate deduplication statistics. These
statics files can be visualized using MultiQC version 1.14 or later, or
inspected directly.

stats.dat
~~~~~~~~~
This is probably the most usefull file, it contains the statistics about the number of reads in each of the following categories.

.. list-table:: stats.dat
  :header-rows: 1

  * - Field
    - Definition
  * - total
    - Total number of input reads
  * - usable
    - Total number of reads that were usable (did not contain N)
  * - unique
    - Total number of distinct input
  * - clusters
    - Total number unique reads after clustering and deduplication

neigh.dat
~~~~~~~~~
This file contains counts for the number of neighbours for each usable read. The first number is the `number of neighbours` and the second number is how many distinct reads have this number of neighbours.

clusters.dat
~~~~~~~~~~~~
This file contains the counts for the size of each cluster. The first number is the `cluster size`, the second number is how many clusters are of the specified size.

counts.dat
~~~~~~~~~~
This file contains the counts for exact duplicates in the usable input reads. The first number is the `number of exact duplicates` and the second number is how many distinct reads have this number of exact duplicates.

