HUMID
^^^^^
----

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

Quickly and easily remove duplicate reads from FastQ files, with or without UMIs.


Installation
------------
You can install HUMID from conda

.. code-block:: bash

    conda install -c bioconda humid

If you want to, you can also install HUMID `from source <https://humid.readthedocs.io/en/latest/install.html#from-source>`_


Usage
-----
Both the input and output of HUMID are plain FastQ files, so it you can simply
remove duplicates as a pre-processing step before starting your analysis:

.. code-block:: bash

    humid forward.fastq.gz reverse.fastq.gz

Please see the `usage <https://humid.readthedocs.io/en/latest/usage.html>`_ section of the documentation for details.
