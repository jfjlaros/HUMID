Installation
============

You can download a static binary directly from GitHub_, or use any of the
alternative installation methods listed below


From Conda_
-----------

::

    conda install humid

From source
-----------

::

    git clone https://github.com/jfjlaros/HUMID.git
    cd HUMID/src
    git submodule update --recursive --init
    make deps

Compilation
~~~~~~~~~~~

::

    make

Static compilation
~~~~~~~~~~~~~~~~~~

::

    make static

.. _Conda: https://anaconda.org/bioconda/humid
.. _GitHub: https://github.com/jfjlaros/HUMID/releases
