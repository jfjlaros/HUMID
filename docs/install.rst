Installation
============

You can download a static binary directly from _github, or use any of the
alternative installation methods listed below

From conda_
-----------

::

    conda install humid

From source
-----------

::

    git clone https://github.com/jfjlaros/HUMID.git
    cd HUMID
    git submodule update --recursive --init

Compilation
~~~~~~~~~~~

::

    cd src
    make

Static compilation
~~~~~~~~~~~~~~~~~~

::

    cd lib/isa-l
    ./autogen.sh
    ./configure
    make
    cd ../../src
    make static

.. _conda: https://anaconda.org/bioconda/humid
.. _github: https://github.com/jfjlaros/HUMID/releases
