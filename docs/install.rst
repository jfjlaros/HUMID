Installation
============

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
