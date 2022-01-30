Installation
============

From source
-----------

::

    git clone https://github.com/jfjlaros/dedup.git
    cd dedup
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
    ./autogen
    ./configure
    make
    cd ../src
    make static
