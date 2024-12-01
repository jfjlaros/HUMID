Troubleshooting
===============

If you believe you have encountered a bug, please open an issue on github_.

Segmentatio fault during cluster calculation
--------------------------------------------
If you see a message similar to the one below, it is possible that HUMID has
run out of stack space::

  Reading data... done. (22m42s)
  Calculating neighbours using Hamming distance... done. (17m10s)
  Calculating maximum clusters... Segmentation fault (core dumped)

To resolve this, you can increase the maximum stack size on your system using
``ulimit -Ss 102400``.

.. _github: https://github.com/jfjlaros/HUMID/issues
