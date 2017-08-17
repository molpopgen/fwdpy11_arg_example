fwdpy11_arg_example
**********************************************************

This project shows what is possible in terms of tracking the ancestral recombination graph (ARG) during a foward-time simulation implemented using fwdpp_ and msprime_.  The interface to the implementation is in terms of a Python package implemented using fwdpy11_ and pybind11_.

License
----------------------------------

GPLv3 or later (see COPYING)

Crude usage instructions
----------------------------------

Currently, one can run simple simulations under very restrictive parameter combinations. There is currently **no** integration into msprime_ for generating the trees.  To be safe, I only recommend running the unit test for now.  While the underlying machinery is hooked up to fwdpy11_'s general scheme to model variation in mutation and recombination rates, the node/edge tracking currently assumes all positions are on the [0,1) interval.  

This has been confirmed to work in a clean conda environment using Python3.

Install the following dependencies using conda:

.. code-block:: bash

    conda install pybind11 gcc
    conda install -c bioconda fwdpy11

Install msprime_ from the current master branch on github. (Right now, you don't even need to do this in order to do a local build).

Make a local build and run the unit tests:

.. code-block:: bash

    #--gcc only needed on OS X.  Does no harm 
    #on Linux.
    python setup.py build_ext -i --gcc
    python -m unittest discover tests

The unit test runs N = 1e3, 4Nr = 1e4, and some sites under selection, for 10N generations.  It takes about 20 seconds on my Linux box.  

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _pybind11: http://github.com/pybind/pybind11
.. _msprime: http://github.com/jeromekelleher/msprime
