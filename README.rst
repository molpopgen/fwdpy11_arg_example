fwdpy11_arg_example
**********************************************************

This project shows what is possible in terms of tracking the ancestral recombination graph (ARG) during a foward-time simulation implemented using fwdpp_ and msprime_.  The interface to the implementation is in terms of a Python package implemented using fwdpy11_ and pybind11_.

License
----------------------------------

GPLv3 or later (see COPYING)

Overview
----------------------------------

We define a C++ class called "ancestry_tracker", which stores nodes and edges as they appear forwards in time.  These data structures, and their updating, are non-intrusive, meaning that they don't care about any of the fwdpp_ internals.  Rather, we simply have to define a new "iterate a generation" function that uses both fwdpp_ machinery and updates an ancestry_tracker as appropriate.

Using pybid11_, we make anestry_trackers visible to Python as an AncestryTracker class.  The Python class has access to the nodes and edges as NumPy structured arrays, which can be viewed "for free", meaning that no copy from C++ to Python is required to look at them.

We define an ArgSimplifier class to bridge the AncestryTracker and msprime.  The __call__ operator of ArgSimplifier will accept an AncestryTracker as an argument and use the msprime API to simplify the input data into a set of trees.


Crude usage instructions
----------------------------------

Currently, one can run simple simulations under very restrictive parameter combinations. There is currently **no** integration into msprime_ for generating the trees.  To be safe, I only recommend running the unit test for now.  While the underlying machinery is hooked up to fwdpy11_'s general scheme to model variation in mutation and recombination rates, the node/edge tracking currently assumes all positions are on the [0,1) interval.  

This has been confirmed to work in a clean conda environment using Python3.

Install the following dependencies using conda:

.. code-block:: bash

    #hdf5 needed for msprime, which we 
    #are installing from github
    conda install pybind11 gcc hdf5
    conda install -c bioconda fwdpy11

Install msprime_ from the current master branch on github. (Right now, you don't even need to do this in order to do a local build).

Make a local build and run the unit tests:

.. code-block:: bash

    #--gcc only needed on OS X.  Does no harm 
    #on Linux.
    python setup.py build_ext -i --gcc
    python -m unittest discover tests

The unit test runs N = 1e3, 4Nr = 1e4, and some sites under selection, for 10N generations.  It takes about 20 seconds on my Linux box.  

"Brute-force" simulation
+++++++++++++++++++++++++++++++++

To run a proof-of-principle example where we do an entire simulation and then have msprime clean up the mess:

.. code-block:: bash

    python test_brute_force.py N 4Nr seed

The output will be the times spent in various steps.


.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _pybind11: http://github.com/pybind/pybind11
.. _msprime: http://github.com/jeromekelleher/msprime
