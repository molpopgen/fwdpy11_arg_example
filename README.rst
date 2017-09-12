fwdpy11_arg_example
**********************************************************

This project shows what is possible in terms of tracking the ancestral recombination graph (ARG) during a foward-time simulation implemented using fwdpp_ and msprime_.  The interface to the implementation is in terms of a Python package implemented using fwdpy11_ and pybind11_.

License
----------------------------------

GPLv3 or later (see COPYING)

Overview
----------------------------------

We define a C++ class called "ancestry_tracker", which stores nodes and edges as they appear forwards in time.  These data structures, and their updating, are non-intrusive, meaning that they don't care about any of the fwdpp_ internals.  Rather, we simply have to define a new "iterate a generation" function that uses both fwdpp_ machinery and updates an ancestry_tracker as appropriate.

Using pybind11_, we make anestry_trackers visible to Python as an AncestryTracker class.  The Python class has access to the nodes and edges as NumPy structured arrays, which can be viewed "for free", meaning that no copy from C++ to Python is required to look at them.

We define an ArgSimplifier class to bridge the AncestryTracker and msprime.  The __call__ operator of ArgSimplifier will accept an AncestryTracker as an argument and use the msprime API to simplify the input data into a set of trees.


Crude usage instructions
----------------------------------

Currently, one can run simple simulations under very restrictive parameter combinations. There is currently **no** integration into msprime_ for generating the trees.  To be safe, I only recommend running the unit test for now.  While the underlying machinery is hooked up to fwdpy11_'s general scheme to model variation in mutation and recombination rates, the node/edge tracking currently assumes all positions are on the [0,1) interval.  

This has been confirmed to work in a clean conda environment using Python3.  **We strongly recommend that this package is installed into a clean conda env.**

Install the following dependencies using conda:

.. code-block:: bash

    #hdf5 needed for msprime, which we 
    #are installing from github
    conda install -c conda-forge pybind11==2.1.1 gcc hdf5
    conda install -c bioconda fwdpy11

.. note::
    This code was developed using pybind11 version 2.1.1.  Version 2.2.0 of that project changes how C++ containers are made "transparent" to Python.  Please make sure you are using 2.1.1!!!

Install msprime_ from the current master branch on github. 

Make a local build and run the unit tests:

.. code-block:: bash

    #--gcc only needed on OS X.  Does no harm 
    #on Linux.
    python setup.py build_ext -i --gcc
    python -m unittest discover tests


Test simulation
+++++++++++++++++++++++++++++++++

To run a proof-of-principle example where we do an entire simulation and then have msprime clean up the mess:

.. code-block:: bash

    python test_evolve.py N 4Nr seed

The output will be the times spent in various steps.

Source code overview
-----------------------------------------

The package consists of a mix of C++ and Python code. All source code is in the fwdpy11_arg_example subdirectory of thie main repository.

C++ code
+++++++++++++++++++++

We define nodes and edges as simple structs, meaning that they are "C-like", consisting only of POD and no constructors or other C++ stuff.  This simple design allows C++ vectors of these structs to be treated as NumPy record arrays visible fom Python without needing to make a copy.

* `node.hpp` defines a node as a simple C-like struct.
* `edge.hpp` defines and edge as a simple C-like struct.
* `ancestry_tracker.hpp` defines a C++ struct/class called ancestry_tracker to accumulate nodes and edges during a simulation.
* `evolve_generation.hpp` handles the details of updating a Wright-Fisher population with an ancestry_tracker.
* `handle_recombination.cc/.hpp` handles the conversion of fwdpp's recombination breakpoints into types use to make edges.
* `wfarg.cc` defines a Python module (called `wfarg`) implemented in C++ via pybind11_.  It exposes our C++ back-end to Python.  The most important user-facing type defined is AncestryTracker, which wraps the C++ ancestry_tracker.

Python code
+++++++++++++++++++++

* `argsimplifier.py` defines `ArgSimplifier`, which is the bridge between the C++ code to evolve a population and the msprime_ functionality to simplify the simulated nodes and edges.
* `evolve_arg.py` defines a function that evolves a population while tracking its ancestry.  It integrates concepts from fwdpy11_ with the types defined in this package.

.. _fwdpy11: http://molpopgen.github.io/fwdpy11
.. _fwdpp: http://molpopgen.github.io/fwdpp
.. _pybind11: http://github.com/pybind/pybind11
.. _msprime: http://github.com/jeromekelleher/msprime
