fwdpy11_arg_example
**********************************************************

This project shows what is possible in terms of tracking the ancestral recombination graph (ARG) during a foward-time simulation implemented using fwdpp_ and msprime_.  The interface to the implementation is in terms of a Python package implemented using fwdpy11_ and pybind11_.

License
----------------------------------

GPLv3 or later (see COPYING)

Overview
----------------------------------

We define a C++ class called "ancestry_tracker", which stores nodes and edges as they appear forwards in time.  These data structures, and their updating, are non-intrusive, meaning that they don't care about any of the fwdpp_ internals.  Rather, we simply have to define a new "iterate a generation" function that uses both fwdpp_ machinery and updates an ancestry_tracker as appropriate.

Using pybind11_, we make ancestry_trackers visible to Python as an AncestryTracker class.  The Python class has access to the nodes and edges as NumPy structured arrays, which can be viewed "for free", meaning that no copy from C++ to Python is required to look at them.

We define an ArgSimplifier class to bridge the AncestryTracker and msprime.  The __call__ operator of ArgSimplifier will accept an AncestryTracker as an argument and use the msprime API to simplify the input data into a set of trees.

Caution
----------------------------------

This is a proof of principle implementation only.  It has been tested that it gives the correct sample properties in distribution and that internal data structures are sane during simulation.  The API is not guaranteed to be neither ideal nor idiomatic.  Further, various safety checks are missing on the C++ side.  We simply do not check for integer overflow or other possible things that may happen if garbage collection occurs too infrequently.  These limitations, and others, will be dealt with (and tested for) when we port these ideas into fwdpy11_.

Crude usage instructions
----------------------------------

Currently, one can run simple simulations under very restrictive parameter combinations. There is currently **no** integration into msprime_ for generating the trees.  To be safe, I only recommend running the unit test for now.  While the underlying machinery is hooked up to fwdpy11_'s general scheme to model variation in mutation and recombination rates, the node/edge tracking currently assumes all positions are on the [0,1) interval.  

This has been confirmed to work in a clean conda environment using Python3.  **We strongly recommend that this package is installed into a clean conda env.**

Instructions for conda on Linux:

.. code-block:: bash

    conda create --name ftprime_ms python
    source activate ftprime_ms
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install gcc fwdpy12==0.1.4 msprime==0.5.0 pybind11==2.2.1 pandas
    https://github.com/molpopgen/fwdpy11_arg_example
    cd fwdpy11_arg_example
    python setup.py build_ext -i

Test simulation
+++++++++++++++++++++++++++++++++

To run a proof-of-principle example where we do an entire simulation and then have msprime clean up the mess:

.. code-block:: bash

    python test_evolve.py N 4Nr seed

The output will be the times spent in various steps.

Running the simulations found in the paper
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The script `benchmarking.py` was used for running the simulations in the Kelleher et al. manuscript.

Source code overview
-----------------------------------------

The package consists of a mix of C++ and Python code. All source code is in the fwdpy11_arg_example subdirectory of thie main repository.

For example, let's run a simulation with the following parameters:

* `N=5e4` diploids 
* a region size of `theta = rho = 1e4` 
* deleterious mutation rate equal to one percent of the neutral mutation rate
* simplify ("GC", or garbage-collect) every 1,000 generations
* apply mutations to a sample of 100 individuals at the end of the simulation
* write the timings to a file called `timings.txt`

The command line for the above is:

.. code-block:: bash

    python benchmarking.py --popsize 50000 --theta 10000 --rho 10000 --pdel 0.01 --gc 1000 --nsam 100 --outfile1 timings.txt.gz \
    --seed $RANDOM

.. note:: The output file is gzip compressed!

Please be mindful of running these simulations on machines with little RAM!  In general, forward simulations are
intended to be run on HPC-strength hardware.  While tree sequence simplification results in very efficient run times, we
are sometimes still using a substantial amount of RAM.

An example of the output is:

.. code-block:: bash

    prepping	sorting	appending	simplifying	fwd_sim_runtime	N	theta	rho	simplify_interval
    0.05370585599999733	4.384206619999995	0.2173980950000004	2.9446647440000016	5.174604999999977	1000	1000.0	1000.0	100

The fields are:

* `prepping`: cumulative time spent preparing data for a copy from the C++ side to msprime
* `sorting`: cumulative time spent sorting tables, which is a requirement for simplification
* `simplifying`: cumulative time spent simplifying tables
* `fwd_sim_runtime`: The total time spent simulating

The remaining four columns are the command-line parameters.

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
