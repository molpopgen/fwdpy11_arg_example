#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <fwdpy11/types.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpp/diploid.hh>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include "ancestry_edge_sets.hpp"
#include "evolve_generation.hpp"
namespace py = pybind11;

PYBIND11_PLUGIN(wfarg)
{
    py::module m("wfarg", "Simple example of Wright-Fisher simulation with "
                          "selection and ARG tracking");

    //Register nodes and edges as NumPy dtypes:
    PYBIND11_NUMPY_DTYPE(node, id, generation, deme);
    PYBIND11_NUMPY_DTYPE(edge, left, right, parent, child);

    //Register vectors of nodes and edges as "opaque"
    //types supporting Python's buffer protocol, creating
    //Python classes that are castable to NumPy structured
    //arrays without a copy.
    py::bind_vector<std::vector<node>>(
        m, "NodeArray", "Container of nodes. This can be cast to a NumPy "
                        "record array without making a copy.",
        py::buffer_protocol());

    py::bind_vector<std::vector<edge>>(
        m, "EdgeArray", "Container of edges.  This can be cast to a NumPy "
                        "record array without making a copy",
        py::buffer_protocol());

    //Make our C++ function callable from Python.
    //This is NOT part of a user-facing Python API.
    //Rather, we need a wrapper to integrate it with
    //the rest of the fwdpy11 world.
    m.def("evolve_singlepop_regions_track_ancestry",
          &evolve_singlepop_regions_track_ancestry);

    return m.ptr();
}
