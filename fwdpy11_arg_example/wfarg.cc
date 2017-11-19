#include <future>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "ancestry_tracker.hpp"
#include "evolve_generation.hpp"
#include "evolve_functions.hpp"

namespace py = pybind11;


//Register vectors of nodes and edges as "opaque"
PYBIND11_MAKE_OPAQUE(std::vector<node>);
PYBIND11_MAKE_OPAQUE(std::vector<edge>);
PYBIND11_MAKE_OPAQUE(std::vector<ancestry_tracker::integer_type>);

PYBIND11_MODULE(wfarg, m)
{
    m.doc() = "Simple example of Wright-Fisher simulation with selection and "
              "ARG tracking";

    //Register nodes and edges as NumPy dtypes:
    PYBIND11_NUMPY_DTYPE(node, id, population, generation);
    PYBIND11_NUMPY_DTYPE(edge, left, right, parent, child);

    //Create Python classes of node/edge containers.
    //These types support Python's buffer protocol, creating
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

    py::bind_vector<std::vector<ancestry_tracker::integer_type>>(
        m, "VecInt32", "Vector of 32-bit, signed integers.  Castable to Numpy "
                       "array without copy.",
        py::buffer_protocol());

    //Expose the C++ ancestry_tracker to Python.
    //We only expose the stuff that a user really needs
    //to see.
    py::class_<ancestry_tracker>(m, "AncestryTracker")
        .def(py::init<decltype(edge::parent), bool, decltype(edge::parent)>(),
             py::arg("N"), py::arg("init_with_TreeSequence"),
             py::arg("next_index"))
        .def_readwrite("nodes", &ancestry_tracker::nodes,
                       "Data for msprime.NodeTable.")
        .def_readwrite("edges", &ancestry_tracker::edges,
                       "Data for msprime.EdgesetTable.")
        .def_readwrite("samples", &ancestry_tracker::offspring_indexes,
                       "Sample indexes.")
        .def_readonly(
            "offspring_generation", &ancestry_tracker::generation,
            "Read-only access to current offspring/children generation.")
        .def_readonly("last_gc_time", &ancestry_tracker::last_gc_time,
                      "Last time point where garbage collection happened.")
        .def("update_indexes", &ancestry_tracker::update_indexes)
        .def("prep_for_gc", &ancestry_tracker::prep_for_gc,
             "Call this immediately before you are going to simplify.");

    //Make our C++ function callable from Python.
    //This is NOT part of a user-facing Python API.
    //Rather, we need a wrapper to integrate it with
    //the rest of the fwdpy11 world.
    m.def("evolve_singlepop_regions_track_ancestry",
          &evolve_singlepop_regions_track_ancestry);
    m.def("evolve_singlepop_regions_track_ancestry_async",
          &evolve_singlepop_regions_track_ancestry_async);
    m.def("evolve_singlepop_regions_track_ancestry_python_queue",
          &evolve_singlepop_regions_track_ancestry_python_queue);
}
