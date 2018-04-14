#include <future>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "ancestry_tracker.hpp"
#include "evolve_generation.hpp"
#include "evolve_functions.hpp"

namespace py = pybind11;

//Register vectors of nodes, edges, and mutations as "opaque"
PYBIND11_MAKE_OPAQUE(std::vector<node>);
PYBIND11_MAKE_OPAQUE(std::vector<edge>);
PYBIND11_MAKE_OPAQUE(std::vector<mutation>);

PYBIND11_MODULE(wfarg, m)
{
    m.doc() = "Simple example of Wright-Fisher simulation with selection and "
              "ARG tracking";

    //Register nodes and edges as NumPy dtypes:
    PYBIND11_NUMPY_DTYPE(node, id, population, generation);
    PYBIND11_NUMPY_DTYPE(edge, left, right, parent, child);
    PYBIND11_NUMPY_DTYPE(mutation, node_id, pos, mutation_id);

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
        
    py::bind_vector<std::vector<mutation>>(
        m, "MutationArray", "Container of mutations.  This can be cast to a NumPy "
                        "record array without making a copy",
        py::buffer_protocol());

    //Expose the C++ ancestry_tracker to Python.
    //We only expose the stuff that a user really needs
    //to see.
    py::class_<ancestry_tracker>(m, "AncestryTracker")
        .def(py::init<decltype(edge::parent), decltype(edge::parent), decltype(edge::parent)>(),
             py::arg("N"), 
             py::arg("next_index"),
             py::arg("total_generations"))
        .def(py::init<ancestry_tracker&>())
        .def_readwrite("nodes", &ancestry_tracker::nodes,
                       "Data for msprime.NodeTable.")
        .def_readwrite("edges", &ancestry_tracker::edges,
                       "Data for msprime.EdgeTable.")
        .def_readwrite("mutations", &ancestry_tracker::mutations,
                       "Data for msprime.MutationTable and msprime.SiteTable.")
        .def_readonly("node_indexes", &ancestry_tracker::node_indexes,
            "Read-only access to current generation's node start-end indexes.")
        .def("preserve_mutations_sample", &ancestry_tracker::preserve_mutations_sample)
        .def("pre_process_gc", &ancestry_tracker::pre_process_gc)
        .def("post_process_gc", &ancestry_tracker::post_process_gc);
    //    .def("release", [](ancestry_tracker& a) {})
    //    .def("acquire", [](ancestry_tracker& a) {});
    //.def("update_indexes", &ancestry_tracker::update_indexes)
    //.def("prep_for_gc", &ancestry_tracker::prep_for_gc,
    //     "Call this immediately before you are going to simplify.")
    //.def("release_spinlock", &ancestry_tracker::release_spinlock,
    //     "Releases the spin lock.  Used in multi-threaded applications of "
    //     "the msprime machinery.");

//     py::class_<ancestry_data>(
//         m, "_AncestryData",
//         "Used internally for data-swapping in multi-threaded scenarios")
//         .def(py::init<>())
//         .def_readwrite("nodes", &ancestry_data::nodes,
//                        "Data for msprime.NodeTable.")
//         .def_readwrite("edges", &ancestry_data::edges,
//                        "Data for msprime.EdgeTable.")
//         .def_readwrite("mutations", &ancestry_data::mutations,
//                        "Data for msprime.MutationTable and msprime.SiteTable.")
//         .def_readwrite("samples", &ancestry_data::samples, "Sample indexes.")
//         .def("release", [](ancestry_data& a) { a.lock_.attr("release")(); })
//         .def("acquire", [](ancestry_data& a) { a.lock_.attr("acquire")(); });

    //Make our C++ functions callable from Python.
    //This is NOT part of a user-facing Python API.
    //Rather, we need a wrapper to integrate it with
    //the rest of the fwdpy11 world.
    m.def("evolve_singlepop_regions_track_ancestry",
          &evolve_singlepop_regions_track_ancestry);
    /*m.def("evolve_singlepop_regions_track_ancestry_async",
          &evolve_singlepop_regions_track_ancestry_async);
    m.def("evolve_singlepop_regions_track_ancestry_python_queue",
          &evolve_singlepop_regions_track_ancestry_python_queue);*/
}
