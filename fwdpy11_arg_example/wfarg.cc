#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <fwdpy11/rules/wf_rules.hpp>
#include <fwdpy11/sim_functions.hpp>
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

void
evolve_singlepop_regions_track_ancestry(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    ancestry_tracker& ancestry, py::function ancestry_processor,
    py::array_t<std::uint32_t> popsizes, const double mu_selected,
    const double recrate, const KTfwd::extensions::discrete_mut_model& mmodel,
    const KTfwd::extensions::discrete_rec_model& rmodel,
    fwdpy11::single_locus_fitness& fitness, const double selfing_rate)
{
    if (pop.generation > 0)
        {
            throw std::runtime_error(
                "this population has already been evolved.");
        }
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    if (recrate < 0.)
        {
            throw std::runtime_error("negative recombination rate: "
                                     + std::to_string(recrate));
        }
    pop.mutations.reserve(
        std::ceil(std::log(2 * pop.N)
                  * (4. * double(pop.N) * (mu_selected)
                     + 0.667 * (4. * double(pop.N) * (mu_selected)))));
    const auto recmap = KTfwd::extensions::bind_drm(
        rmodel, pop.gametes, pop.mutations, rng.get(), recrate);
    const auto mmodels = KTfwd::extensions::bind_dmm(
        mmodel, pop.mutations, pop.mut_lookup, rng.get(), 0.0, mu_selected,
        &pop.generation);
    ++pop.generation;
    auto rules = fwdpy11::wf_rules();

    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    auto wbar = rules.w(pop, fitness_callback);
    //ancestry_tracker ancestry(pop.N);
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            //By hacking the API, we can show that it is
            //passing nodes/edges into Python that is
            //causin a slowdown.  But I *should* be able
            //to do this w/o a copy.  Must investigate!
            py::bool_ processor_rv = ancestry_processor(
                pop.generation, ancestry); //.nodes, ancestry.edges);
            bool gc = processor_rv.cast<bool>();
            if (gc)
                {
                    // ancestry.nodes.clear();
                    // ancestry.edges.clear();
                }
            ancestry.offspring_indexes.clear();
            const auto N_next = popsizes.at(generation);
            evolve_generation(
                rng, pop, N_next, mu_selected, mmodels, recmap,
                std::bind(&fwdpy11::wf_rules::pick1, &rules,
                          std::placeholders::_1, std::placeholders::_2),
                std::bind(&fwdpy11::wf_rules::pick2, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, selfing_rate),
                std::bind(&fwdpy11::wf_rules::update, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, std::placeholders::_4,
                          std::placeholders::_5),
                ancestry, std::true_type());
            pop.N = N_next;
            fwdpy11::update_mutations_wrapper()(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N);
            fitness.update(pop);
            wbar = rules.w(pop, fitness_callback);
        }
    --pop.generation;
    //py::print("Ending simulation with ", ancestry.nodes.size(), "nodes and ",
    //          ancestry.edges.size(), " edges.");
}

//Register vectors of nodes and edges as "opaque"
PYBIND11_MAKE_OPAQUE(std::vector<node>);
PYBIND11_MAKE_OPAQUE(std::vector<edge>);
PYBIND11_MAKE_OPAQUE(std::vector<ancestry_tracker::integer_type>);

PYBIND11_PLUGIN(wfarg)
{
    py::module m("wfarg", "Simple example of Wright-Fisher simulation with "
                          "selection and ARG tracking");

    //Register nodes and edges as NumPy dtypes:
    PYBIND11_NUMPY_DTYPE(node, id, population, generation);
    PYBIND11_NUMPY_DTYPE(edge, left, right, parent, child);

    //py::class_<edge>(m,"Edge");
    //py::class_<node>(m,"Node");

    //Create Python classes of node/edgec containers.
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

    py::class_<ancestry_tracker>(m, "AncestryTracker")
        .def(py::init<KTfwd::uint_t>(), py::arg("N"))
        .def_readwrite("nodes", &ancestry_tracker::nodes,
                       "Data for msprime.NodeTable.")
        .def_readwrite("edges", &ancestry_tracker::edges,
                       "Data for msprime.EdgesetTable.")
        .def_readwrite("samples", &ancestry_tracker::offspring_indexes,
                       "Sample indexes.");
    //Make our C++ function callable from Python.
    //This is NOT part of a user-facing Python API.
    //Rather, we need a wrapper to integrate it with
    //the rest of the fwdpy11 world.
    m.def("evolve_singlepop_regions_track_ancestry",
          &evolve_singlepop_regions_track_ancestry);

    return m.ptr();
}
