#include <future>
#include <chrono>
#include <pybind11/chrono.h>
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
#include "ancestry_tracker.hpp"
#include "evolve_generation.hpp"

namespace py = pybind11;

// This function runs the simulation itself.
// The details of a generation are in the file
// evolve_generation.hpp. This function gets
// exposed to Python, but it is not called
// directly by a user. The function "evolve_track"
// in evolve_arg.py takes the fwdpy11 objects and
// passes them along to this function. The main purpose
// of this function is to check parameters and send
// data to msprime when needed.
// The argument ancestry_processor is a Python callable.
// It is to handle the GC/simplification step via msprime.
// It should be an instance of ARGsimplifier.
// The return value is the time spent simulating.
double
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

    double time_simulating = 0.0;
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            //Ask if we need to garbage collect:
            py::tuple processor_rv
                = ancestry_processor(pop.generation, ancestry);
            //If we did GC, then the ancestry_tracker has
            //some cleaning upto do:
            ancestry.post_process_gc(processor_rv);

            //This is not great API design, but
            //we need to clear the offspring indexes here:
            ancestry.offspring_indexes.clear();
            const auto N_next = popsizes.at(generation);
            auto start = std::clock();
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
            fwdpy11::update_mutations(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N, true);
            fitness.update(pop);
            wbar = rules.w(pop, fitness_callback);
            auto stop = std::clock();
            auto dur = (stop - start) / (double)CLOCKS_PER_SEC;
            time_simulating += dur;
        }
    --pop.generation;
    return time_simulating;
}

struct wf_rules_async_fitness : public fwdpy11::wf_rules
{
    int nthreads;
    std::vector<std::future<double>> fitness_futures;
    struct accumulate_fitnesses
    // This version directly fills in the fitness vector.
    // The potential issue with this approach is contention
    // due to false sharing at the boundaries of ranges
    // updated by different threads.
    {
        template <typename itr, typename gcont, typename mcont, typename wfunc>
        inline double
        operator()(const gcont& gametes, const mcont& mutations,
                   std::vector<double>& w, std::size_t i, itr beg, itr end,
                   const wfunc& ff)
        {
            double sum = 0.0;
            for (; beg < end; ++beg, ++i)
                {
                    double result = ff(*beg, gametes, mutations);
                    beg->w = result;
                    beg->g = result;
                    w[i] = result;
                    sum += result;
                }
            return sum;
        }
    };

    wf_rules_async_fitness(int nthreads_ = 1)
        : wf_rules(), nthreads(nthreads_), fitness_futures{}
    {
        if (nthreads > 1)
            {
                fitness_futures.resize(nthreads - 1);
            }
    }

    double
    w(fwdpy11::singlepop_t& pop, const fwdpy11::single_locus_fitness_fxn& ff)
    {
        double wbar = 0.0;
        if (nthreads > 1)
            {
                auto N_curr = pop.diploids.size();
				fitnesses.resize(N_curr);
                std::size_t increment = N_curr / nthreads + 1;
                std::size_t offset = increment;
                auto first_diploid = pop.diploids.begin() + offset;
                for (std::size_t i = 0;
                     i < static_cast<std::size_t>(nthreads) - 1; ++i)
                    {
                        auto last_diploid_in_range = std::min(
                            pop.diploids.end(), first_diploid + increment);
                        // operator()(const gcont& gametes, const mcont& mutations,
                        //            std::vector<double>& w, std::size_t i, itr beg, itr end,
                        //            const wfunc& ff)
                        fitness_futures[i] = std::async(
                            std::launch::async,
                            std::bind(
                                accumulate_fitnesses(), std::cref(pop.gametes),
                                std::cref(pop.mutations), std::ref(fitnesses),
                                static_cast<std::size_t>(std::distance(
                                    pop.diploids.begin(), first_diploid)),
                                first_diploid, last_diploid_in_range, ff));
						first_diploid=last_diploid_in_range;
                    }
                wbar = accumulate_fitnesses()(
                    pop.gametes, pop.mutations, fitnesses, 0,
                    pop.diploids.begin(), pop.diploids.begin() + offset, ff);
                for (auto& f : fitness_futures)
                    {
                        wbar += f.get();
                    }
                wbar /= static_cast<double>(N_curr);
                lookup = KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr(
                    gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
            }
        else
            {
                fwdpy11::wf_rules::w(pop, ff);
            }
        return wbar;
    }
};

double
evolve_singlepop_regions_track_ancestry_async(
    const fwdpy11::GSLrng_t& rng, fwdpy11::singlepop_t& pop,
    ancestry_tracker& ancestry, py::function ancestry_processor,
    const KTfwd::uint_t gc_interval, py::array_t<std::uint32_t> popsizes,
    const double mu_selected, const double recrate,
    const KTfwd::extensions::discrete_mut_model& mmodel,
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
    auto rules = wf_rules_async_fitness(2);

    auto fitness_callback = fitness.callback();
    fitness.update(pop);
    auto wbar = rules.w(pop, fitness_callback);

    std::future<py::object> msprime_future;
    ancestry_tracker local_ancestry_tracker(ancestry);
    double time_simulating = 0.0;
    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            if (pop.generation > 0 && pop.generation % gc_interval == 0.)
                {
                    if (msprime_future.valid())
                        {
                            msprime_future.wait();
                            auto result = msprime_future.get();
                            auto result_tuple = result.cast<py::tuple>();
                            ancestry.post_process_gc(result_tuple, false);
                        }
                    {
                        py::gil_scoped_acquire acquire;
                        py::print("We got our future:", pop.generation,
                                  ancestry.nodes[0].id,
                                  ancestry.nodes.back().id,
                                  ancestry.first_parental_index,
                                  ancestry.next_index);
                    }
                    ancestry.exchange_for_async(local_ancestry_tracker);
                    //auto async_data = ancestry.prep_for_async();
                    py::print(local_ancestry_tracker.nodes.size(),
                              local_ancestry_tracker.edges.size(),
                              local_ancestry_tracker.offspring_indexes.size(),
                              ancestry.nodes.size(), ancestry.edges.size(),
                              ancestry.nodes.capacity());
                    msprime_future = std::async(
                        std::launch::async, ancestry_processor, pop.generation,
                        std::ref(local_ancestry_tracker));
                    //msprime_future.wait();
                    //auto result = msprime_future.get();
                    //auto result_tuple = result.cast<py::tuple>();
                    //ancestry.post_process_gc(result_tuple);

                    //If we did GC, then the ancestry_tracker has
                    //some cleaning upto do:
                    //ancestry.post_process_gc(
                    //    processor_rv_future.get().cast<py::tuple>());
                }
            //This is not great API design, but
            //we need to clear the offspring indexes here:
            ancestry.offspring_indexes.clear();
            const auto N_next = popsizes.at(generation);
            auto start = std::clock();
            evolve_generation(
                rng, pop, N_next, mu_selected, mmodels, recmap,
                std::bind(&decltype(rules)::pick1, &rules,
                          std::placeholders::_1, std::placeholders::_2),
                std::bind(&decltype(rules)::wf_rules::pick2, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, selfing_rate),
                std::bind(&decltype(rules)::update, &rules,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, std::placeholders::_4,
                          std::placeholders::_5),
                ancestry, std::true_type());
            pop.N = N_next;
            fwdpy11::update_mutations(
                pop.mutations, pop.fixations, pop.fixation_times,
                pop.mut_lookup, pop.mcounts, pop.generation, 2 * pop.N, true);
            fitness.update(pop);
            wbar = rules.w(pop, fitness_callback);
            auto stop = std::clock();
            auto dur = (stop - start) / (double)CLOCKS_PER_SEC;
            time_simulating += dur;
        }
    py::print("leaving sim. future state: ", msprime_future.valid());
    --pop.generation;
    if (msprime_future.valid())
        {
            msprime_future.wait();
            auto result = msprime_future.get();
            auto result_tuple = result.cast<py::tuple>();
            ancestry.post_process_gc(result_tuple, false);
        }

    return time_simulating;
}
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
}
