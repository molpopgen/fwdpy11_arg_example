#ifndef WF_RULES_ASYNC_HPP__
#define WF_RULES_ASYNC_HPP__

#include <vector>
#include <future>
#include <cstddef>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpy11/rules/wf_rules.hpp>


struct wf_rules_async_fitness : public fwdpy11::wf_rules
{
    int nthreads;
    std::vector<std::future<double>> fitness_futures;
    struct accumulate_fitnesses
    // Based on ideas from https://github.com/molpopgen/fwdpp_ideas/blob/master/parallel.md
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
                        first_diploid = last_diploid_in_range;
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

#endif
