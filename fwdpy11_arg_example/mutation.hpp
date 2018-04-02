#ifndef ANCESTRY_MUT_HPP__
#define ANCESTRY_MUT_HPP__

#include <cstdint>

struct mutation
{
    std::int32_t node_id;
    double pos;
    std::uint32_t mutation_id;
};

#endif
