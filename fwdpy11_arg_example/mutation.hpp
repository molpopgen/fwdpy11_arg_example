#ifndef ANCESTRY_MUT_HPP__
#define ANCESTRY_MUT_HPP__

#include <cstdint>

struct mutation
{
    double position, origin_generation;
    std::int32_t node_id;
};

inline mutation
make_mutation(double position, double origin_generation, std::int32_t node_id)
{
    mutation m;
    m.position = position;
    m.origin_generation = origin_generation;
    m.node_id = node_id;
    return m;
}

#endif
