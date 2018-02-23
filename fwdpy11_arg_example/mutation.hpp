#ifndef ANCESTRY_MUT_HPP__
#define ANCESTRY_MUT_HPP__

#include <cstdint>

struct mutation
{
    std::int32_t node_id;
    std::uint32_t mutation_id;
};

inline mutation
make_mutation(std::int32_t node_id, std::uint32_t mutation_id)
{
    mutation m;
    m.node_id = node_id;
    m.mutation_id = mutation_id;
    return m;
}

#endif
