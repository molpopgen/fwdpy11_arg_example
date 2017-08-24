#ifndef ANCESTRY_EDGE_HPP__
#define ANCESTRY_EDGE_HPP__

#include <cstdint>
#include <tuple>

struct edge
{
    double left, right;
    std::int32_t parent, child;
};

inline edge
make_edge(double left, double right, std::int32_t parent, std::int32_t child)
{
    edge e;
    e.left = left;
    e.right = right;
    e.parent = parent;
    e.child = child;
    return e;
}

inline auto
get_tied_edge(const edge& e)
    -> decltype(std::tie(e.child, e.parent, e.left, e.right))
{
    return std::tie(e.child, e.parent, e.left, e.right);
}

inline bool
operator<(const edge& lhs, const edge& rhs)
{
    return get_tied_edge(lhs) < get_tied_edge(rhs);
};

inline bool
operator>(const edge& lhs, const edge& rhs)
{
    return get_tied_edge(lhs) > get_tied_edge(rhs);
};
#endif
