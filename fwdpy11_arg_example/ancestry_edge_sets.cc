#include <vector>
#include <utility>
#include <tuple>

inline bool
overlap_detail(const std::tuple<double, double>& lhs,
               const std::tuple<double, double>& rhs)
{
    // return not (x[1] < y[0] or y[1] < x[0])
    return !((std::get<1>(lhs) <= std::get<0>(rhs))
             || (std::get<1>(rhs) <= std::get<0>(lhs)));
}

bool
overlap(const std::tuple<double, double>& lhs,
        const std::tuple<double, double>& rhs)
{
    if (std::get<0>(lhs) < std::get<0>(rhs))
        {
            return overlap_detail(lhs, rhs);
        }
    return overlap_detail(rhs, lhs);
}

