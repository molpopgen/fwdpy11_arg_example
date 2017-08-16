#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <fwdpy11/types.hpp>

#include "ancestry_edge_sets.hpp"

namespace py = pybind11;

PYBIND11_PLUGIN(wfarg)
{
    py::module m("wfarg", "Simple example of Wright-Fisher simulation with "
                          "selection and ARG tracking");

	return m.ptr();
}
