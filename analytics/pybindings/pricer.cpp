#include <pybind11/pybind11.h>
#include "hull_white_1f.hpp" 

namespace py = pybind11;

PYBIND11_MODULE(mc_pricer_py, m) {
    m.doc() = "Monte Carlo pricer for an interest rate cap";

    m.def("price", &price_cap_monte_carlo,
        py::arg("start_date"),
        py::arg("end_date"),
        py::arg("strike"),
         py::arg("notional"),
          py::arg("num_paths"),
           py::arg("a"),
           py::arg("sigma"),
           py::arg("curve_node_dates"),
           py::arg("curve_node_values"),
        "Compute the option price via Monte Carlo");

    m.def("price_black", &price_cap_black,
        py::arg("start_date"),
        py::arg("end_date"),
        py::arg("strike"),
        py::arg("notional"),
        py::arg("num_paths"),
        py::arg("a"),
        py::arg("sigma"),
        py::arg("curve_node_dates"),
        py::arg("curve_node_values"),
        "Compute the option price via Black's formula");
}