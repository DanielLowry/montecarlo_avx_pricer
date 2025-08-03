#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hull_white_1f.hpp"

namespace py = pybind11;

double price_cap_monte_carlo(int start_date, int end_date, double strike, double notional,
                             int num_paths, double a, double sigma,
                             std::vector<double> curve_node_dates,
                             std::vector<double> curve_node_values) {
    std::shared_ptr<discount_curve> curve = std::make_shared<discount_curve>(
        double_to_date_vector(curve_node_dates), curve_node_values);
    hull_white_1f pricer(curve, a, sigma);
    return pricer.price_cap_monte_carlo(date_from_int(start_date), date_from_int(end_date), strike, notional, num_paths);
}

double price_cap_black(int start_date, int end_date, double strike, double notional,
                       int num_paths, double a, double sigma, std::vector<double> curve_node_dates,
                       std::vector<double> curve_node_values) {
    std::shared_ptr<discount_curve> curve = std::make_shared<discount_curve>(
        double_to_date_vector(curve_node_dates), curve_node_values);
    hull_white_1f pricer(curve, a, sigma);
    return pricer.price_cap_black(date_from_int(start_date), date_from_int(end_date), strike, notional);
}

PYBIND11_MODULE(mc_pricer_py, m) {
    m.doc() = "Monte Carlo pricer for an interest rate cap";

    m.def("price_cap_mc", &price_cap_monte_carlo, py::arg("start_date"), py::arg("end_date"),
          py::arg("strike"), py::arg("notional"), py::arg("num_paths"), py::arg("a"),
          py::arg("sigma"), py::arg("curve_node_dates"), py::arg("curve_node_values"),
          "Compute the option price via Monte Carlo");

    m.def("price_cap_black", &price_cap_black, py::arg("start_date"), py::arg("end_date"),
          py::arg("strike"), py::arg("notional"), py::arg("num_paths"), py::arg("a"),
          py::arg("sigma"), py::arg("curve_node_dates"), py::arg("curve_node_values"),
          "Compute the option price via Black's formula");
}