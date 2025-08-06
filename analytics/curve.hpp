#pragma once

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "date.hpp"

// Very simple discount curve class - not vectorised for now!
class discount_curve {
public:
    discount_curve(std::vector<date_t> node_dates, std::vector<double> discount_factors);

    date_t valuation_date() const { return node_dates_.front(); }

    // Get discount factor for a given date (linear interpolation)
    double df(const date_t& d) const;

    // Forward discount factor between d1 and d2
    double fwd_df(const date_t& d1, const date_t& d2) const;

    // Forward rate between d1 and d2 (act/365)
    double fwd(const date_t& d1, const date_t& d2) const;

    // Instantaneous forward rate at a given date
    // Approximation using finite difference
    double inst_fwd(const date_t& d) const;

private:
    std::vector<date_t> node_dates_;
    std::vector<double> discount_factors_;

    // Functions taking double instead of date_t arguments.
    // Mostly used for instantaneous forward rate calculations.
    double df(const double d) const;
};
