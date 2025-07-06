#pragma once

#include <vector>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <cmath>


// chrono::sys_days is a new type in C++20 that represents a date
// It is not time-aware, so can only be queried on a date
using date_t = std::chrono::sys_days;

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

    // Helper to get double from date
    static double to_double(const date_t& d);

private:
    std::vector<date_t> node_dates_;
    std::vector<double> discount_factors_;
};
