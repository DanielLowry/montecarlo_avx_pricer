#pragma once

#include <vector>
#include <chrono>
#include <stdexcept>
#include <algorithm>
#include <cmath>

// A very simple vol surface class.
// Not vectorised
// Industry standard would be to use some sort of SABR model (extended or generalized), but this works for our example
class simple_vol_surface {
public:
    using date_t = std::chrono::sys_days;

    simple_vol_surface(std::vector<date_t> expiries,
                std::vector<date_t> tenors,
                std::vector<double> strikes,
                std::vector<std::vector<std::vector<double>>> vols);

    // fetch the implied vol for a given expiry/tenor/strike
    double implied_vol(const date_t& expiry,
                       const date_t& tenor,
                       double strike) const;

private:
    std::vector<date_t> expiries_, tenors_;
    std::vector<double> strikes_;

    // Stored by expiry/tenor/strike
    std::vector<std::vector<std::vector<double>>> vols_;
};