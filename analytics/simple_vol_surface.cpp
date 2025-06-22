#include "simple_vol_surface.hpp"

simple_vol_surface::simple_vol_surface(std::vector<date_t> expiries,
                                         std::vector<date_t> tenors,
                                         std::vector<double> strikes,
                                         std::vector<std::vector<std::vector<double>>> vols)
    : expiries_(std::move(expiries)), tenors_(std::move(tenors)), strikes_(std::move(strikes)), vols_(std::move(vols)) {

    if (expiries_.empty() || tenors_.empty() || strikes_.empty() || vols_.empty()) 
        throw std::invalid_argument("Vol surface data cannot be empty");
    
}

double simple_vol_surface::implied_vol(const date_t& expiry,
                                                const date_t& tenor,
                                                double strike) const {

    auto exp_it = std::find(expiries_.begin(), expiries_.end(), expiry);
    if (exp_it == expiries_.end()) 
        throw std::invalid_argument("Expiry not found in vol surface");

    size_t exp_index = std::distance(expiries_.begin(), exp_it);

    auto ten_it = std::find(tenors_.begin(), tenors_.end(), tenor);
    if (ten_it == tenors_.end()) 
        throw std::invalid_argument("Tenor not found in vol surface");

    size_t ten_index = std::distance(tenors_.begin(), ten_it);

    auto strike_it = std::find(strikes_.begin(), strikes_.end(), strike);
    if (strike_it == strikes_.end()) 
        throw std::invalid_argument("Strike not found in vol surface");

    size_t strike_index = std::distance(strikes_.begin(), strike_it);

    return vols_[exp_index][ten_index][strike_index];
}