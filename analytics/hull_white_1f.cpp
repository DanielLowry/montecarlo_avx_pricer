
#include "hull_white_1f.hpp"

#include <random>

hull_white_1f::hull_white_1f(const std::shared_ptr<const discount_curve>& curve, double a,
                             double sigma)
    : curve_(curve), a_(a), sigma_(sigma) {}

namespace {  // Quick and dirty standard normal RV generator
// Not vectorised - obvioulsy many improvements can be made here
inline double random_normal() {
    static std::mt19937_64 rng{std::random_device{}()};
    static std::normal_distribution<double> dist{0.0, 1.0};
    return dist(rng);
}

inline double operator-(date_t lhs, date_t rhs) {
    return discount_curve::to_double(lhs) - discount_curve::to_double(rhs);
}

inline date_t operator+(date_t t, double dt) { return t + dt; }

inline date_t operator+(double dt, date_t t) { return t + dt; }

inline date_t operator-(date_t t, double dt) { return t + (-dt); }

}  // namespace

double hull_white_1f::price_cap_monte_carlo(date_t start_date, date_t end_date, double strike,
                                            double notional, int num_paths) const {
    double dcf = (end_date - start_date) / 365.0;  // start date only needed to calculate the DCF
    double df = curve_->df(end_date);

    date_t val_date = curve_->valuation_date();

    double sum_payoffs = 0.0;
    for (int i = 0; i < num_paths; ++i) {
        // For now we just assume daily time steps
        size_t num_steps = end_date - val_date;

        // Run the monte carlo simulation
        double r = curve_->fwd(val_date, val_date + (end_date - start_date));
        for (size_t j = 0; j < num_steps; ++j) r = evolve(val_date + i, r, 1);

        // Caplet payoff
        sum_payoffs += std::max(r - strike, 0.0);
    }

    return df * notional * sum_payoffs / num_paths;
}

double hull_white_1f::price_cap_black(date_t start_date, date_t end_date, double strike,
                                      double notional) const {
    // Simple Black-Scholes formula for an IR caplet
    double fwd = curve_->fwd(start_date, end_date);

    double df = curve_->df(end_date);

    double T = (end_date - start_date) / 365.0;

    if (T <= 0.0 || sigma_ <= 0.0) {
        // No time or no vol, payoff is max(fwd - strike, 0)
        return df * notional * std::max(fwd - strike, 0.0) * T;
    }

    // Black's formula for caplet
    double stddev = sigma_ * std::sqrt(T);
    double d1 = (std::log(fwd / strike) + 0.5 * stddev * stddev) / stddev;
    double d2 = d1 - stddev;

    // Standard normal CDF
    auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };

    double caplet_price = df * notional * T * (fwd * norm_cdf(d1) - strike * norm_cdf(d2));
    return caplet_price;
}

double hull_white_1f::evolve(date_t ti,        // Current time
                             double ri,        // Rate associated with the current time
                             double dt) const  // Time step
{
    // We set theta to be such that the model fits the discount curve
    // It's a reasonably complicated derivation, which leads to the following formula:
    // theta(t) = dfwd(0,t)/dt + a * fwd(0,t) + sigma^2 / (2 * a) * (1 - exp(-2 * a * t))
    // We calculate dfwd(0,t)/dt by finite difference
    constexpr double epsilon = 1e-4;  // for finite difference
    double dfwd_dt =
        (curve_->fwd(ti + epsilon, ti + 2 * epsilon) - curve_->fwd(ti, ti + epsilon)) / epsilon;
    double theta =
        dfwd_dt + a_ * curve_->fwd(ti, ti + epsilon) +
        sigma_ * sigma_ / (2.0 * a_) * (1.0 - std::exp(-2.0 * a_ * discount_curve::to_double(ti)));

    // Standard 1F hull white
    double rn = random_normal();
    return ri + (theta - a_ * ri) * dt + sigma_ * std::sqrt(dt) * rn;
}

std::vector<date_t> double_to_date_vector(const std::vector<double>& double_dates) {
    std::vector<date_t> date_vector;
    date_vector.reserve(double_dates.size());
    
    for (const auto& d : double_dates) {
        // Convert double (days since epoch) to date_t
        date_vector.emplace_back(std::chrono::sys_days{std::chrono::days{static_cast<int>(d)}});
    }
    
    return date_vector;
}


double price_cap_monte_carlo(date_t start_date, date_t end_date, double strike, double notional,
                             int num_paths, double a, double sigma,
                             std::vector<double> curve_node_dates,
                             std::vector<double> curve_node_values) {
    std::shared_ptr<discount_curve> curve =
        std::make_shared<discount_curve>(double_to_date_vector(curve_node_dates), curve_node_values);
    hull_white_1f pricer(curve, a, sigma);
    return pricer.price_cap_monte_carlo(start_date, end_date, strike, notional, num_paths);
}

double price_cap_black(date_t start_date, date_t end_date, double strike, double notional,
                       int num_paths, double a, double sigma, std::vector<double> curve_node_dates,
                       std::vector<double> curve_node_values) {
    std::shared_ptr<discount_curve> curve =
        std::make_shared<discount_curve>(double_to_date_vector(curve_node_dates), curve_node_values);
    hull_white_1f pricer(curve, a, sigma);
    return pricer.price_cap_black(start_date, end_date, strike, notional);
}