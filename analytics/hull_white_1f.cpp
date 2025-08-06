
#include "hull_white_1f.hpp"

#include <random>

#include "logging.hpp"

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
    return date_to_double(lhs) - date_to_double(rhs);
}

inline date_t operator+(date_t t, double dt) {
    // Add dt days to t
    return t + std::chrono::days(static_cast<int>(dt));
}

inline date_t operator+(double dt, date_t t) { return t + dt; }

inline date_t operator-(date_t t, double dt) { return t + (-dt); }

}  // namespace

double hull_white_1f::price_cap_black(date_t start_date, date_t end_date, double strike,
                                      double notional) const {
    // Simple Black-Scholes formula for an IR caplet
    double fwd = curve_->fwd(start_date, end_date);
    LOG("Black pricing: fwd=" << fwd);

    // We don't yet handle negative rates
    if (fwd <= 0.0) {
        throw std::runtime_error("Negative forward rates are not yet supported");
    }

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
    LOG("Black pricing: d1=" << d1 << ", d2=" << d2);

    // Standard normal CDF
    auto norm_cdf = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
    double n_d1 = norm_cdf(d1);
    double n_d2 = norm_cdf(d2);
    LOG("Black pricing: N(d1)=" << n_d1 << ", N(d2)=" << n_d2);

    double caplet_price = df * notional * T * (fwd * n_d1 - strike * n_d2);
    return caplet_price;
}

double hull_white_1f::price_cap_monte_carlo(date_t start_date, date_t end_date,
                                            double strike, double notional,
                                            int num_paths) const {
    double dcf =
        (end_date - start_date) / 365.0;  // start date only needed to calculate the DCF
    double df = curve_->df(end_date);

    date_t val_date = curve_->valuation_date();
    if (end_date <= val_date)
        throw std::runtime_error("end_date must be after valuation_date");

    size_t num_steps = start_date - val_date;
    double initial_rate = curve_->inst_fwd(val_date);
    //double initial_rate = curve_->fwd(val_date, val_date + (end_date - start_date));
    LOG("MC pricing: initial_rate=" << initial_rate);

    double sum_payoffs = 0.0;
    // Log just the first path to see rate evolution
    double r = initial_rate;
    for (size_t j = 0; j < num_steps; ++j) {
        r = evolve(val_date + j, r, 1.0 / 365.0);
    }
    LOG("MC pricing: first path final_rate=" << r
                                             << ", payoff=" << std::max(r - strike, 0.0));

    // Run remaining paths
    for (int i = 1; i < num_paths; ++i) {
        r = initial_rate;
        for (size_t j = 0; j < num_steps; ++j)
            r = evolve(val_date + j, r, 1.0 / 365.0 );
        sum_payoffs += std::max(r - strike, 0.0);
    }

    double avg_payoff = sum_payoffs / num_paths;
    LOG("MC pricing: avg_payoff=" << avg_payoff);
    return df * dcf * notional * avg_payoff;
}

double hull_white_1f::evolve(date_t ti,        // Current time
                             double ri,        // Rate associated with the current time
                             double dt) const  // Time step
{
    // We set theta to be such that the model fits the discount curve
    // It's a reasonably complicated derivation, which leads to the following formula:
    // theta(t) = dfwd(0,t)/dt + a * fwd(0,t) + sigma^2 / (2 * a) * (1 - exp(-2 * a * t))

    constexpr int epsilon_days = 1;
    constexpr double epsilon = epsilon_days / 365.0;
    double f0 = curve_->inst_fwd(ti); 
    double f1 = curve_->inst_fwd(ti+epsilon_days);
    double dfwd_dt = (f1 - f0) / epsilon;

    double theta = dfwd_dt + a_ * curve_->inst_fwd(ti) +
                   sigma_ * sigma_ / (2.0 * a_) *
                       (1.0 - std::exp(-2.0 * a_ * date_to_double(ti) / 365.0));

    // Standard 1F hull white
    double rn = random_normal();
    return ri + (theta - a_ * ri) * dt + sigma_ * std::sqrt(dt) * rn;
}

double hull_white_1f::black_vol_to_hw_vol(double black_vol, double a, double dcf,
                                          double T) {
    return black_vol * std::sqrt((2 * a * T) / (1 - std::exp(-2 * a * T))) * (a * dcf) /
           (1 - std::exp(-a * dcf));
}
