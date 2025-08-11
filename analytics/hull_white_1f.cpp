
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
    double fwd_rate = curve_->fwd(start_date, end_date);
    LOG("Black pricing: fwd_rate=" << fwd_rate);

    // We don't yet handle negative rates
    if (fwd_rate <= 0.0) {
        throw std::runtime_error("Negative forward rates are not yet supported");
    }

    // Discount to payment date and separate accrual vs. option maturity
    double discount_factor = curve_->df(end_date);
    const date_t val_date = curve_->valuation_date();
    const double tau = (end_date - start_date) / 365.0;

    if (start_date <= val_date) {
        // No time to option maturity or no vol: intrinsic value at reset discounted to
        // payment
        LOG("Black pricing: no time to option maturity or no vol: intrinsic value at "
            "reset discounted to payment");
        return discount_factor * notional * tau * std::max(fwd_rate - strike, 0.0);
    }

    // Black's formula for caplet with time to reset
    double t = (start_date - val_date) / 365.0;
    double volatility_stddev = sigma_ * std::sqrt(t);
    double d1 =
        (std::log(fwd_rate / strike) + 0.5 * volatility_stddev * volatility_stddev) /
        volatility_stddev;
    double d2 = d1 - volatility_stddev;
    LOG("Black pricing: d1=" << d1 << ", d2=" << d2);

    // Standard normal CDF
    auto normal_cdf = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
    double n_d1 = normal_cdf(d1);
    double n_d2 = normal_cdf(d2);
    LOG("Black pricing: N(d1)=" << n_d1 << ", N(d2)=" << n_d2);

    double caplet_price =
        discount_factor * notional * tau * (fwd_rate * n_d1 - strike * n_d2);
    return caplet_price;
}

double hull_white_1f::price_cap_monte_carlo(date_t start_date, date_t end_date,
                                            double strike, double notional,
                                            int num_paths) const {
    // Accrual factor for the caplet period [start_date, end_date]
    const double dcf = (end_date - start_date) / 365.0;

    // Valuation date (t=0)
    const date_t val_date = curve_->valuation_date();
    if (end_date <= val_date) {
        throw std::runtime_error("end_date must be after valuation_date");
    }

    // Number of Euler steps (1 day granularity)
    const size_t steps_to_start = start_date - val_date;
    const size_t steps_to_end = end_date - val_date;
    const double time_step = 1.0 / 365.0;

    // Initial short rate approximation from the curve
    const double initial_rate = curve_->inst_fwd(val_date);
    LOG("MC pricing: initial_rate=" << initial_rate);

    double total_discounted_payoffs = 0.0;

    for (int i = 0; i < num_paths; ++i) {
        short_rate_path_info path =
            simulate_short_rate_to(start_date, end_date, initial_rate, time_step);

        const double libor_rate = path.libor(end_date);
        const double payoff_at_start =
            std::max(libor_rate - strike, 0.0) * dcf * notional;
        const double discount_to_end = path.discount_0_to_end();
        total_discounted_payoffs += discount_to_end * payoff_at_start;

        if (i == 0) {
            LOG("MC pricing: first path r_start=" << path.r_at_start
                                                  << ", L=" << libor_rate
                                                  << ", payoff_start=" << payoff_at_start
                                                  << ", disc(0,end)=" << discount_to_end);
        }
    }

    const double price = total_discounted_payoffs / static_cast<double>(num_paths);
    LOG("MC pricing: avg discounted payoff=" << price);
    return price;
}

double hull_white_1f::hw_B(double t, double T) const {
    const double x = a_ * (T - t);
    return (1.0 - std::exp(-x)) / a_;
}

double hull_white_1f::hw_a(double t, date_t t_date, date_t maturity_date) const {
    const date_t val_date = curve_->valuation_date();
    const double p0_t = curve_->df(t_date);
    const double p0_maturity = curve_->df(maturity_date);
    const double f0t = curve_->inst_fwd(t_date);
    const double b = hw_B(t, (maturity_date - val_date) / 365.0);
    // Var-term in A(t,T) contains an integral over B(s,T)^2 from 0..t
    // For HW with constant a and sigma: \int_0^t B(s,T)^2 ds = (b^2) * (1 - e^{-2 a t}) /
    // (2 a) leading to (sigma^2 / 2) * \int_0^t B(s,T)^2 ds = (sigma^2 / (4 a)) * (1 -
    // e^{-2 a t}) * b^2
    const double variance_term =
        (sigma_ * sigma_) / (4.0 * a_) * (1.0 - std::exp(-2.0 * a_ * t)) * (b * b);
    return (p0_maturity / p0_t) * std::exp(b * f0t - variance_term);
}

auto hull_white_1f::simulate_short_rate_to(date_t start_date, date_t end_date, double r0,
                                           double dt) const -> short_rate_path_info {
    const date_t val_date = curve_->valuation_date();
    const size_t steps_to_start = start_date - val_date;
    const size_t steps_to_end = end_date - val_date;

    double current_rate = r0;
    double integral_to_start = 0.0;
    double integral_to_end = 0.0;
    double rate_at_start = r0;

    for (size_t j = 0; j < steps_to_end; ++j) {
        if (j < steps_to_start)
            integral_to_end += current_rate * dt;
        current_rate = evolve(val_date + j, current_rate, dt);
        if (j + 1 == steps_to_start)
            rate_at_start = current_rate;
    }

    return short_rate_path_info{val_date,      start_date,      end_date,
                                rate_at_start, integral_to_end, *this};
}

// Convert the pathwise short rate information at start_date (S) into the
// simple-compounded forward LIBOR L(S,E) for the period [S, E].
// Rationale:
//   Under HW1F, the zero-coupon bond price at S for maturity E is
//     P(S,E) = A(S,E) * exp(-B(S,E) * r(S)).
//   Given the path's short rate r(S), we compute P(S,E), then convert the
//   discount factor to the simple forward via L(S,E) = (1/P(S,E) - 1) / dcf.
double hull_white_1f::short_rate_path_info::libor(date_t end_date) const {
    // Times in years from valuation date (t=0)
    const double t_start = (start_date - valuation_date) / 365.0;
    const double t_end = (end_date - valuation_date) / 365.0;
    const double dcf = (end_date - start_date) / 365.0;  // accrual year fraction [S,E]

    // HW1F bond pricing coefficients evaluated at S for maturity E
    const double b = model.hw_B(t_start, t_end);                      // B(S,E)
    const double a_term = model.hw_a(t_start, start_date, end_date);  // A(S,E)

    // Pathwise discount factor P(S,E)
    const double p_start_end = a_term * std::exp(-b * r_at_start);

    // Simple-compounded forward from P(S,E)
    return (1.0 / p_start_end - 1.0) / dcf;
}

double hull_white_1f::evolve(date_t ti,        // Current time
                             double ri,        // Rate associated with the current time
                             double dt) const  // Time step
{
    // We set theta to be such that the model fits the discount curve
    // It's a reasonably complicated derivation, which leads to the following formula:
    // theta(t) = dfwd(0,t)/dt + a * fwd(0,t) + sigma^2 / (2 * a) * (1 - exp(-2 * a * t))

    constexpr int epsilon_days_local = 1;
    constexpr double epsilon_time = epsilon_days_local / 365.0;
    double inst_fwd_rate = curve_->inst_fwd(ti);
    double inst_fwd_rate_next = curve_->inst_fwd(ti + epsilon_days_local);
    double dfwd_dt = (inst_fwd_rate_next - inst_fwd_rate) / epsilon_time;

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
