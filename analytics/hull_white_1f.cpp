
#include <random>
#include "hull_white_1f.hpp"

hull_white_1f::hull_white_1f(
    const std::shared_ptr<const discount_curve>& curve,
    double a,
    double sigma
)
: curve_(curve)
, a_(a)
, sigma_(sigma)
{}

namespace
{// Quick and dirty standard normal RV generator
    // Not vectorised - obvioulsy many improvements can be made here
    inline double random_normal() {
        static std::mt19937_64 rng{ std::random_device{}() };
        static std::normal_distribution<double> dist{ 0.0, 1.0 };
        return dist(rng);
    }

    using date_t = discount_curve::date_t;
inline date_t operator+(date_t t, double dt) {
    return t + dt;
}

inline date_t operator+(double dt, date_t t) {
    return t + dt;
}

inline date_t operator-(date_t t, double dt) {
    return t + (-dt);
}
}

double hull_white_1f::evolve(
        date_t ti,
		double ri,
		double dt) const
{
    // We set theta to be such that the model fits the discount curve
    // It's a reasonably complicated derivation, which leads to the following formula:
    // theta(t) = dfwd(0,t)/dt + a * fwd(0,t) + sigma^2 / (2 * a) * (1 - exp(-2 * a * t))
    // We calculate dfwd(0,t)/dt by finite difference
    constexpr double epsilon = 1e-4; // for finite difference
    double dfwd_dt = (curve_->fwd(ti + epsilon, ti + 2 * epsilon) - curve_->fwd(ti, ti + epsilon)  ) / epsilon;
    double theta = dfwd_dt + a_ * curve_->fwd(ti, ti + epsilon) + sigma_ * sigma_ / (2.0 * a_) * (1.0 - std::exp(-2.0 * a_ * discount_curve::to_double(ti)));

    // Standard 1F hull white
    return ri + (theta - a_ * ri) * dt + sigma_ * std::sqrt(dt) * random_normal();
}