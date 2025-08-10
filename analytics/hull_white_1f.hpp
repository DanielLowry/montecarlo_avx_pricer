
#include <memory>

#include "curve.hpp"
#include "simple_vol_surface.hpp"

/// Simple 1 factor hull white model

// Hull-White 1-factor model: dr(t) = [theta(t) - a*r(t)]dt + sigma*dW(t)
// a: mean reversion, sigma: volatility, theta(t): time-dependent drift,
// used to ensure the model fits the discount curve

class hull_white_1f {
   public:
    // Remove using date_t = discount_curve::date_t; and use global date_t

    // If you already know the vol, you can construct the model directly
    // Otherwise, you can calibrate the model from a discount curve and vol
    // surface
    hull_white_1f(const std::shared_ptr<const discount_curve>& curve, double a, double sigma);

    // We want to be able to calibrate the hull white 1F model from a
    // discount curve and vol surface No implementation yet
    static std::shared_ptr<hull_white_1f> calibrate(
        const std::shared_ptr<const discount_curve>& curve,
        const std::shared_ptr<simple_vol_surface>& vol_surf, double a = 0.01, double tol = 1e-8);

    double a() const { return a_; }
    double sigma() const { return sigma_; }

    // Monte Carlo caplet pricing under 1F Hull-White
    // - start_date: accrual period start (reset) date S
    // - end_date: accrual period end (payment) date E
    // - strike: caplet strike K (simple-compounded)
    // - notional: notional N
    // - num_paths: number of Monte Carlo paths to simulate
    double price_cap_monte_carlo(date_t start_date, date_t end_date, double strike, double notional,
                                 int num_paths) const;

    double price_cap_black(date_t start_date, date_t end_date, double strike,
                           double notional) const;

    // Convert Black caplet vol to equivalent Hull-White vol
    // - black_vol: Black vol for tenor dcf with maturity T
    // - a: mean-reversion
    // - dcf: accrual year fraction of the caplet period
    // - T: time to option maturity in years
    static double black_vol_to_hw_vol(double black_vol, double a, double dcf, double T);

   private:
    std::shared_ptr<const discount_curve> curve_;
    double a_, sigma_;

    // Evolve short rate one Euler step
    // - ti: current date
    // - ri: short rate r(ti)
    // - dt: time step in years
    double evolve(date_t ti, double ri, double dt) const;

    // Hull-White analytic helpers for bond price P(t,T) = A(t,T) * exp(-B(t,T) * r_t)
    // - t: time in years from valuation date to the evaluation time
    // - T: time in years from valuation date to maturity
    double hw_B(double t, double T) const;
    
    // - t_date: date corresponding to time t
    // - T_date: maturity date
    double hw_a(double t, date_t t_date, date_t maturity_date) const;

    // Struct holding pathwise info up to a given target date
    struct short_rate_path_info {
        date_t valuation_date;      // t=0
        date_t start_date;          // target date S
        double r_at_start;          // r(S) along the path
        double integral_r_dt;       // \int_0^S r(s) ds along the path
        const hull_white_1f& model; // non-owning reference to model for A/B/curve

        // Simple-compounded forward LIBOR L(S,E) implied by the path at S
        // - end_date: period end date E
        double libor(date_t end_date) const;

        // Pathwise discount factor from 0 to S: exp(-\int_0^S r(s) ds)
        double discount_0_to_start() const { return std::exp(-integral_r_dt); }
    };

    // Simulate short rate up to target_date, returning path info
    // - target_date: date S to reach
    // - r0: initial short rate at valuation date
    // - dt: time step in years
    short_rate_path_info simulate_short_rate_to(date_t target_date, double r0, double dt) const;
};


                             