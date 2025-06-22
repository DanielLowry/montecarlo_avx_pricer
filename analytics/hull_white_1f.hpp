
#include <memory>

#include "curve.hpp"
#include "simple_vol_surface.hpp"

/// Simple 1 factor hull white model

// Hull-White 1-factor model: dr(t) = [theta(t) - a*r(t)]dt + sigma*dW(t)
// a: mean reversion, sigma: volatility, theta(t): time-dependent drift, used to ensure the model fits the discount curve


class hull_white_1f {
public:
	using date_t = discount_curve::date_t;

	// If you already know the vol, you can construct the model directly
	// Otherwise, you can calibrate the model from a discount curve and vol surface
	hull_white_1f(
		const std::shared_ptr<const discount_curve>& curve,
		double a,
		double sigma
	);

	// We want to be able to calibrate the hull white 1F model from a discount curve and vol surface
	// No implementation yet
	static std::shared_ptr<hull_white_1f> calibrate(
		const std::shared_ptr<const discount_curve>& curve,
		const std::shared_ptr<simple_vol_surface>& vol_surf,
		double a = 0.01,
		double tol = 1e-8
	);

	double a() const { return a_; }
	double sigma() const { return sigma_; }



private:

	std::shared_ptr<const discount_curve> curve_;
	double a_, sigma_;

	// evolve method to help simulate monte carlo paths
	double evolve(
        date_t ti,
		double ri,
		double dt
	) const;

};