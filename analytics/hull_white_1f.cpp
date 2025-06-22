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

double hull_white_1f::evolve()