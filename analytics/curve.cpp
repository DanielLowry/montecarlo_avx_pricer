#include "curve.hpp"
#include "logging.hpp"

discount_curve::discount_curve(std::vector<date_t> node_dates, std::vector<double> discount_factors)
    : node_dates_(std::move(node_dates)), discount_factors_(std::move(discount_factors))
{
    if (node_dates_.size() != discount_factors_.size() || node_dates_.empty()) 
        throw std::invalid_argument("Node dates and discount factors must have the same non-zero size.");
    
    if (!std::is_sorted(node_dates_.begin(), node_dates_.end())) 
        throw std::invalid_argument("Node dates must be sorted in ascending order.");    

    if(discount_factors_.front() != 1.0)
        throw std::invalid_argument("First discount factor must be 1.0.");
}



double discount_curve::df(const date_t& d) const
{
    if (d <= node_dates_.front()) 
        return discount_factors_.front();    

    if (d >= node_dates_.back()) 
        return discount_factors_.back();    

    const auto it = std::upper_bound(node_dates_.begin(), node_dates_.end(), d);
    size_t idx = std::distance(node_dates_.begin(), it) - 1;

    const auto& d0 = node_dates_[idx];
    const auto& d1 = node_dates_[idx + 1];

    double t0 = date_to_double(d0);
    double t1 = date_to_double(d1);
    double t = date_to_double(d);

    double w = (t - t0) / (t1 - t0);
    double v0 = discount_factors_[idx];
    double v1 = discount_factors_[idx + 1];

    // Geometric (implies constant continuous rate between nodes)
    return std::exp(std::log(v0) + w * (std::log(v1) - std::log(v0)));
}

double discount_curve::fwd_df(const date_t& d1, const date_t& d2) const
{
    return df(d2) / df(d1);
}

double discount_curve::fwd(const date_t& d1, const date_t& d2) const
{
    // LOG("Calculating forward rate for " << date_to_string(d1) << " to " 
    //     << date_to_string(d2));

    double df1 = df(d1);
    double df2 = df(d2);

    if (df1 <= 0.0 || df2 <= 0.0) 
        throw std::domain_error("Discount factors must be positive.");
    
    double days = std::chrono::duration_cast<std::chrono::days>(d2 - d1).count();
    if (days <= 0) 
        throw std::invalid_argument("d2 must be after d1.");
    
    // Assume ACT/365 convention
    double dcf = days / 365.0;
    return (df1 / df2 - 1.0) / dcf;
}

double discount_curve::df(const double d) const
{
    // Interpolate discount factor for a time expressed as double (days since epoch)
    double t = d;

    // Get node times as doubles
    double t_first = date_to_double(node_dates_.front());
    double t_last = date_to_double(node_dates_.back());

    if (t <= t_first)
        return discount_factors_.front();

    if (t >= t_last)
        return discount_factors_.back();

    // Find the interval [t0, t1] such that t0 <= t < t1
    // Linear search is fine for small number of nodes
    size_t idx = 0;
    for (size_t i = 0; i < node_dates_.size() - 1; ++i) {
        double t0 = date_to_double(node_dates_[i]);
        double t1 = date_to_double(node_dates_[i + 1]);
        if (t >= t0 && t < t1) {
            idx = i;
            break;
        }
    }

    double t0 = date_to_double(node_dates_[idx]);
    double t1 = date_to_double(node_dates_[idx + 1]);
    double w = (t - t0) / (t1 - t0);

    double v0 = discount_factors_[idx];
    double v1 = discount_factors_[idx + 1];

    // Geometric interpolation (constant continuous rate between nodes)
    return std::exp(std::log(v0) + w * (std::log(v1) - std::log(v0)));
}

double discount_curve::inst_fwd(const date_t& d) const
{
    constexpr double epsilon = 1e-4;
    double t = date_to_double(d);
    double t = date_to_double(d);

    // If at or before the first node, use forward difference
    if (d == node_dates_.front()) {
        double df1 = df(t);
        double df2 = df(t + epsilon);
        LOG("df1=" << df1 << ", df2=" << df2);
        if (df1 <= 0.0 || df2 <= 0.0)
            throw std::domain_error("Discount factors must be positive.");
        constexpr double dcf = epsilon / 365.0;
        return (df1 / df2 - 1.0) / dcf;
    } else {
        // Central difference
        double df_minus = df(t - epsilon);
        double df_plus = df(t + epsilon);
        if (df_minus <= 0.0 || df_plus <= 0.0)
            throw std::domain_error("Discount factors must be positive.");
        constexpr double dcf = 2 * epsilon / 365.0;
        return (df_minus / df_plus - 1.0) / dcf;
    }
}