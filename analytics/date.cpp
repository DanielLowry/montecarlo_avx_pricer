#include "date.hpp"


date_t date_from_int(int dt) {
    return date_t(std::chrono::days(dt));
}

namespace
{
    double sys_date_to_double(const date_t& d)
    {
        return std::chrono::duration_cast<std::chrono::days>(d.time_since_epoch()).count();
    }
}

double date_to_double(const date_t& d) {
    return sys_date_to_double(d);
}

std::string date_to_string(const date_t& d)
{
    return std::to_string(date_to_double(d));
}

std::vector<date_t> double_to_date_vector(const std::vector<double>& double_dates)
{
    std::vector<date_t> date_vector;
    date_vector.reserve(double_dates.size());
    
    for (const auto& d : double_dates) {
        // Convert double (days since epoch) to date_t
        date_vector.emplace_back(std::chrono::sys_days{std::chrono::days{static_cast<int>(d)}});
    }
    
    return date_vector;
}