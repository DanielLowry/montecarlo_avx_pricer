#include <chrono>

// chrono::sys_days is a new type in C++20 that represents a date
// It is not time-aware, so can only be queried on a date
using date_t = std::chrono::sys_days;

date_t date_from_int(int dt);

// Helper to get double from date
double date_to_double(const date_t& d);

std::string date_to_string(const date_t& d);

std::vector<date_t> double_to_date_vector(const std::vector<double>& double_dates);