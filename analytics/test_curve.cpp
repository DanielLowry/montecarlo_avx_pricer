#include "curve.hpp"
#include <cmath>
#include <iostream>
#include <cassert>
#include <iomanip>

// Simple test framework
#define TEST_CASE(name) std::cout << "\nRunning test case: " << name << std::endl
#define CHECK(condition) do { \
    if (!(condition)) { \
        std::cerr << "Check failed at line " << __LINE__ << ": " << #condition << std::endl; \
        ++g_failed_tests; \
    } else { \
        ++g_passed_tests; \
    } \
} while(0)

#define CHECK_THROWS(expr) do { \
    try { \
        expr; \
        std::cerr << "Expected exception not thrown at line " << __LINE__ << std::endl; \
        ++g_failed_tests; \
    } catch (...) { \
        ++g_passed_tests; \
    } \
} while(0)

#define CHECK_APPROX(a, b, epsilon) do { \
    if (std::abs((a) - (b)) > (epsilon)) { \
        std::cerr << "Approx check failed at line " << __LINE__ << ": " << \
            #a << " (" << std::fixed << std::setprecision(10) << (a) << ") != " << #b << " (" << std::fixed << std::setprecision(10) << (b) << ")" << std::endl; \
        ++g_failed_tests; \
    } else { \
        ++g_passed_tests; \
    } \
} while(0)

static int g_failed_tests = 0;
static int g_passed_tests = 0;

// Helper function to create a simple test curve
std::shared_ptr<discount_curve> create_test_curve(bool flat = true) {
    std::vector<date_t> dates;
    std::vector<double> dfs;
    
    date_t base_date = date_from_int(19723);  // 2024-01-01
    dates.push_back(base_date);  // t0
    dates.push_back(base_date + std::chrono::days(365));   // t1
    dates.push_back(base_date + std::chrono::days(730));   // t2
    dates.push_back(base_date + std::chrono::days(1095));  // t3
    dates.push_back(base_date + std::chrono::days(1460));  // t4

    if (flat) {
        // Flat 5% curve
        dfs.push_back(1.0);
        dfs.push_back(std::exp(-0.05 * 1.0));
        dfs.push_back(std::exp(-0.05 * 2.0));
        dfs.push_back(std::exp(-0.05 * 3.0));
        dfs.push_back(std::exp(-0.05 * 4.0));
    } else {
        // Non-flat curve: e.g., rates of 2%, 3%, 4%, 5% at each node
        dfs.push_back(1.0);  // t0
        dfs.push_back(std::exp(-0.02 * 1.0));  // 2% at t1
        dfs.push_back(std::exp(-0.03 * 2.0));  // 3% at t2
        dfs.push_back(std::exp(-0.04 * 3.0));  // 4% at t3
        dfs.push_back(std::exp(-0.05 * 4.0));  // 5% at t4
    }

    return std::make_shared<discount_curve>(dates, dfs);
}

void test_curve_construction() {
    TEST_CASE("Discount curve construction and validation");
    
    // Valid construction
    try {
        auto curve = create_test_curve();
        ++g_passed_tests;
    } catch (...) {
        std::cerr << "Valid construction failed" << std::endl;
        ++g_failed_tests;
    }

    // Invalid construction - mismatched sizes
    {
        std::vector<date_t> dates{date_from_int(19723)};  // 2024-01-01
        std::vector<double> dfs{1.0, 0.95};
        CHECK_THROWS(discount_curve(dates, dfs));
    }

    // Invalid construction - empty vectors
    {
        std::vector<date_t> dates;
        std::vector<double> dfs;
        CHECK_THROWS(discount_curve(dates, dfs));
    }

    // Invalid construction - unsorted dates
    {
        std::vector<date_t> dates{
            date_from_int(19724),  // 2024-01-02
            date_from_int(19723)   // 2024-01-01
        };
        std::vector<double> dfs{1.0, 0.95};
        CHECK_THROWS(discount_curve(dates, dfs));
    }

    // Invalid construction - first df not 1.0
    {
        std::vector<date_t> dates{
            date_from_int(19723),  // 2024-01-01
            date_from_int(19724)   // 2024-01-02
        };
        std::vector<double> dfs{0.99, 0.95};
        CHECK_THROWS(discount_curve(dates, dfs));
    }
}

void test_discount_factors() {
    TEST_CASE("Discount factor calculations");
    
    auto curve = create_test_curve();
    date_t base_date = curve->valuation_date();

    // Exact node points
    CHECK_APPROX(curve->df(base_date), 1.0, 1e-10);
    CHECK_APPROX(curve->df(base_date + std::chrono::days(365)), 
                 std::exp(-0.05 * 1.0), 1e-10);
    CHECK_APPROX(curve->df(base_date + std::chrono::days(730)), 
                 std::exp(-0.05 * 2.0), 1e-10);

    // Interpolated points
    date_t mid_point = base_date + std::chrono::days(182);
    double expected_df = std::exp(-0.05 * 0.5);
    CHECK_APPROX(curve->df(mid_point), expected_df, 0.001);

    // Before first date
    date_t early_date = base_date - std::chrono::days(10);
    CHECK_APPROX(curve->df(early_date), 1.0, 1e-10);

    // After last date
    date_t late_date = base_date + std::chrono::days(2000);
    CHECK_APPROX(curve->df(late_date), std::exp(-0.05 * 4.0), 1e-10);
}

void test_forward_rates() {
    TEST_CASE("Forward rate calculations");
    
    auto curve = create_test_curve();
    date_t base_date = curve->valuation_date();

    // Forward rate between nodes
    date_t start = base_date;
    date_t end = base_date + std::chrono::days(365);
    double expected_fwd = (1.0 / curve->df(end) - 1.0) / 1.0;  // Calculate based on DF
    CHECK_APPROX(curve->fwd(start, end), expected_fwd, 1e-4);

    // Forward rate for interpolated period (approx 1 year from 182 to 547 days)
    date_t start2 = base_date + std::chrono::days(182);  // 6M
    date_t end2 = base_date + std::chrono::days(547);    // 18M, approx 1 year from start2
    double dcf2 = std::chrono::duration_cast<std::chrono::days>(end2 - start2).count() / 365.0;
    double expected_fwd2 = (curve->df(start2) / curve->df(end2) - 1.0) / dcf2;  // Calculate based on DF
    CHECK_APPROX(curve->fwd(start2, end2), expected_fwd2, 1e-3);

    // Invalid forward rate calculation
    CHECK_THROWS(curve->fwd(base_date, base_date));
}

void test_forward_dfs() {
    TEST_CASE("Forward discount factor calculations");
    
    auto curve = create_test_curve();
    date_t base_date = curve->valuation_date();

    // Forward DF between nodes
    date_t start = base_date;
    date_t end = base_date + std::chrono::days(365);
    double expected_fwd_df = curve->df(end) / curve->df(start);
    CHECK_APPROX(curve->fwd_df(start, end), expected_fwd_df, 1e-10);

    // Forward DF for interpolated period
    date_t start2 = base_date + std::chrono::days(182);
    date_t end2 = base_date + std::chrono::days(547);
    double expected_fwd_df2 = curve->df(end2) / curve->df(start2);
    CHECK_APPROX(curve->fwd_df(start2, end2), expected_fwd_df2, 1e-10);
}

void test_instantaneous_forwards(std::shared_ptr<discount_curve> curve, const std::string& curve_type) {
    TEST_CASE("Instantaneous forward rate calculations for " + curve_type);
    date_t base_date = curve->valuation_date();

    // Instantaneous forward at nodes
    if (curve_type == "flat") {
        CHECK_APPROX(curve->inst_fwd(base_date), 0.05, 1e-3);
        CHECK_APPROX(curve->inst_fwd(base_date + std::chrono::days(365)), 0.05, 1e-3);
    } else {  // Non-flat
        // For non-flat, check against expected rates based on the curve
        CHECK_APPROX(curve->inst_fwd(base_date), 0.02, 1e-3);  // First node's rate
        CHECK_APPROX(curve->inst_fwd(base_date + std::chrono::days(365)), 0.03, 1e-3);  // Second node's rate, approximately
    }

    // Instantaneous forward at interpolated points
    date_t mid_point = base_date + std::chrono::days(182);
    if (curve_type == "flat") {
        CHECK_APPROX(curve->inst_fwd(mid_point), 0.05, 1e-3);
    } else {
        // For non-flat, it should be interpolated, so check it's between expected values
        double inst_rate_mid = curve->inst_fwd(mid_point);
        CHECK(inst_rate_mid > 0.02 && inst_rate_mid < 0.03);  // Between 2% and 3%
    }

    // Consistency with forward rate
    date_t test_date = base_date + std::chrono::days(500);
    double inst_rate = curve->inst_fwd(test_date);
    double fwd_rate = curve->fwd(test_date, test_date + std::chrono::days(1));
    CHECK_APPROX(fwd_rate, inst_rate, 1e-3);  // This should hold for both
}

void test_non_flat_curve() {
    TEST_CASE("Non-flat curve tests");
    auto curve = create_test_curve(false);  // Non-flat curve
    date_t base_date = curve->valuation_date();

    // Test discount factors at nodes
    CHECK_APPROX(curve->df(base_date + std::chrono::days(365)), std::exp(-0.02 * 1.0), 1e-10);
    CHECK_APPROX(curve->df(base_date + std::chrono::days(730)), std::exp(-0.03 * 2.0), 1e-10);
    CHECK_APPROX(curve->df(base_date + std::chrono::days(1095)), std::exp(-0.04 * 3.0), 1e-10);

    // Test interpolated discount factor (e.g., at 6 months)
    date_t mid_point = base_date + std::chrono::days(182);  // Approx halfway to first node
    double interpolated_df = curve->df(mid_point);  // Should be between 1.0 and exp(-0.02*1.0)
    CHECK(interpolated_df > std::exp(-0.02 * 1.0) && interpolated_df < 1.0);

    // Test forward rate between nodes (should reflect increasing rates)
    date_t start = base_date + std::chrono::days(365);  // From t1 to t2
    date_t end = base_date + std::chrono::days(730);
    double expected_fwd = (curve->df(start) / curve->df(end) - 1.0) / 1.0;  // Approx average of 2% and 3%
    CHECK_APPROX(curve->fwd(start, end), expected_fwd, 1e-3);

    // Test forward rate for interpolated period
    date_t start2 = base_date + std::chrono::days(182);
    date_t end2 = base_date + std::chrono::days(547);
    double dcf2 = std::chrono::duration_cast<std::chrono::days>(end2 - start2).count() / 365.0;
    double expected_fwd2 = (curve->df(start2) / curve->df(end2) - 1.0) / dcf2;
    CHECK_APPROX(curve->fwd(start2, end2), expected_fwd2, 1e-3);
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    test_curve_construction();
    test_discount_factors();
    test_forward_rates();
    test_forward_dfs();
    test_instantaneous_forwards(create_test_curve(true), "flat");  // Test flat curve
    test_instantaneous_forwards(create_test_curve(false), "non-flat");  // Test non-flat curve
    test_non_flat_curve();  // Keep this for additional non-flat tests

    std::cout << "\nTest summary:" << std::endl;
    std::cout << "Passed: " << g_passed_tests << std::endl;
    std::cout << "Failed: " << g_failed_tests << std::endl;

    return g_failed_tests > 0 ? 1 : 0;
} 