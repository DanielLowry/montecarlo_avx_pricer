import unittest
from mc_pricer_py import price_cap_mc, price_cap_black
import numpy as np

class TestCapPricing(unittest.TestCase):
    def setUp(self):
        # Base case parameters (kept for backward compatibility)
        self.start_date = 1  # in days
        self.end_date = 500  # in days
        self.strike = 0.02
        self.notional = 1000000
        self.num_paths = 5000  # Reduced from 10000
        self.a = 0.008
        self.sigma = 0.01  # black vol
        self.curve_node_dates = [1, 500]
        self.curve_node_values = [1.0, 0.98]

    def _run_mc_vs_black_comparison(self, start_date, end_date, strike, notional,
                                  num_paths, a, sigma, curve_node_dates, curve_node_values,
                                  rel_tol=1e-2):
        """Helper method to compare Monte Carlo vs Black prices"""
        result_black = price_cap_black(
            start_date, end_date, strike, notional,
            num_paths, a, sigma,
            curve_node_dates, curve_node_values
        )

        result_mc = price_cap_mc(
            start_date, end_date, strike, notional,
            num_paths, a, sigma,
            curve_node_dates, curve_node_values
        )

        # Handle zero and near-zero price cases
        if abs(result_black) < 1e-10 and abs(result_mc) < 1e-10:
            # Both prices are effectively zero - test passes
            return result_mc, result_black
        elif abs(result_black) < 1e-10:
            # Black price is zero but MC isn't - check if MC is small enough
            self.assertLess(abs(result_mc), 1.0, 
                          f"MC price ({result_mc}) should be near zero when Black price is zero")
            return result_mc, result_black
        else:
            # Normal case - use relative tolerance
            self.assertAlmostEqual(result_mc/result_black, 1.0, delta=rel_tol,
                                 msg=f"MC price: {result_mc}, Black price: {result_black}")
            return result_mc, result_black

    def test_base_case(self):
        """Original test case kept for backward compatibility"""
        self._run_mc_vs_black_comparison(
            self.start_date, self.end_date, self.strike, self.notional,
            self.num_paths, self.a, self.sigma,
            self.curve_node_dates, self.curve_node_values
        )

    def test_time_horizons(self):
        """Test different time horizons"""
        horizons = [
            (1, 30),    # 1 month
            (1, 90),    # 3 months
            (1, 180),   # 6 months
            (1, 365),   # 1 year
            (1, 730),   # 2 years
            (90, 455),  # Starting in 3 months
            (180, 545), # Starting in 6 months
        ]

        for start, end in horizons:
            with self.subTest(f"Time horizon: {start} to {end} days"):
                self._run_mc_vs_black_comparison(
                    start, end, self.strike, self.notional,
                    self.num_paths, self.a, self.sigma,
                    [start, end], [1.0, 0.99]  # Simple curve for each horizon
                )

    def test_model_parameters(self):
        """Test different combinations of mean reversion and volatility"""
        params = [
            (0.001, 0.005),  # Low mean reversion, low vol
            (0.05, 0.005),   # High mean reversion, low vol
            (0.001, 0.02),   # Low mean reversion, high vol
            (0.05, 0.02),    # High mean reversion, high vol
        ]

        for a, sigma in params:
            with self.subTest(f"Parameters: a={a}, sigma={sigma}"):
                self._run_mc_vs_black_comparison(
                    self.start_date, self.end_date, self.strike, self.notional,
                    self.num_paths, a, sigma,
                    self.curve_node_dates, self.curve_node_values
                )

    def test_complex_curves(self):
        """Test various curve shapes with more node points"""
        # Upward sloping curve (all positive rates)
        up_dates = [1, 90, 180, 365, 730]
        # Approximately 1% continuous rate increasing to 2%
        up_values = [1.0, 0.9925, 0.985, 0.97, 0.955]  

        # Downward sloping curve (all positive rates)
        down_dates = [1, 90, 180, 365, 730]
        # Approximately 0.5% continuous rate decreasing to 0.25%
        down_values = [1.0, 0.9988, 0.9975, 0.995, 0.9925]  

        # Humped curve (all positive rates)
        hump_dates = [1, 90, 180, 365, 730]
        # Rate starts at 0.5%, peaks at 1%, ends at 0.75%
        hump_values = [1.0, 0.9988, 0.9965, 0.992, 0.985]  

        curves = [
            (up_dates, up_values, "Upward sloping"),
            (down_dates, down_values, "Downward sloping"),
            (hump_dates, hump_values, "Humped")
        ]

        for dates, values, curve_type in curves:
            with self.subTest(f"Curve type: {curve_type}"):
                self._run_mc_vs_black_comparison(
                    self.start_date, 730, self.strike, self.notional,
                    self.num_paths, self.a, self.sigma,
                    dates, values,
                    rel_tol=0.02  # Slightly wider tolerance for complex curves
                )

    def test_negative_rates(self):
        """Test that negative forward rates throw an exception"""
        # Curve that will generate negative forwards
        dates = [1, 365]
        values = [1.0, 1.1]  # This will create negative forwards
        
        with self.assertRaises(RuntimeError) as context:
            price_cap_black(
                self.start_date, self.end_date, self.strike, self.notional,
                self.num_paths, self.a, self.sigma,
                dates, values
            )
        
        self.assertTrue("Negative forward rates are not yet supported" in str(context.exception))

    def test_convergence(self):
        """Test convergence with increasing number of paths"""
        path_counts = [100, 250, 500, 1000]  # Reduced path counts for faster testing
        results = []
        rel_errors = []  # Track relative errors vs Black price
        
        # Use a curve with meaningful forward rates and a middle point
        test_dates = [1, 365, 730]  # 1Y and 2Y points
        test_values = [1.0, 0.98, 0.96]  # Approximately 2% continuous rate, consistent across curve
        
        # Use parameters that will generate meaningful option prices with more volatility
        test_strike = 0.02  # At-the-money strike (close to forward rate)
        test_sigma = 0.02   # Higher volatility to create more variance
        
        # First calculate the Black price (our target)
        black_price = price_cap_black(
            self.start_date, self.end_date, test_strike, self.notional,
            self.num_paths, self.a, test_sigma,
            test_dates, test_values
        )

        print(f"\nBlack price: {black_price}")  # Print reference price for debugging

        for num_paths in path_counts:
            with self.subTest(f"Number of paths: {num_paths}"):
                mc_price, _ = self._run_mc_vs_black_comparison(
                    self.start_date, self.end_date, test_strike, self.notional,
                    num_paths, self.a, test_sigma,
                    test_dates, test_values,
                    rel_tol=0.05  # Wider tolerance for lower path counts
                )
                results.append(mc_price)
                
                # Calculate relative error vs Black price
                rel_error = abs(mc_price - black_price) / abs(black_price)
                rel_errors.append(rel_error)
                
                print(f"Paths: {num_paths}, MC price: {mc_price}, Rel error: {rel_error}")  # Debug output
        
        # Check that variance between MC runs decreases
        variances = []
        for i in range(len(path_counts)-1):
            variance = abs(results[i] - results[i+1]) / abs(results[i])
            variances.append(variance)
            print(f"Variance between {path_counts[i]} and {path_counts[i+1]} paths: {variance}")  # Debug output
        
        # Verify decreasing variance between MC runs
        for i in range(len(variances)-1):
            self.assertGreaterEqual(variances[i], variances[i+1],
                                  "Variance between consecutive MC runs should decrease with more paths")
        
        # Verify convergence to Black price
        for i in range(len(rel_errors)-1):
            self.assertGreaterEqual(rel_errors[i], rel_errors[i+1],
                                  "Relative error vs Black price should decrease with more paths")
            
        # Check final convergence is within reasonable bounds
        self.assertLess(rel_errors[-1], 0.05,  # 5% relative error for reduced paths
                       f"Final MC price ({results[-1]}) should be within 5% of Black price ({black_price})")

    def test_strike_sensitivity(self):
        """Test different strike rates"""
        strikes = [0.01, 0.02, 0.03, 0.04, 0.05]
        
        for strike in strikes:
            with self.subTest(f"Strike: {strike}"):
                self._run_mc_vs_black_comparison(
                    self.start_date, self.end_date, strike, self.notional,
                    self.num_paths, self.a, self.sigma,
                    self.curve_node_dates, self.curve_node_values
                )

if __name__ == "__main__":
    unittest.main()
