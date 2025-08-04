import unittest
from mc_pricer_py import price_cap_mc, price_cap_black

class TestCapPricing(unittest.TestCase):
    def setUp(self):
        # Use serial integer dates
        self.start_date = 1
        self.end_date = 366
        self.strike = 0.02
        self.notional = 1000000
        self.num_paths = 10
        self.a = 0.01
        self.sigma = 0.01
        self.curve_node_dates = [1, 366]
        self.curve_node_values = [1.0, 0.95]

    def test_price_cap_mc(self):
        result = price_cap_mc(
            self.start_date, self.end_date, self.strike, self.notional,
            self.num_paths, self.a, self.sigma,
            self.curve_node_dates, self.curve_node_values
        )
        self.assertIsInstance(result, float)

    def test_price_cap_black(self):
        result = price_cap_black(
            self.start_date, self.end_date, self.strike, self.notional,
            self.num_paths, self.a, self.sigma,
            self.curve_node_dates, self.curve_node_values
        )
        self.assertIsInstance(result, float)

if __name__ == "__main__":
    unittest.main()
