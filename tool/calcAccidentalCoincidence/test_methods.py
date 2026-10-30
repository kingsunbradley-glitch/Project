"""新增方法的文献例子、独立解析极限、搜索算法与统计边界检查。"""

import importlib.util
import itertools
import json
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2

from common import binomial_interval, poisson_log_tail


BASE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name.lower(), BASE / name / "calculate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


LAZAREV, BERKELEY, DUBNA, MC = [load(name) for name in ("Lazarev", "Berkeley", "Dubna", "MonteCarlo")]


class ProbabilityTests(unittest.TestCase):
    def test_lazarev_hand_calculation_and_two_event_limit(self):
        result = LAZAREV.calculate({"pair_counts": [12, 5, 2], "intermediate_counts": [1000, 800]})
        self.assertAlmostEqual(result["n_b"], 0.00015)
        self.assertAlmostEqual(result["P_at_least_observed_chains"], -math.expm1(-0.00015))
        self.assertAlmostEqual(LAZAREV.calculate({"pair_counts": [2.5]})["n_b"], 2.5)

    def test_berkeley_factorial_and_exact_versus_tail(self):
        config = {"alpha_rate_hz": 0.2, "window_s": 10, "position_bins": 2,
                  "n_alpha": 3, "implantations": 100}
        exact = BERKELEY.calculate(config)
        expected = math.exp(-1) / 6
        self.assertAlmostEqual(exact["P_per_implantation"], expected)
        self.assertAlmostEqual(exact["n_b"], 100 * expected)
        tail = BERKELEY.calculate({**config, "count_mode": "at_least"})
        self.assertAlmostEqual(tail["P_per_implantation"], 1 - math.exp(-1) * 2.5)
        self.assertGreater(tail["P_per_implantation"], expected)

    def test_berkeley_sequential_user_parameters_and_screenshot(self):
        config = {"model": "sequential", "step_background_counts": [65110, 154000, 154000],
                  "background_time_s": 655.2 * 3600, "position_bins": 6144,
                  "step_windows_s": [0.1, 100, 100], "implantations": 13250000}
        result = BERKELEY.calculate(config)
        self.assertAlmostEqual(result["n_b"] / 6.715226456307328e-6, 1)
        self.assertEqual(result["P_per_step"][1], result["P_per_step"][2])
        # 截图数值可由 750 h、alpha1 取整 65000、ER 取整 1.3e7 重建。
        figure = BERKELEY.calculate({**config, "background_time_s": 750 * 3600,
                                    "step_background_counts": [65000, 154000, 154000],
                                    "implantations": 13000000})
        self.assertAlmostEqual(figure["P_per_step"][1] / 9.2791e-4, 1, places=5)
        self.assertAlmostEqual(figure["P_per_implantation"] / 3.3737e-13, 1, places=5)
        self.assertAlmostEqual(figure["n_b"] / 4.38581e-6, 1, places=5)

    def test_berkeley_sequential_zero_rate_and_underflow(self):
        config = {"model": "sequential", "step_rates_hz": [1e-200, 1e-200],
                  "step_windows_s": [1e-200, 1], "implantations": 1}
        self.assertAlmostEqual(BERKELEY.calculate(config)["log10_P_per_implantation"], -600)
        zero = BERKELEY.calculate({**config, "step_rates_hz": [0, 1e-200]})
        self.assertEqual(zero["P_per_implantation"], 0)
        with self.assertRaises(ValueError):
            BERKELEY.calculate({**config, "step_windows_s": [1]})

    def test_dubna_published_counts(self):
        result = DUBNA.counts_test({"window_s": 34 * 60, "categories": [
            {"rate_hz": 1 / 3600, "count": 3}, {"rate_hz": 1.3 / 3600, "count": 1}]})
        self.assertAlmostEqual(result["P_exact_configuration"], 0.006068373369326093)

    def test_dubna_published_times_use_four_intervals(self):
        result = DUBNA.time_test({"background_rate_hz": 2.3 / 3600,
                                  "event_times_s": [0, 30, 954, 1050, 2040]})
        self.assertEqual(result["number_of_gaps"], 4)
        self.assertEqual(result["observed_mean_gap_s"], 510)
        self.assertAlmostEqual(result["background_expected_min_gap_s"] / 60, 6.52173913)
        self.assertAlmostEqual(result["background_expected_max_gap_s"] / 60, 54.34782608696)
        self.assertAlmostEqual(result["chain_fitted_expected_max_gap_s"] / 60, 17.7083333333)

    def test_decay_one_step_and_repeated_rate_erlang(self):
        half_life = math.log(2) / 0.2
        for time in (0, 1e-3, 2, 20):
            self.assertAlmostEqual(DUBNA.completion_cdf([half_life], time), -math.expm1(-0.2 * time))
            expected = 1 - math.exp(-0.2 * time) * (1 + 0.2 * time + (0.2 * time) ** 2 / 2)
            self.assertAlmostEqual(DUBNA.completion_cdf([half_life] * 3, time), expected)

    def test_decay_two_step_independent_numerical_convolution(self):
        time = 2.3
        expected = quad(lambda t: 0.2 * math.exp(-0.2 * t) * -math.expm1(-0.7 * (time - t)), 0, time)[0]
        self.assertAlmostEqual(DUBNA.completion_cdf([math.log(2) / 0.2, math.log(2) / 0.7], time), expected)

    def test_decay_bounds_enclose_interior_and_efficiency(self):
        result = DUBNA.decay_test({"time_s": 10, "half_life_bounds_s": [[1, 2], [3, 7]],
                                   "detection_efficiency": 0.8})
        for first in (1, 1.5, 2):
            for second in (3, 5, 7):
                probability = DUBNA.completion_cdf([first, second], 10)
                self.assertLessEqual(result["P_completed_min"], probability)
                self.assertGreaterEqual(result["P_completed_max"], probability)
        self.assertAlmostEqual(result["detected_probabilities"]["P_completed_max"], result["P_completed_max"] * 0.8)

    def test_position_correlations_and_identity_limit(self):
        config = {"standardized_differences": [1, 1, 1, 1],
                  "correlation_matrix": np.eye(4).tolist(), "samples": 100000, "seed": 17}
        independent = DUBNA.position_test(config)
        self.assertAlmostEqual(independent["variance_S"], 8)
        lower, upper = independent["tail_probability_CI95_MC"]
        self.assertLess(lower, chi2.sf(4, 4))
        self.assertGreater(upper, chi2.sf(4, 4))
        matrix = np.eye(4) + np.diag([-0.5] * 3, 1) + np.diag([-0.5] * 3, -1)
        correlated = DUBNA.position_test({**config, "correlation_matrix": matrix.tolist()})
        self.assertAlmostEqual(correlated["variance_S"], 11)
        self.assertNotEqual(correlated["central_reference_interval"], independent["central_reference_interval"])

    def test_poisson_extreme_tail_and_zero(self):
        self.assertEqual(poisson_log_tail(-math.inf, 0), 0)
        self.assertEqual(poisson_log_tail(-math.inf, 1), -math.inf)
        self.assertAlmostEqual(poisson_log_tail(-1000, 2), -2000 - math.log(2))
        result = BERKELEY.calculate({"alpha_rate_hz": 0, "window_s": 1, "n_alpha": 3, "implantations": 100})
        self.assertEqual(result["n_b"], 0)
        self.assertEqual(result["P_per_implantation"], 0)
        empty = DUBNA.counts_test({"window_s": 0, "categories": [{"rate_hz": 0, "count": 0}]})
        self.assertEqual(empty["P_exact_configuration"], 1)

    def test_binomial_zero_hits_has_upper_limit(self):
        lower, upper = binomial_interval(0, 1000)
        self.assertEqual(lower, 0)
        self.assertAlmostEqual(upper, -math.expm1(math.log(0.025) / 1000))
        result = MC.calculate({"mode": "counts", "inserted_events": 1000, "correlated_events": 0})
        self.assertGreater(result["P_per_insertion_interval"][1], 0)

    def test_search_matches_exhaustive_combinations(self):
        rng = np.random.default_rng(12)
        for anchor in range(3):
            for _ in range(50):
                lists = [np.unique(rng.integers(0, 12, 5)).astype(float) for _ in range(3)]
                anchor_time = float(rng.integers(0, 12))
                lists[anchor] = np.array([anchor_time])
                windows = [[0, 3], [1, 5]]
                oracle = any(all(t[i+1] > t[i] and lo <= t[i+1] - t[i] <= hi
                                 for i, (lo, hi) in enumerate(windows))
                             for t in itertools.product(*lists))
                self.assertEqual(MC.has_chain(lists, anchor_time, anchor, windows), oracle)

    def test_mc_demo_matches_analytic_acceptance(self):
        config = json.loads((BASE / "MonteCarlo/input.json").read_text())
        result = MC.calculate(config, BASE / "MonteCarlo")
        lower, upper = result["P_per_insertion_interval"]
        self.assertLess(lower, 0.075)
        self.assertGreater(upper, 0.075)

    def test_reject_invalid_inputs(self):
        cases = [
            (LAZAREV.calculate, {"pair_counts": [1, 2], "intermediate_counts": [0]}),
            (LAZAREV.calculate, {"pair_counts": [1, 2], "intermediate_counts": [2, 3]}),
            (BERKELEY.calculate, {"alpha_rate_hz": 1, "window_s": 1, "implantations": 2, "n_alpha": 0}),
            (DUBNA.time_test, {"background_rate_hz": 1, "event_times_s": [0, 2, 1]}),
            (DUBNA.decay_test, {"time_s": 1, "half_life_bounds_s": [[2, 1]]}),
            (DUBNA.position_test, {"standardized_differences": [1, 2], "correlation_matrix": [[1, 2], [2, 1]]}),
            (MC.calculate, {"mode": "counts", "inserted_events": 10, "correlated_events": 11}),
            (MC.calculate, {"mode": "counts", "inserted_events": 0, "correlated_events": 0}),
            (MC.calculate, {"mode": "counts", "inserted_events": 10, "correlated_events": 1, "seed": 1}),
        ]
        for function, config in cases:
            with self.subTest(function=function, config=config), self.assertRaises(ValueError):
                function(config)


if __name__ == "__main__":
    unittest.main()
