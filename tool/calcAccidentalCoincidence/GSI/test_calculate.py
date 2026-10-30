"""论文回归值、独立积分与数值边界检查。"""

import copy
import math
import unittest

from scipy.integrate import quad

from tool.calcAccidentalCoincidence.GSI.calculate import EXAMPLES, calculate, poisson_tail


class CoincidenceTests(unittest.TestCase):
    def test_repeated_category_uses_unique_total_rate(self):
        config = {"mode": "fixed", "rates_hz": [0.2, 0.3, 0.1],
                  "measurement_time_s": 1000, "sequence": [0, 1, 2, 2],
                  "windows_s": [0.5, 2, 2]}
        # 三个真实类别的下一事件概率密度独立积分，第三类发生两次。
        expected = 200
        for rate, dt in [(0.3, 0.5), (0.1, 2), (0.1, 2)]:
            expected *= rate / 0.6 * -math.expm1(-0.6 * dt)
        result = calculate(config)
        self.assertAlmostEqual(result["n_b"], expected)
        self.assertAlmostEqual(result["n_b_low_rate"], 1000 * 0.2 * 0.3 * 0.1**2 * 0.5 * 2 * 2)
        with self.assertRaises(ValueError):
            calculate({**config, "mode": "partial"})

    def test_paper_th214(self):
        for name, mean_low, probability, bound in [
            ("th214-5s", 0.875, 1.3879261702e-10, 1e-9),
            ("th214-40s", 7.0, 2.8451964606e-5, 1e-3),
        ]:
            with self.subTest(name=name):
                result = calculate(EXAMPLES[name])
                self.assertAlmostEqual(result["n_b_low_rate"], mean_low)
                self.assertAlmostEqual(result["P_err"] / probability, 1, places=9)
                self.assertLess(result["P_err_low_rate"], bound)

    def test_paper_element109(self):
        result = calculate(EXAMPLES["element109"])
        self.assertAlmostEqual(result["n_b_low_rate"] / 2.086933333333333e-18, 1)
        self.assertAlmostEqual(result["P_err"] / 2.0868390726e-18, 1, places=9)

    def test_two_event_density_integral_and_order_modes(self):
        config = {"rates_hz": [0.2, 0.3], "measurement_time_s": 1000,
                  "windows_s": [2], "observed_chains": 1}
        # 对原文概率密度独立数值积分。
        expected = 200 * quad(lambda t: 0.3 * math.exp(-0.5 * t), 0, 2)[0]
        for mode in ("fixed", "partial"):
            result = calculate({**config, "mode": mode})
            self.assertAlmostEqual(result["n_b"], expected)

    def test_position_bins_equal_sum_of_identical_regions(self):
        config = {"rates_hz": [0.2, 0.3, 0.1], "measurement_time_s": 1000,
                  "windows_s": [2, 5], "position_bins": 10}
        for mode in ("fixed", "partial"):
            total = calculate({**config, "mode": mode})
            single = calculate({**config, "mode": mode, "position_bins": 1,
                                "rates_hz": [0.02, 0.03, 0.01]})
            self.assertAlmostEqual(total["n_b"], 10 * single["n_b"])

    def test_counts_rates_and_measurement_exposure(self):
        config = copy.deepcopy(EXAMPLES["element109"])
        before = copy.deepcopy(config)
        first = calculate(config)
        self.assertEqual(config, before)
        del config["background_counts"], config["background_time_s"]
        config["rates_hz"] = first["rates_hz"]
        self.assertEqual(first["n_b"], calculate(config)["n_b"])
        config["measurement_time_s"] *= 2
        self.assertAlmostEqual(calculate(config)["n_b"] / first["n_b"], 2)

    def test_poisson_tail_and_underflow(self):
        for mean in (1e-20, 0.1, 5):
            probability, _ = poisson_tail(math.log(mean), 1)
            self.assertAlmostEqual(probability / -math.expm1(-mean), 1)
        probability, log_p = poisson_tail(math.log(1e-100), 4)
        self.assertEqual(probability, 0)
        self.assertAlmostEqual(log_p, 4 * math.log(1e-100) - math.log(24))
        # 即使期望数本身下溢，对数概率也应保留。
        result = calculate({"measurement_time_s": 1, "rates_hz": [1e-200, 1e-200],
                            "windows_s": [1], "observed_chains": 1})
        self.assertAlmostEqual(result["log10_P_err"], -400)

    def test_zero_and_observation_monotonicity(self):
        config = {"measurement_time_s": 100, "rates_hz": [0.1, 0.2],
                  "windows_s": [1]}
        probabilities = [calculate({**config, "observed_chains": n})["P_err"]
                         for n in range(5)]
        self.assertEqual(probabilities[0], 1)
        self.assertEqual(probabilities, sorted(probabilities, reverse=True))
        self.assertEqual(calculate({**config, "windows_s": [0]})["P_err"], 0)
        zero = calculate({**config, "rates_hz": [0, 0]})
        self.assertEqual(zero["P_err"], 0)
        self.assertTrue(zero["warnings"])

    def test_reject_invalid_input(self):
        base = EXAMPLES["th214-5s"]
        for override in [
            {"mode": "unknown"}, {"position_bins": 0}, {"position_bins": True},
            {"measurement_time_s": float("nan")}, {"measurement_time_s": -1},
            {"windows_s": [float("inf")]}, {"windows_s": [1, 2]},
            {"windows_s": [-1]}, {"observed_chains": 1.5}, {"typo": 1},
            {"background_counts": [1, 2.5]}, {"background_time_s": [1]},
            {"background_time_s": 0}, {"rates_hz": [0.1, 0.2]},
        ]:
            with self.subTest(override=override), self.assertRaises(ValueError):
                calculate({**base, **override})


if __name__ == "__main__":
    unittest.main()
