"""半衰期计算的解析、积分及命令行交叉验证；仅使用标准库。"""

import math
from pathlib import Path
import shutil
from statistics import NormalDist
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import calc_half_life as calc


class HalfLifeTests(unittest.TestCase):
    def test_single_event_analytic_interval(self):
        lower, upper = calc.definition_interval(1)
        self.assertAlmostEqual(lower, -1 / math.log(0.16), places=13)
        self.assertAlmostEqual(upper, -1 / math.log(0.84), places=12)
        n, tau, half, low, high = calc.calculate([10.0])
        self.assertEqual((n, tau), (1, 10.0))
        self.assertAlmostEqual(half, 10 * math.log(2))
        self.assertAlmostEqual(low / half, lower)
        self.assertAlmostEqual(high / half, upper)

    def test_paper_table_one_at_one_sigma(self):
        confidence = math.erf(1 / math.sqrt(2))
        for n, expected in [(1, (0.543, 5.79)), (2, (0.606, 2.82))]:
            low, high = calc.definition_interval(n, confidence)
            self.assertAlmostEqual(low, expected[0], delta=0.0005)
            self.assertAlmostEqual(high, expected[1], delta=0.005)

    def test_definition_coverage_by_independent_density_integration(self):
        # 对 Gamma 密度作 Simpson 积分，独立验证有限和求根得到的中心覆盖率。
        steps = 4000
        for n in range(1, 21):
            low, high = calc.definition_interval(n)
            a, b = n / high, n / low
            h = (b - a) / steps

            def density(x):
                return math.exp((n - 1) * math.log(x) - x - math.lgamma(n))

            integral = h / 3 * math.fsum(
                (1 if k in (0, steps) else (4 if k % 2 else 2)) * density(a + k * h)
                for k in range(steps + 1)
            )
            self.assertAlmostEqual(integral, 0.68, places=10)

    def test_likelihood_endpoints_using_direct_log_likelihood(self):
        threshold = NormalDist().inv_cdf(0.84) ** 2 / 2
        for n in (21, 22, 100, 10000):
            low, high = calc.likelihood_interval(n)
            self.assertLess(low, 1.0)
            self.assertGreater(high, 1.0)
            for tau in (low, high):
                log_likelihood = -n * math.log(tau) - n / tau
                self.assertAlmostEqual(-n - log_likelihood, threshold, places=10)

    def test_method_switch_at_21(self):
        for n in (20, 21):
            selected = "definition_interval" if n < 21 else "likelihood_interval"
            with patch.object(calc, selected, wraps=getattr(calc, selected)) as method:
                calc.calculate([1.0] * n)
                method.assert_called_once_with(n)

    def test_time_unit_scaling(self):
        for n in (1, 20, 21, 100):
            times = [float(i + 1) for i in range(n)]
            original = calc.calculate(times)
            scaled = calc.calculate([t * 1000 for t in times])
            for a, b in zip(original[1:], scaled[1:]):
                self.assertAlmostEqual(b / a, 1000)

    def test_invalid_times(self):
        for times in ([], [0], [-1], [math.nan], [math.inf]):
            with self.subTest(times=times), self.assertRaises(ValueError):
                calc.calculate(times)

    def test_input_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.dat"
            path.write_text("\ufeff# comment\n1.2D+1 # event\n\n3e-1\n", encoding="utf-8")
            self.assertEqual(calc.read_times(path), [12.0, 0.3])
            for text in ("", "# empty", "1 2", "nan", "inf", "0", "-1", "bad"):
                path.write_text(text, encoding="utf-8")
                with self.subTest(text=text), self.assertRaises(ValueError):
                    calc.read_times(path)

    def test_cli_uses_script_directory_and_reports_selected_method(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scripts = root / "scripts"
            scripts.mkdir()
            script = scripts / "calc_half_life.py"
            shutil.copyfile(calc.__file__, script)
            for n, expected in ((20, "论文定义法"), (21, "极大似然估计")):
                (scripts / "input.dat").write_text("1\n" * n, encoding="utf-8")
                result = subprocess.run(
                    [sys.executable, str(script), "--unit", "ms"],
                    cwd=root, capture_output=True, text=True, check=True,
                )
                self.assertIn(expected, result.stdout)
                self.assertIn(f"N = {n}", result.stdout)
                self.assertIn("C.L. = 0.68", result.stdout)
                self.assertIn("ms", result.stdout)


if __name__ == "__main__":
    unittest.main()
