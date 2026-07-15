from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "check_expected_outputs.py"
)
SPEC = importlib.util.spec_from_file_location("check_expected_outputs", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checker)


class CsvReferenceComparisonTests(unittest.TestCase):
    def test_serialization_tolerance_matches_recorded_precision(self) -> None:
        self.assertAlmostEqual(
            checker.serialization_abs_tol("0.000313876"), 5e-10, places=15
        )
        self.assertAlmostEqual(
            checker.serialization_abs_tol("2.78341e-09"), 5e-15, places=20
        )

    def test_full_precision_csv_matches_with_declared_numeric_tolerance(self) -> None:
        spec = {
            "columns": ["BR_B__Xs_gamma [0, 0]", "BR_Bs__mu_mu [0, 0]"],
            "rows": 1,
            "minimum": 0.0,
            "maximum": 1.0,
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.csv"
            output = root / "output.csv"
            for path, row in (
                (reference, ["0.00031387579123833616", "2.7834103000000001e-09"]),
                (output, ["0.00031387579123833620", "2.7834103000000005e-09"]),
            ):
                with path.open("w", newline="") as handle:
                    writer = csv.writer(handle, lineterminator="\n")
                    writer.writerow(spec["columns"])
                    writer.writerow(row)

            errors = checker.compare_csv(reference, output, spec, 0.0, 1e-12)
            self.assertEqual(errors, [])

    def test_reference_outside_rounding_interval_fails(self) -> None:
        spec = {"columns": ["x"], "rows": 1}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.csv"
            output = root / "output.csv"
            reference.write_text("x\n0.000313876\n")
            output.write_text("x\n0.000313874\n")
            errors = checker.compare_csv(reference, output, spec, 0.0, 1e-8)
            self.assertTrue(errors)


class FrozenReferenceTests(unittest.TestCase):
    def test_r5_reference_has_200_rows(self) -> None:
        root = Path(__file__).resolve().parents[2]
        csv_path = (
            root / "reproducibility" / "expected_outputs" / "statistics_samples.csv"
        )
        with csv_path.open(newline="") as handle:
            rows = list(csv.reader(handle))
        self.assertEqual(len(rows) - 1, 200)


if __name__ == "__main__":
    unittest.main()
