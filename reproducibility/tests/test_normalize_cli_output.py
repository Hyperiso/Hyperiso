from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "scripts" / "normalize_cli_output.py"
)
SPEC = importlib.util.spec_from_file_location("normalize_cli_output", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
normalizer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(normalizer)


class NormalizeCliOutputTests(unittest.TestCase):
    def test_diagnostic_error_includes_the_original_line(self) -> None:
        text = (
            "[2026-07-15_01] [WARN] LHA reader: Unknown block THDM "
            "encountered. Skipping.\n"
            "\n== Wilson summary ==\norder=NNLO\n"
        )
        with self.assertRaisesRegex(ValueError, "Unknown block THDM"):
            normalizer.normalize(text)

    def test_clean_summary_is_preserved(self) -> None:
        text = "startup text\n\n== Wilson summary ==\norder=NNLO\n"
        self.assertEqual(
            normalizer.normalize(text),
            "== Wilson summary ==\norder=NNLO\n",
        )


if __name__ == "__main__":
    unittest.main()
