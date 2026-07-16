"""Structural tests for spectrum-based THDM and SUSY reference cases."""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REPRO = ROOT / "reproducibility"


class BsmReferenceManifestTests(unittest.TestCase):
    """Ensure optional spectrum generators are not required by R6 and R7."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((REPRO / "manifest.json").read_text())
        cls.examples = {item["id"]: item for item in cls.manifest["examples"]}
        cls.metadata = json.loads(
            (REPRO / "expected_outputs" / "reference_metadata.json").read_text()
        )

    def test_spectrum_cases_are_present(self) -> None:
        """R6 and R7 must be part of the mandatory manifest."""
        self.assertIn("R6", self.examples)
        self.assertIn("R7", self.examples)

    def test_spectrum_cases_do_not_invoke_generators(self) -> None:
        """Archived spectra must be selected explicitly on the CLI."""
        for case_id in ("R6", "R7"):
            command = self.examples[case_id]["command"]
            self.assertIn("--spectrum true", command)
            self.assertIn("--contribution BSM", command)
            self.assertNotIn("softsusy", command.lower())
            self.assertNotIn("2hdmc", command.lower())
            self.assertNotIn("marty", command.lower())

    def test_archived_input_hashes_match_metadata(self) -> None:
        """Input files must match the hashes recorded with the references."""
        input_metadata = self.metadata.get("inputs")
        self.assertIsInstance(input_metadata, dict)
        expected_inputs = {
            path.name for path in (REPRO / "inputs").iterdir() if path.is_file()
        }
        self.assertEqual(set(input_metadata), expected_inputs)
        for filename, item in input_metadata.items():
            path = REPRO / "inputs" / filename
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(item["sha256"], digest)


if __name__ == "__main__":
    unittest.main()
