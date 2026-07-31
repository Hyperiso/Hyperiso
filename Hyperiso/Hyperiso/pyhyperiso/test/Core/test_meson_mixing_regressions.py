from __future__ import annotations

import runpy
from pathlib import Path


def test_meson_mixing_static_regressions() -> None:
    root = Path(__file__).resolve().parents[5]
    runpy.run_path(
        str(root / "tools" / "check_meson_mixing_regressions.py"),
        run_name="__main__",
    )
