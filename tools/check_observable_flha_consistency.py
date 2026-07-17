#!/usr/bin/env python3
"""Validate duplicated observable FLHA maps and project-defined type ids."""

from __future__ import annotations

import importlib.util
from collections import defaultdict
from pathlib import Path
from types import ModuleType
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
INPUT_MAP_PATH = ROOT / "Hyperiso/InputHelper/observable_maps.py"
GUI_MAP_PATH = (
    ROOT
    / "GHyperiso/HyperisoDashGUI/pyhyperiso_dash/latex_data/observable_flha_to_latex_map.py"
)

EXPECTED_IDS = {
    "F_T_B0__KSTAR0_MU_MU": (511, 932, 3, 313, 13, -13),
    "ALPHA_K_B0__KSTAR0_MU_MU": (511, 933, 3, 313, 13, -13),
    "P_TAU_B0__D_TAU_NU": (511, 9212, 3, 411, -15, 16),
    "P_TAU_B__DSTAR0_TAU_NU": (521, 9212, 3, 423, -15, 16),
    "P_D_B0__DSTAR_TAU_NU": (511, 9221, 3, 413, -15, 16),
    "P_D_B__DSTAR0_TAU_NU": (521, 9221, 3, 423, -15, 16),
}


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalized(mapping: dict[str, Any]) -> dict[str, tuple[int, ...]]:
    return {name: tuple(int(part) for part in value) for name, value in mapping.items()}


def duplicate_ids(
    mapping: dict[str, tuple[int, ...]],
) -> dict[tuple[int, ...], list[str]]:
    reverse: defaultdict[tuple[int, ...], list[str]] = defaultdict(list)
    for name, flha_id in mapping.items():
        reverse[flha_id].append(name)
    return {flha_id: names for flha_id, names in reverse.items() if len(names) > 1}


def main() -> int:
    input_module = load_module(INPUT_MAP_PATH, "hyperiso_input_observable_maps")
    gui_module = load_module(GUI_MAP_PATH, "hyperiso_gui_observable_maps")

    input_map = normalized(input_module.OBSERVABLE_FLHA_MAPPING)
    gui_map = normalized(gui_module.OBSERVABLE_ENUM_TO_FLHA_MAP)

    errors: list[str] = []

    duplicates = duplicate_ids(input_map)
    if duplicates:
        errors.append(f"InputHelper contains duplicate FLHA ids: {duplicates}")

    if input_map != gui_map:
        missing_gui = sorted(input_map.keys() - gui_map.keys())
        missing_input = sorted(gui_map.keys() - input_map.keys())
        mismatched = sorted(
            name
            for name in input_map.keys() & gui_map.keys()
            if input_map[name] != gui_map[name]
        )
        errors.append(
            "GUI/InputHelper FLHA maps differ: "
            f"missing_gui={missing_gui}, missing_input={missing_input}, mismatched={mismatched}"
        )

    for name, expected in EXPECTED_IDS.items():
        actual = input_map.get(name)
        if actual != expected:
            errors.append(f"{name}: expected {expected}, got {actual}")

    if any(92015 in flha_id or 921423 in flha_id for flha_id in input_map.values()):
        errors.append("Legacy polarization type ids 92015/921423 are still present")

    if getattr(gui_module, "OBSERVABLE_COLLIDING_FLHA_IDS", None):
        errors.append(
            "GUI generated map still reports FLHA collisions: "
            f"{gui_module.OBSERVABLE_COLLIDING_FLHA_IDS}"
        )

    legacy_tau = (511, 92015, 3, 411, -15, 16)
    canonical_tau = EXPECTED_IDS["P_TAU_B0__D_TAU_NU"]
    if input_module.canonicalize_observable_flha_id(legacy_tau) != canonical_tau:
        errors.append("InputHelper does not canonicalize legacy fermion polarization")
    if gui_module.get_observable_enum_names(legacy_tau) != ("P_TAU_B0__D_TAU_NU",):
        errors.append("GUI lookup does not accept legacy fermion polarization")

    legacy_dstar = (511, 921423, 3, 413, -15, 16)
    canonical_dstar = EXPECTED_IDS["P_D_B0__DSTAR_TAU_NU"]
    if input_module.canonicalize_observable_flha_id(legacy_dstar) != canonical_dstar:
        errors.append("InputHelper does not canonicalize legacy vector polarization")
    if gui_module.get_observable_enum_names(legacy_dstar) != ("P_D_B0__DSTAR_TAU_NU",):
        errors.append("GUI lookup does not accept legacy vector polarization")

    if errors:
        for error in errors:
            print(f"[FAIL] {error}")
        return 1

    print(
        f"[OK] {len(input_map)} observable FLHA ids are unique and synchronized; "
        "932/933 and 92ij conventions are consistent."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
