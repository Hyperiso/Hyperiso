#!/usr/bin/env python3
"""Audit neutral-meson-mixing observable exposure across all project layers."""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OBS = {
    "A_FS_D": (511, 74, 1, -511),
    "A_FS_S": (531, 74, 1, -531),
    "DELTA_M_D": (421, 7, 1, -421),
}


def text(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def load(rel: str, name: str):
    path = ROOT / rel
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    errors: list[str] = []
    cpp_enum = text("Hyperiso/Hyperiso/core/src/Common/GeneralEnum.h")
    mapper = text("Hyperiso/Hyperiso/core/src/Common/Mapper/Map.cpp")
    binding = text("Hyperiso/Hyperiso/core/src/Python/src/common/binding_common.cpp")
    wrapper = text("Hyperiso/Hyperiso/pyhyperiso/core/Common/GeneralEnum.py")
    mixing = text("Hyperiso/Hyperiso/core/src/BusinessLogic/domain/Decays/M0_Mixing.cpp")

    for name in OBS:
        for label, content in (
            ("C++ enum", cpp_enum),
            ("mapper", mapper),
            ("pybind", binding),
            ("Python enum", wrapper),
            ("M0_Mixing switch", mixing),
        ):
            if name not in content:
                errors.append(f"{name} missing from {label}")

    input_maps = load("Hyperiso/InputHelper/observable_maps.py", "mix_input_maps")
    gui_maps = load(
        "GHyperiso/HyperisoDashGUI/pyhyperiso_dash/latex_data/observable_flha_to_latex_map.py",
        "mix_gui_maps",
    )
    for name, expected in OBS.items():
        got_input = tuple(input_maps.OBSERVABLE_FLHA_MAPPING.get(name, ()))
        got_gui = tuple(gui_maps.OBSERVABLE_ENUM_TO_FLHA_MAP.get(name, ()))
        if got_input != expected:
            errors.append(f"InputHelper {name}: expected {expected}, got {got_input}")
        if got_gui != expected:
            errors.append(f"GUI {name}: expected {expected}, got {got_gui}")
        if not gui_maps.OBSERVABLE_ENUM_TO_LATEX_MAP.get(name):
            errors.append(f"GUI LaTeX missing for {name}")

    params = json.loads(text("Hyperiso/Hyperiso/pyhyperiso/assets/default/parameters.json"))["M0_Mix"]
    required_params = {
        "7_1", "7_2", "8", "9", "10", "11", "12", "13", "14", "15"
    }
    missing = required_params - set(params)
    if missing:
        errors.append(f"M0_Mix default parameters missing: {sorted(missing)}")

    # The family scales must follow the new observable-parameter numbering.
    for index in (10, 11, 12):
        needle = f'"M0_Mix", {index}'
        if needle not in mixing:
            errors.append(f"M0_Mixing scale index {index} is not used")

    if errors:
        for error in errors:
            print(f"[FAIL] {error}")
        return 1

    print("[OK] A_FS_D, A_FS_S and DELTA_M_D are synchronized across C++, maps, bindings, Python, defaults and GUI.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
