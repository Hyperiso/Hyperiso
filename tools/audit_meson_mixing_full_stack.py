#!/usr/bin/env python3
"""Static audit of the complete neutral-meson-mixing implementation."""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def status(ok: bool, message: str) -> None:
    print(f"[{'PASS' if ok else 'FAIL'}] {message}")
    if not ok:
        raise AssertionError(message)


def main() -> None:
    m0 = read("Hyperiso/Hyperiso/core/src/BusinessLogic/domain/Decays/M0_Mixing.cpp")
    running = read("Hyperiso/Hyperiso/core/src/PhysicalModel/domain/WGroup/MesonMixingWilsonGroup.cpp")
    wids = read("Hyperiso/Hyperiso/core/src/Common/Mapper/wcoef_ids.hpp")
    enums = read("Hyperiso/Hyperiso/core/src/Common/GeneralEnum.h")
    maps = read("Hyperiso/Hyperiso/core/src/Common/Mapper/Map.cpp")
    bindings = read("Hyperiso/Hyperiso/core/src/Python/src/common/binding_common.cpp")
    py_enums = read("Hyperiso/Hyperiso/pyhyperiso/core/Common/GeneralEnum.py")
    deps = read("Hyperiso/Hyperiso/core/src/Common/DependenciesHelper.cpp")

    families = ("BD", "BS", "SD", "CU")
    coefficients = [
        f"{prefix}_{family}_{index}"
        for family in families
        for prefix, indices in (("C", range(1, 6)), ("CT", range(1, 4)))
        for index in indices
    ]
    status(all(name in wids for name in coefficients),
           "all 32 meson-mixing Wilson coefficients are in the PhysicalModel group")
    status(
        "getFR(WGroup::MESON_MIXING" in m0
        and "ContributionType::BSM" in m0,
        "M0_Mixing reads the running BSM coefficients from PhysicalModel",
    )
    status(
        all(token in running for token in (
            'is_B_family ? "UM_MATRIX_5" : "UM_MATRIX_4"',
            "MMRP::SUSY_to_BMU",
            "MMRP::BMU_to_SUSY",
            "/ (4. * PI)",
        )),
        "family-dependent QCD running and SUSY/BMU basis conversion are present",
    )

    observables = (
        "PHI_D", "DELTA_M_BD", "PHI_S", "DELTA_M_BS",
        "A_FS_D", "A_FS_S", "DELTA_M_K", "ABS_EPSILON_K",
        "X_D", "DELTA_M_D",
    )
    for name in observables:
        status(
            all(name in source for source in (enums, maps, bindings, py_enums, m0)),
            f"{name} is synchronized across enum, mapper, binding, Python and BusinessLogic",
        )

    expected_flha = {
        "PHI_D": "{511,   71,     1, -511}",
        "DELTA_M_BD": "{511,   7,      1, -511}",
        "PHI_S": "{531,   71,     1, -531}",
        "DELTA_M_BS": "{531,   7,      1, -531}",
        "A_FS_D": "{511,   74,     1, -511}",
        "A_FS_S": "{531,   74,     1, -531}",
        "DELTA_M_K": "{311,   7,      1, -311}",
        "ABS_EPSILON_K": "{311,   75,     1, -311}",
        "X_D": "{421,   72,     1, -421}",
        "DELTA_M_D": "{421,   7,      1, -421}",
    }
    for name, identifier in expected_flha.items():
        status(
            re.search(rf"\{{{name},\s*{re.escape(identifier)}\}}", maps) is not None,
            f"{name} has the expected FLHA identifier",
        )

    for pdg in (511, 531, 311):
        for index in range(1, 6):
            status(
                f'"FBAG", {{{pdg}, {index}}}' in deps,
                f"FBAG[{pdg},{index}] is propagated as a mixing nuisance",
            )

    assets = ROOT / "Hyperiso/Hyperiso/pyhyperiso/assets/default"
    obs = json.loads((assets / "observables.json").read_text())["FOBS"]["DEFAULT"]
    required_exp = {
        "511_7_0_0_0_0_0_0_1_-511",
        "531_7_0_0_0_0_0_0_1_-531",
        "311_75_0_0_0_0_0_0_1_-311",
        "511_74_0_0_0_0_0_0_1_-511",
        "531_74_0_0_0_0_0_0_1_-531",
    }
    status(required_exp <= set(obs), "safe default mixing likelihood is present")

    print("\n[WARN] PHI_D/PHI_S are arg(M12), not direct experimental interference phases.")
    print("[WARN] DELTA_M_K and D-mixing observables are not safe default likelihood entries.")
    print("[WARN] Gamma12 inputs and the ad-hoc B-mixing SM shifts still require physics review.")


if __name__ == "__main__":
    main()
