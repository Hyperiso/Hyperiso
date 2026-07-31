#!/usr/bin/env python3
"""Static non-regression checks for the neutral-meson-mixing pipeline.

This checker intentionally avoids importing the compiled extension. It protects
physics/convention fixes that are easy to regress during refactors:

* synchronized observable exposure (C++, maps, bindings, Python and GUI);
* M0_Mix parameter numbering and family-specific low scales;
* chiral enhancement of Q2-Q5;
* alpha_s/(4*pi), U5/U4 selection and BMU->SUSY return;
* complex FWCOEF + i IMFWCOEF composition;
* tree-level Z' normalization of C1, CT1 and C5 in all four families.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_observables() -> None:
    cpp_enum = text("Hyperiso/Hyperiso/core/src/Common/GeneralEnum.h")
    cpp_map = text("Hyperiso/Hyperiso/core/src/Common/Mapper/Map.cpp")
    binding = text("Hyperiso/Hyperiso/core/src/Python/src/common/binding_common.cpp")
    py_enum = text("Hyperiso/Hyperiso/pyhyperiso/core/Common/GeneralEnum.py")
    gui = text("GHyperiso/HyperisoDashGUI/pyhyperiso_dash/latex_data/observable_flha_to_latex_map.py")
    m0 = text("Hyperiso/Hyperiso/core/src/BusinessLogic/domain/Decays/M0_Mixing.cpp")

    for name in ("A_FS_D", "A_FS_S", "DELTA_M_D"):
        for label, source in (
            ("C++ enum", cpp_enum), ("C++ map", cpp_map),
            ("pybind", binding), ("Python enum", py_enum),
            ("GUI map", gui), ("M0_Mixing", m0),
        ):
            require(name in source, f"{name} missing from {label}")

    require("A_FS = A_FS_S" in cpp_enum, "C++ A_FS compatibility alias missing")
    require('"A_FS": "A_FS_S"' in gui, "GUI A_FS compatibility alias missing")


def check_parameters_and_m0() -> None:
    params = json.loads(text("Hyperiso/Hyperiso/pyhyperiso/assets/default/parameters.json"))["M0_Mix"]
    expected = {
        "9": 0.9154, "10": 4.16, "11": 3.0, "12": 3.0,
        "13": 1.87, "14": 0.496, "15": 5.293e9,
    }
    for key, value in expected.items():
        got = float(params[key]["central_value"])
        require(abs(got - value) <= 1e-12 * max(1.0, abs(value)),
                f"M0_Mix[{key}]={got}, expected {value}")
    for key in ("7_1", "7_2", "8"):
        require(key in params, f"M0_Mix[{key}] missing")

    m0 = text("Hyperiso/Hyperiso/core/src/BusinessLogic/domain/Decays/M0_Mixing.cpp")
    require("if (i > 0)" in m0, "Q2 chiral enhancement regression: expected i > 0")
    require('"SCALE_NUIS", 1' in m0,
            "M0_Mixing must derive its SM electroweak scale from SCALE_NUIS")
    require('cache.mu_W = (*p)(ParamId{ParameterType::WILSON, "EW_SCALE", 1}' not in m0,
            "SM mixing must not reuse the BSM EW_SCALE matching scale")
    require("std::pow(2.0, x_W) * cache.m_W" in m0,
            "SM mixing scale must remain tied to m_W")
    require('"B_SCALE", 1}, mu_B' in m0, "B_SCALE not installed before coefficient reads")
    require('"K_SCALE", 1}, mu_K' in m0, "K_SCALE not installed before coefficient reads")
    require('"D_SCALE", 1}, mu_D' in m0, "D_SCALE not installed before coefficient reads")
    pos_d = m0.index('"D_SCALE", 1}, mu_D')
    pos_k_read = m0.index("populate_C(cache.C_K")
    require(pos_d < pos_k_read, "D_SCALE must be installed before SD running is read")
    for idx, name in ((10,"mu_B"),(11,"mu_K"),(12,"mu_D"),(13,"eta_cc"),(14,"eta_ct"),(15,"delta_MK_exp")):
        require(f'"M0_Mix", {idx}' in m0, f"M0_Mix index {idx} for {name} missing")


def check_running() -> None:
    source = text("Hyperiso/Hyperiso/core/src/PhysicalModel/domain/WGroup/MesonMixingWilsonGroup.cpp")
    require('is_B_family ? "UM_MATRIX_5" : "UM_MATRIX_4"' in source,
            "family-specific U5/U4 selection missing")
    require('/ (4. * PI)' in source, "NLO prefactor must be alpha_s/(4*pi)")
    require('n == 2 ? 6 : 5' in source, "K/D alpha_s low-scale selection missing")
    require(source.count("MMRP::BMU_to_SUSY") >= 2,
            "BMU->SUSY matching/running return missing")


def check_complex_input() -> None:
    source = text("Hyperiso/Hyperiso/core/src/PhysicalModel/domain/WilsonManager.cpp")
    require("WilsonBlockNames::imfwcoef()" in source, "IMFWCOEF block is not read")
    require("scalar_t(real_part, imag_part)" in source,
            "FWCOEF and IMFWCOEF are not combined into a complex coefficient")


def check_templates() -> None:
    base = ROOT / "Hyperiso/Hyperiso/pyhyperiso/assets/template/MARTY"
    for family in ("BD", "BS", "SD", "CU"):
        c1 = (base / f"C_{family}_1.cpp").read_text(encoding="utf-8")
        ct1 = (base / f"CT_{family}_1.cpp").read_text(encoding="utf-8")
        c5 = (base / f"C_{family}_5.cpp").read_text(encoding="utf-8")
        require("DiracCoupling::VL, mty::DiracCoupling::VL" in c1,
                f"C_{family}_1 is not VLL")
        require(re.search(r"Expr C\s*=\s*2\s*\*", c1) is not None,
                f"C_{family}_1 identical-current factor 2 missing")
        require("DiracCoupling::VR, mty::DiracCoupling::VR" in ct1,
                f"CT_{family}_1 is not VRR")
        require(re.search(r"Expr C\s*=\s*2\s*\*", ct1) is not None,
                f"CT_{family}_1 identical-current factor 2 missing")
        require("DiracCoupling::VL" in c5 and "DiracCoupling::VR" in c5,
                f"C_{family}_5 does not project vector LR")
        require(re.search(r"Expr C\s*=\s*-8\s*\*", c5) is not None,
                f"C_{family}_5 BMU->SUSY normalization -8 missing")


def main() -> None:
    check_observables()
    check_parameters_and_m0()
    check_running()
    check_complex_input()
    check_templates()
    print("[OK] meson-mixing observables, running, complex input and MARTY templates are synchronized.")


if __name__ == "__main__":
    main()
