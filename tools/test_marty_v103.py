#!/usr/bin/env python3
"""Validation suite for the HyperIso 1.0.3 MARTY patch.

The default mode is safe and dependency-light: it checks repository assets,
release metadata, the Z-prime mapping, critical C9/C10 template ABI markers,
SM-path portability, and the source-level SM/BSM/TOTAL/cache fixes.

Optional modes:

* ``--ctest-build-dir BUILD`` runs C++ MARTY/Wilson tests from an existing build.
* ``--integration --marty-install PATH`` runs a real Python/MARTY model,
  validates finite C1..C10 values and TOTAL = SM + BSM, checks that C2 is not
  doubled, compares the pure BSM result across both SM providers, and repeats
  the build to detect analytical cache rewrites.
* ``--thdm-type 2`` rewrites MINPAR(24) in a temporary LHA file for a THDM-II
  regression run.
* ``--zprime-mass-scan ...`` checks that a selected BSM coefficient decouples
  at the heavy-mass endpoint (C9 by default).

The integration mode intentionally runs in child processes because HyperIso owns
process-global runtime state. This also tests cache reuse across Python processes.
"""

from __future__ import annotations

import argparse
import compileall
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = ROOT / "Hyperiso" / "Hyperiso"
PYTHON_PACKAGE = PACKAGE_ROOT / "pyhyperiso"
PACKAGED_ASSETS = PYTHON_PACKAGE / "assets"
CANONICAL_TEMPLATES = PACKAGED_ASSETS / "template" / "MARTY"
MIRROR_TEMPLATES = ROOT / "Assets" / "template" / "MARTY"
CHILD_RESULT_PREFIX = "HYPERISO_MARTY_TEST_RESULT="


@dataclass
class CheckResult:
    name: str
    status: str
    detail: str = ""
    duration_s: float = 0.0


class ValidationError(RuntimeError):
    """Raised when a validation step fails."""


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValidationError(f"cannot read {path}: {exc}") from exc


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def _require_contains(path: Path, needles: Iterable[str]) -> None:
    text = _read(path)
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise ValidationError(f"{path}: missing expected markers: {missing}")


def _run_command(
    command: Sequence[str],
    *,
    cwd: Path = ROOT,
    env: dict[str, str] | None = None,
    timeout: float | None = None,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    if completed.returncode != 0:
        output = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
        raise ValidationError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{output.strip()}"
        )
    return completed


def check_repository_layout() -> str:
    required = [
        PACKAGE_ROOT / "pyproject.toml",
        CANONICAL_TEMPLATES / "C9.cpp",
        CANONICAL_TEMPLATES / "C10.cpp",
        PACKAGED_ASSETS / "input_files" / "marty_model" / "sm.h",
        PACKAGED_ASSETS / "input_files" / "marty_model" / "ZPrime.h",
        PACKAGED_ASSETS / "input_files" / "marty_mapping" / "zprime.json",
        ROOT / "tools" / "check_marty_template_sync.py",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    _require(not missing, f"missing required repository files: {missing}")
    return f"{len(required)} required files found"


def check_version_consistency() -> str:
    tool = ROOT / "tools" / "check_version_consistency.py"
    completed = _run_command([sys.executable, str(tool), "--root", str(ROOT)])
    output = completed.stdout.strip()
    _require("1.0.3" in output, f"expected release 1.0.3, got: {output}")
    return output


def check_template_sync() -> str:
    tool = ROOT / "tools" / "check_marty_template_sync.py"
    completed = _run_command([sys.executable, str(tool)])
    return completed.stdout.strip()


def check_template_abi() -> str:
    expected = {
        "C9.cpp": (
            "HYPERISO_MARTY_OPERATOR_NORM_ABI",
            "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-c9-tree-first-split-regprop",
            "HyperisoMartyC9LinkerSelection::NonPhotonVector",
            "HyperisoMartyC9LinkerSelection::PhotonOnly",
            "DiagramParticleType::External",
            "DiagramParticleType::Mediator",
        ),
        "C10.cpp": (
            "HYPERISO_MARTY_OPERATOR_NORM_ABI",
            "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-c10",
        ),
        "CP9.cpp": (
            "HYPERISO_MARTY_OPERATOR_NORM_ABI",
            "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-cp9-tree-first-split-regprop",
            "DiagramParticleType::External",
            "DiagramParticleType::Mediator",
        ),
        "CP10.cpp": (
            "HYPERISO_MARTY_OPERATOR_NORM_ABI",
            "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-cp10-tree-first-split-regprop",
            "DiagramParticleType::External",
            "DiagramParticleType::Mediator",
        ),
    }
    for name, markers in expected.items():
        _require_contains(CANONICAL_TEMPLATES / name, markers)
    return f"validated ABI markers in {', '.join(expected)}"


def check_zprime_mapping() -> str:
    path = PACKAGED_ASSETS / "input_files" / "marty_mapping" / "zprime.json"
    try:
        mapping = json.loads(_read(path))
    except json.JSONDecodeError as exc:
        raise ValidationError(f"invalid JSON in {path}: {exc}") from exc

    for parameter in ("m_X", "m_Z_X"):
        _require(parameter in mapping, f"{path}: missing {parameter}")
        entry = mapping[parameter]
        _require(entry.get("block") == "MASS", f"{parameter} must map to MASS")
        _require(str(entry.get("pdgCode")) == "32", f"{parameter} must map to PDG 32")
    return "m_X and m_Z_X both map to MASS(32)"


def check_packaged_model_portability() -> str:
    model_dir = PACKAGED_ASSETS / "input_files" / "marty_model"
    offending: list[str] = []
    forbidden = ("/project/Assets", "Third_party/MARTY/MARTY_INSTALL/include")
    for path in sorted(model_dir.glob("*.h")):
        text = _read(path)
        for token in forbidden:
            if token in text:
                offending.append(f"{path.relative_to(ROOT)} contains {token!r}")
    _require(not offending, "non-portable packaged model paths found: " + "; ".join(offending))
    return f"checked {len(list(model_dir.glob('*.h')))} packaged model headers"


def check_sm_path_resolution() -> str:
    path = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "Core"
        / "adapters"
        / "MartyAdapter.cpp"
    )
    _require_contains(
        path,
        (
            "case MartyPath::SM_MODEL_FILE",
            "APIPath::ASSETS_ROOT",
            '"input_files" / "marty_model" / "sm.h"',
            "Packaged MARTY SM model header is missing",
        ),
    )
    _require("/project/Assets" not in _read(path), f"{path} still embeds /project/Assets")
    return "SM model resolves from the runtime packaged assets root"


def check_contribution_composition_sources() -> str:
    registry = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "PhysicalModel"
        / "domain"
        / "WilsonCoefficientRegistry.cpp"
    )
    builder = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "PhysicalModel"
        / "adapters"
        / "WilsonBuilder.cpp"
    )
    manager = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "PhysicalModel"
        / "domain"
        / "WilsonManager.cpp"
    )
    _require_contains(
        registry,
        (
            "ctx.contrib == ContributionType::BSM",
            "bsm_split_generation = true",
            "ctx.backend == Backend::Marty",
        ),
    )
    builder_text = _read(builder)
    _require(
        builder_text.count("? ContributionType::SM\n            : ContributionType::BSM;") >= 2,
        "WilsonBuilder must store non-SM MARTY calculations as pure BSM in build() and add()",
    )
    _require(
        "ContributionType::TOTAL : ContributionType::BSM" not in builder_text,
        "WilsonBuilder still selects TOTAL for the MARTY target model",
    )
    _require_contains(
        manager,
        (
            "src.get_val(pid_sm) + src.get_val(pid_bsm)",
            "calculation_is_total",
        ),
    )
    return "normal MARTY path stores pure BSM and composes TOTAL = SM + BSM"


def check_tree_first_policy() -> str:
    modifier = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "ExternalIntegration"
        / "MartyInterface"
        / "GeneralModelModifier.cpp"
    )
    c10 = CANONICAL_TEMPLATES / "C10.cpp"
    _require_contains(
        modifier,
        (
            "HYPERISO_MARTY_BSM_SPLIT_ABI: model-split-v26",
            "hyperiso_marty_require_non_sm_diagram_particle",
            "DiagramParticleType::External",
            "mty::Order::TreeLevel",
            'trimmed == "mty::Order::TreeLevel,"',
            'line.replace(first, last - first + 1, "hyperiso_marty_order,")',
            "hyperiso_marty_use_tree_level",
            "if (!hyperiso_marty_use_tree_level)",
            "mty::Order::OneLoop",
            "selected order=",
        ),
    )
    _require_contains(
        c10,
        (
            "semileptonic-c10-tree-first-full-4f-v9",
            "if (C10_tree == CSL_0)",
            "mty::Order::OneLoop",
            "[MARTY C10] selected order=",
        ),
    )
    return "C9/CP9/CP10 and C10 use tree-first matching with one-loop fallback"


def check_cache_source_signatures() -> str:
    interface = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "ExternalIntegration"
        / "MartyInterface"
        / "MartyInterface.cpp"
    )
    modifier = (
        ROOT
        / "Hyperiso"
        / "Hyperiso"
        / "core"
        / "src"
        / "ExternalIntegration"
        / "MartyInterface"
        / "GeneralModelModifier.cpp"
    )
    _require_contains(
        interface,
        (
            "HYPERISO_MARTY_CACHE_ABI: pyhyperiso-1.0.3-v2",
            "HYPERISO_MARTY_GENERATION_MODE",
            "has_cache_abi = has_cache_abi ||",
            "has_model_signature",
            "has_template_signature",
            "has_generation_mode",
            "while (std::getline(in, line))",
        ),
    )
    _require_contains(modifier, ("FNV-1a", "modelSignature"))
    return "cache key covers ABI, full model content, template content and generation mode"


def check_python_compilation() -> str:
    with tempfile.TemporaryDirectory(prefix="hyperiso-pyc-") as pycache:
        old_prefix = sys.pycache_prefix
        sys.pycache_prefix = pycache
        try:
            targets = [ROOT / "tools", PYTHON_PACKAGE]
            ok = all(compileall.compile_dir(str(path), quiet=1, force=True) for path in targets)
        finally:
            sys.pycache_prefix = old_prefix
    _require(ok, "Python byte-compilation failed")
    return "tools and pyhyperiso Python sources compile"


def check_import_smoke() -> str:
    env = os.environ.copy()
    current = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(PACKAGE_ROOT) + (os.pathsep + current if current else "")
    code = """
import pyhyperiso
from importlib.resources import files
p = files('pyhyperiso').joinpath('assets', 'template', 'MARTY', 'C9.cpp')
assert p.is_file()
print(pyhyperiso.__version__)
print(p)
"""
    try:
        completed = _run_command([sys.executable, "-c", code], env=env, timeout=60)
    except ValidationError as exc:
        message = str(exc)
        if "phyperiso" in message or "ModuleNotFoundError" in message or "ImportError" in message:
            raise FileNotFoundError(
                "native pyhyperiso extension is not installed/built; import smoke test skipped"
            ) from exc
        raise
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    _require(lines and lines[0] == "1.0.3", f"unexpected imported version: {lines}")
    return f"installed/importable pyhyperiso {lines[0]} uses packaged C9.cpp"


def run_ctest(build_dir: Path) -> str:
    build_dir = build_dir.resolve()
    _require(build_dir.is_dir(), f"CMake build directory does not exist: {build_dir}")
    _require(shutil.which("ctest") is not None, "ctest is not available")
    completed = _run_command(
        [
            "ctest",
            "--test-dir",
            str(build_dir),
            "--output-on-failure",
            "-R",
            "Marty|Wilson",
        ],
        timeout=900,
    )
    summary = next(
        (line.strip() for line in reversed(completed.stdout.splitlines()) if "tests passed" in line),
        "ctest MARTY/Wilson selection passed",
    )
    return summary



_LHA_BLOCK_RE = re.compile(r"^\s*BLOCK\s+(\S+)", re.IGNORECASE)


def _parse_float_list(raw: str) -> list[float]:
    values: list[float] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        value = float(token)
        _require(math.isfinite(value) and value > 0.0, f"invalid positive scan value: {token}")
        values.append(value)
    _require(len(values) >= 2, "a decoupling scan needs at least two positive values")
    return values


def _rewrite_lha_scalar(source: Path, destination: Path, block: str, code: int, value: float) -> None:
    """Rewrite a one-index Les Houches entry while preserving the rest of the card."""
    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    current_block = ""
    replaced = False
    rendered = f"{value:.16e}"

    for index, line in enumerate(lines):
        match = _LHA_BLOCK_RE.match(line)
        if match:
            current_block = match.group(1).upper()
            continue
        if current_block != block.upper():
            continue

        body, marker, comment = line.partition("#")
        tokens = body.split()
        if len(tokens) < 2:
            continue
        try:
            entry_code = int(tokens[0])
        except ValueError:
            continue
        if entry_code != code:
            continue

        tokens[1] = rendered
        indent = line[: len(line) - len(line.lstrip())]
        rebuilt = indent + "    ".join(tokens)
        if marker:
            rebuilt += "  #" + comment.rstrip("\n")
        lines[index] = rebuilt + ("\n" if line.endswith("\n") else "")
        replaced = True
        break

    _require(replaced, f"could not find {block}({code}) in {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines), encoding="utf-8")


def _coefficient_from_payload(payload: dict[str, object], name: str, contribution: str) -> complex:
    coefficients = payload.get("coefficients")
    _require(isinstance(coefficients, dict), "invalid integration coefficient payload")
    raw_coefficient = coefficients.get(name)
    _require(isinstance(raw_coefficient, dict), f"coefficient {name} is absent from integration payload")
    raw_value = raw_coefficient.get(contribution)
    _require(isinstance(raw_value, dict), f"{name} {contribution} is absent from integration payload")
    return _decode_complex(raw_value)

def _complex_payload(value: object) -> dict[str, float]:
    z = complex(value)  # Scalar implements __complex__.
    return {"real": z.real, "imag": z.imag}


def _integration_child(args: argparse.Namespace) -> int:
    # Keep imports in the child so the default static suite does not require the
    # compiled native extension.
    sys.path.insert(0, str(PACKAGE_ROOT))

    from pyhyperiso.Common import (  # type: ignore[import-not-found]
        ContributionType,
        Model,
        QCDOrder,
        WCoeff,
        WGroup,
    )
    from pyhyperiso.Core import (  # type: ignore[import-not-found]
        ExternalFlag,
        HyperisoConfig,
        HyperisoMaster,
    )
    from pyhyperiso.Wilson import (  # type: ignore[import-not-found]
        WilsonBuildConfig,
        WilsonInterface,
    )

    hyp_as_sm = args._child_hyp_as_sm == "true"
    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: hyp_as_sm,
        },
        model=Model.MARTY,
        mty_model_name=args.model_name,
        mty_model_path=args.model_file,
        mty_bsm_mapping_path=str(args.mapping_file),
    )

    master = HyperisoMaster()
    master.pre_init_set_marty_path(str(args.marty_install))
    master.init(lha_file=str(args.lha_file), config=config)

    interface = WilsonInterface()
    interface.build(
        WilsonBuildConfig(
            groups={WGroup.B},
            matching_scale=args.matching_scale,
            hadronic_scale=args.hadronic_scale,
            order=QCDOrder.LO,
        )
    )

    coefficients: dict[str, dict[str, dict[str, float]]] = {}
    for index in range(1, 11):
        coefficient = getattr(WCoeff, f"C{index}")
        values: dict[str, dict[str, float]] = {}
        for contribution in (
            ContributionType.SM,
            ContributionType.BSM,
            ContributionType.TOTAL,
        ):
            values[contribution.name] = _complex_payload(
                interface.get_M(
                    WGroup.B,
                    coefficient,
                    QCDOrder.LO,
                    contribution,
                )
            )
        coefficients[coefficient.name] = values

    payload = {
        "hyp_as_sm_marty": hyp_as_sm,
        "model_name": args.model_name,
        "coefficients": coefficients,
    }
    print(CHILD_RESULT_PREFIX + json.dumps(payload, sort_keys=True))
    return 0


def _decode_complex(payload: dict[str, float]) -> complex:
    return complex(float(payload["real"]), float(payload["imag"]))


def _parse_child_payload(stdout: str) -> dict[str, object]:
    for line in reversed(stdout.splitlines()):
        if line.startswith(CHILD_RESULT_PREFIX):
            return json.loads(line.removeprefix(CHILD_RESULT_PREFIX))
    raise ValidationError("integration child returned no structured result")


def _cache_snapshot(cache_dir: Path, model_name: str) -> dict[str, tuple[int, int]]:
    if not cache_dir.is_dir():
        return {}

    candidates: set[Path] = set()
    for output_model in (model_name, "SM"):
        candidates.update(cache_dir.glob(f"generated_{output_model}_*.cpp"))
        candidates.update(
            path
            for path in cache_dir.glob(f"generated_{output_model}_*")
            if path.is_file() and path.suffix != ".cpp"
        )
        candidates.update(cache_dir.glob(f"libs/*_{output_model}/bin/*.x"))
        candidates.update(cache_dir.glob(f"libs/*_{output_model}/script/example_*.cpp"))

    snapshot: dict[str, tuple[int, int]] = {}
    for path in sorted(candidates):
        try:
            stat = path.stat()
        except FileNotFoundError:
            continue
        snapshot[str(path.relative_to(cache_dir))] = (stat.st_mtime_ns, stat.st_size)
    return snapshot


def _clean_model_cache(cache_dir: Path, model_name: str) -> list[str]:
    removed: list[str] = []
    if not cache_dir.is_dir():
        return removed
    for output_model in (model_name, "SM"):
        patterns = (
            f"generated_{output_model}_*",
            f"{output_model}_wilson.csv",
            f"libs/*_{output_model}",
        )
        for pattern in patterns:
            for path in cache_dir.glob(pattern):
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink(missing_ok=True)
                removed.append(str(path))
    return removed


def _validate_integration_payload(
    payload: dict[str, object],
    *,
    abs_tol: float,
    rel_tol: float,
    max_c2_total: float,
) -> str:
    coefficients = payload.get("coefficients")
    _require(isinstance(coefficients, dict), "invalid integration coefficient payload")

    worst_residual = 0.0
    for name, raw_values in coefficients.items():
        _require(isinstance(raw_values, dict), f"invalid values for {name}")
        values = raw_values
        sm = _decode_complex(values["SM"])
        bsm = _decode_complex(values["BSM"])
        total = _decode_complex(values["TOTAL"])
        for label, value in (("SM", sm), ("BSM", bsm), ("TOTAL", total)):
            _require(
                math.isfinite(value.real) and math.isfinite(value.imag),
                f"{name} {label} is not finite: {value}",
            )
        residual = abs(total - (sm + bsm))
        scale = max(abs(total), abs(sm) + abs(bsm), 1.0)
        worst_residual = max(worst_residual, residual)
        _require(
            residual <= abs_tol + rel_tol * scale,
            f"{name}: TOTAL != SM + BSM; residual={residual:.3e}, "
            f"SM={sm}, BSM={bsm}, TOTAL={total}",
        )

    c2_total = _decode_complex(coefficients["C2"]["TOTAL"])
    _require(
        abs(c2_total.real) <= max_c2_total,
        f"C2 TOTAL appears doubled: Re(C2)={c2_total.real:.8g} > {max_c2_total}",
    )
    mode = payload.get("hyp_as_sm_marty")
    return (
        f"HYP_AS_SM_MARTY={mode}: C1..C10 finite, composition OK; "
        f"max residual={worst_residual:.3e}, C2_TOTAL={c2_total.real:.8g}"
    )


def _integration_command(args: argparse.Namespace, mode: str, lha_file: Path) -> list[str]:
    return [
        sys.executable,
        str(Path(__file__).resolve()),
        "--_integration-child",
        "--_child-hyp-as-sm",
        mode,
        "--root",
        str(ROOT),
        "--marty-install",
        str(args.marty_install),
        "--model-file",
        str(args.model_file),
        "--mapping-file",
        str(args.mapping_file),
        "--lha-file",
        str(lha_file),
        "--model-name",
        args.model_name,
        "--matching-scale",
        str(args.matching_scale),
        "--hadronic-scale",
        str(args.hadronic_scale),
    ]


def _run_integration_child(args: argparse.Namespace, mode: str, lha_file: Path) -> dict[str, object]:
    completed = _run_command(
        _integration_command(args, mode, lha_file),
        timeout=args.integration_timeout,
    )
    # Surface only the compact analytical-generation diagnostics.  In
    # particular these lines state whether C9/C10 selected TreeLevel or OneLoop
    # without flooding the validation output with the complete HyperIso banner.
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("[MARTY "):
            print(f"[HYP_AS_SM_MARTY={mode}] {stripped}")
    return _parse_child_payload(completed.stdout)



def _parse_coefficient_names(raw: str) -> list[str]:
    names: list[str] = []
    for item in raw.split(","):
        name = item.strip().upper()
        if not name:
            continue
        _require(
            re.fullmatch(r"C(?:10|[1-9])", name) is not None,
            f"unsupported coefficient in --report-coefficients: {item}",
        )
        if name not in names:
            names.append(name)
    return names


def _render_selected_coefficients(
    payload: dict[str, object], names: Sequence[str]
) -> str:
    rendered: list[str] = []
    for name in names:
        sm = _coefficient_from_payload(payload, name, "SM")
        bsm = _coefficient_from_payload(payload, name, "BSM")
        total = _coefficient_from_payload(payload, name, "TOTAL")
        rendered.append(
            f"{name}[SM={sm.real:.8g}{sm.imag:+.3g}j,"
            f" BSM={bsm.real:.8g}{bsm.imag:+.3g}j,"
            f" TOTAL={total.real:.8g}{total.imag:+.3g}j]"
        )
    return "; ".join(rendered)


def _validate_bsm_mode_invariance(
    first: dict[str, object],
    second: dict[str, object],
    *,
    abs_tol: float,
    rel_tol: float,
) -> str:
    worst = 0.0
    for index in range(1, 11):
        name = f"C{index}"
        left = _coefficient_from_payload(first, name, "BSM")
        right = _coefficient_from_payload(second, name, "BSM")
        residual = abs(left - right)
        scale = max(abs(left), abs(right), 1.0)
        worst = max(worst, residual)
        _require(
            residual <= abs_tol + rel_tol * scale,
            f"{name} BSM depends on HYP_AS_SM_MARTY: true={left}, false={right}, "
            f"residual={residual:.3e}",
        )
    return f"BSM is independent of the SM provider for C1..C10 (max residual={worst:.3e})"

def _validate_decoupling_scan(
    points: list[tuple[float, complex]],
    *,
    coefficient: str,
    max_ratio: float,
    abs_tol: float,
) -> str:
    _require(len(points) >= 2, "decoupling scan returned fewer than two points")
    magnitudes = [abs(value) for _, value in points]
    first = magnitudes[0]
    last = magnitudes[-1]
    allowed = max(abs_tol, max_ratio * first)
    _require(
        last <= allowed,
        f"{coefficient} BSM does not decouple: |C|={last:.6e} at mass={points[-1][0]:g}, "
        f"expected <= {allowed:.6e} from |C|={first:.6e} at mass={points[0][0]:g}",
    )

    # Allow small numerical/non-monotonic fluctuations, but reject a scan whose
    # heavy end systematically grows.  This is deliberately weaker than strict
    # monotonicity because loop cancellations can produce local extrema.
    tail_min = min(magnitudes[:-1])
    _require(
        last <= max(abs_tol, 1.25 * tail_min),
        f"{coefficient} BSM heavy-mass endpoint grows above the scan minimum: "
        f"last={last:.6e}, previous minimum={tail_min:.6e}",
    )
    rendered = ", ".join(f"{mass:g}:{abs(value):.3e}" for mass, value in points)
    return f"{coefficient} BSM decoupling [{rendered}]"


def run_integration(args: argparse.Namespace) -> str:
    for path, label in (
        (args.marty_install, "MARTY install"),
        (args.model_file, "MARTY model"),
        (args.mapping_file, "MARTY mapping"),
        (args.lha_file, "LHA input"),
    ):
        _require(path.exists(), f"{label} does not exist: {path}")

    cache_dir = args.cache_dir.expanduser().resolve()
    if args.clean_model_cache:
        removed = _clean_model_cache(cache_dir, args.model_name)
        print(f"[INFO] removed {len(removed)} model-specific cache entries from {cache_dir}")

    modes = {
        "true": ["true"],
        "false": ["false"],
        "both": ["true", "false"],
    }[args.hyp_as_sm_marty]

    scan_masses = _parse_float_list(args.zprime_mass_scan) if args.zprime_mass_scan else []
    report_coefficients = _parse_coefficient_names(args.report_coefficients)
    coefficient_name = args.decoupling_coefficient.upper()
    _require(
        re.fullmatch(r"C(?:10|[1-9])", coefficient_name) is not None,
        f"unsupported decoupling coefficient: {args.decoupling_coefficient}",
    )

    summaries: list[str] = []
    payloads_by_mode: dict[str, dict[str, object]] = {}
    with tempfile.TemporaryDirectory(prefix="hyperiso-marty-v103-") as tmp_raw:
        tmp = Path(tmp_raw)
        base_lha = args.lha_file.resolve()
        if args.thdm_type is not None:
            typed_lha = tmp / f"thdm_type_{args.thdm_type}.lha"
            _rewrite_lha_scalar(base_lha, typed_lha, "MINPAR", 24, float(args.thdm_type))
            base_lha = typed_lha

        for mode in modes:
            first_payload = _run_integration_child(args, mode, base_lha)
            payloads_by_mode[mode] = first_payload
            summary = _validate_integration_payload(
                first_payload,
                abs_tol=args.abs_tol,
                rel_tol=args.rel_tol,
                max_c2_total=args.max_c2_total,
            )
            if report_coefficients:
                summary += "; " + _render_selected_coefficients(
                    first_payload, report_coefficients
                )

            if not args.no_cache_check:
                first_snapshot = _cache_snapshot(cache_dir, args.model_name)
                _require(first_snapshot, f"no analytical cache artifacts found in {cache_dir}")
                # Filesystems with coarse timestamps benefit from a small boundary.
                time.sleep(1.05)
                second_payload = _run_integration_child(args, mode, base_lha)
                _validate_integration_payload(
                    second_payload,
                    abs_tol=args.abs_tol,
                    rel_tol=args.rel_tol,
                    max_c2_total=args.max_c2_total,
                )
                second_snapshot = _cache_snapshot(cache_dir, args.model_name)
                changed = {
                    path: (first_snapshot.get(path), second_snapshot.get(path))
                    for path in sorted(set(first_snapshot) | set(second_snapshot))
                    if first_snapshot.get(path) != second_snapshot.get(path)
                }
                _require(
                    not changed,
                    "analytical MARTY cache changed on identical second run: "
                    + json.dumps(changed, indent=2),
                )
                summary += f"; {len(second_snapshot)} analytical cache artifacts reused"

            if scan_masses:
                points: list[tuple[float, complex]] = []
                for index, mass in enumerate(scan_masses):
                    scan_lha = tmp / f"mass32_{mode}_{index}_{mass:g}.lha"
                    _rewrite_lha_scalar(base_lha, scan_lha, "MASS", 32, mass)
                    payload = _run_integration_child(args, mode, scan_lha)
                    _validate_integration_payload(
                        payload,
                        abs_tol=args.abs_tol,
                        rel_tol=args.rel_tol,
                        max_c2_total=args.max_c2_total,
                    )
                    points.append(
                        (
                            mass,
                            _coefficient_from_payload(payload, coefficient_name, "BSM"),
                        )
                    )
                summary += "; " + _validate_decoupling_scan(
                    points,
                    coefficient=coefficient_name,
                    max_ratio=args.decoupling_max_ratio,
                    abs_tol=args.decoupling_abs_tol,
                )

            summaries.append(summary)

    if "true" in payloads_by_mode and "false" in payloads_by_mode:
        summaries.append(
            _validate_bsm_mode_invariance(
                payloads_by_mode["true"],
                payloads_by_mode["false"],
                abs_tol=args.abs_tol,
                rel_tol=args.rel_tol,
            )
        )

    return " | ".join(summaries)


def _default_paths(root: Path) -> tuple[Path, Path, Path]:
    assets = root / "Hyperiso" / "Hyperiso" / "pyhyperiso" / "assets"
    return (
        assets / "input_files" / "marty_model" / "ZPrime.h",
        assets / "input_files" / "marty_mapping" / "zprime.json",
        assets / "lha" / "zprime_input.flha",
    )


def parse_args() -> argparse.Namespace:
    default_model, default_mapping, default_lha = _default_paths(ROOT)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help=argparse.SUPPRESS)
    parser.add_argument(
        "--static-only",
        action="store_true",
        help="run only dependency-light repository checks",
    )
    parser.add_argument(
        "--require-import",
        action="store_true",
        help="fail instead of skipping when the compiled pyhyperiso extension cannot be imported",
    )
    parser.add_argument(
        "--ctest-build-dir",
        type=Path,
        help="existing CMake build directory; runs ctest -R 'Marty|Wilson'",
    )
    parser.add_argument(
        "--integration",
        action="store_true",
        help="run the real Z-prime MARTY coefficient and cache checks",
    )
    parser.add_argument(
        "--marty-install",
        type=Path,
        default=Path(os.environ["MARTY_INSTALL"]) if "MARTY_INSTALL" in os.environ else None,
        help="MARTY installation prefix (or set MARTY_INSTALL)",
    )
    parser.add_argument("--model-file", type=Path, default=default_model)
    parser.add_argument("--mapping-file", type=Path, default=default_mapping)
    parser.add_argument("--lha-file", type=Path, default=default_lha)
    parser.add_argument("--model-name", default="ZPrime_Model")
    parser.add_argument(
        "--thdm-type",
        type=int,
        choices=(1, 2, 3, 4),
        help="override MINPAR(24) for a THDM MARTY integration run",
    )
    parser.add_argument(
        "--zprime-mass-scan",
        help="comma-separated MASS(32) values for a BSM decoupling scan, e.g. 100,300,1000,3000",
    )
    parser.add_argument(
        "--decoupling-coefficient",
        default="C9",
        help="C1..C10 coefficient checked by --zprime-mass-scan (default: C9)",
    )
    parser.add_argument(
        "--report-coefficients",
        default="C2,C7,C8,C9,C10",
        help="comma-separated coefficients printed for each integration mode",
    )
    parser.add_argument(
        "--decoupling-max-ratio",
        type=float,
        default=0.2,
        help="maximum allowed |C_BSM(last)| / |C_BSM(first)| in the mass scan",
    )
    parser.add_argument(
        "--decoupling-abs-tol",
        type=float,
        default=1e-8,
        help="absolute decoupling tolerance for the heavy-mass endpoint",
    )
    parser.add_argument(
        "--hyp-as-sm-marty",
        choices=("true", "false", "both"),
        default="both",
        help="integration scenarios to test",
    )
    parser.add_argument("--matching-scale", type=float, default=81.0)
    parser.add_argument("--hadronic-scale", type=float, default=2.0)
    parser.add_argument("--abs-tol", type=float, default=1e-10)
    parser.add_argument("--rel-tol", type=float, default=1e-8)
    parser.add_argument(
        "--max-c2-total",
        type=float,
        default=1.5,
        help="guard against the historical C2 ~= 2 double counting",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path.home() / ".cache" / "pyhyperiso" / "MartyTemp",
    )
    parser.add_argument(
        "--clean-model-cache",
        action="store_true",
        help="remove only the selected model and SM generated MARTY cache entries before integration",
    )
    parser.add_argument(
        "--no-cache-check",
        action="store_true",
        help="do not repeat the integration run to verify analytical cache reuse",
    )
    parser.add_argument("--integration-timeout", type=float, default=1800.0)
    parser.add_argument("--json-report", type=Path, help="write machine-readable results")

    # Internal child-process arguments.
    parser.add_argument("--_integration-child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--_child-hyp-as-sm",
        choices=("true", "false"),
        default="true",
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def main() -> int:
    global ROOT, PACKAGE_ROOT, PYTHON_PACKAGE, PACKAGED_ASSETS, CANONICAL_TEMPLATES, MIRROR_TEMPLATES

    args = parse_args()
    ROOT = args.root.resolve()
    PACKAGE_ROOT = ROOT / "Hyperiso" / "Hyperiso"
    PYTHON_PACKAGE = PACKAGE_ROOT / "pyhyperiso"
    PACKAGED_ASSETS = PYTHON_PACKAGE / "assets"
    CANONICAL_TEMPLATES = PACKAGED_ASSETS / "template" / "MARTY"
    MIRROR_TEMPLATES = ROOT / "Assets" / "template" / "MARTY"

    if args._integration_child:
        if args.marty_install is None:
            print("--marty-install is required", file=sys.stderr)
            return 2
        return _integration_child(args)

    checks: list[tuple[str, Callable[[], str]]] = [
        ("repository layout", check_repository_layout),
        ("release version", check_version_consistency),
        ("MARTY template mirrors", check_template_sync),
        ("MARTY template ABI", check_template_abi),
        ("Z-prime mapping", check_zprime_mapping),
        ("packaged model portability", check_packaged_model_portability),
        ("packaged SM path", check_sm_path_resolution),
        ("SM/BSM/TOTAL source composition", check_contribution_composition_sources),
        ("tree-first semileptonic policy", check_tree_first_policy),
        ("MARTY cache signatures", check_cache_source_signatures),
        ("Python compilation", check_python_compilation),
    ]
    if not args.static_only:
        checks.append(("pyhyperiso import", check_import_smoke))
    if args.ctest_build_dir:
        checks.append(("C++ MARTY/Wilson tests", lambda: run_ctest(args.ctest_build_dir)))
    if args.integration:
        _require(args.marty_install is not None, "--integration requires --marty-install or MARTY_INSTALL")
        checks.append(("Z-prime MARTY integration", lambda: run_integration(args)))

    results: list[CheckResult] = []
    for name, callback in checks:
        started = time.monotonic()
        try:
            detail = callback()
        except FileNotFoundError as exc:
            status = "FAIL" if args.require_import and name == "pyhyperiso import" else "SKIP"
            result = CheckResult(name, status, str(exc), time.monotonic() - started)
        except (ValidationError, OSError, subprocess.SubprocessError, ValueError) as exc:
            result = CheckResult(name, "FAIL", str(exc), time.monotonic() - started)
        except Exception as exc:  # Defensive reporting for native/runtime failures.
            result = CheckResult(
                name,
                "FAIL",
                f"{type(exc).__name__}: {exc}",
                time.monotonic() - started,
            )
        else:
            result = CheckResult(name, "PASS", detail, time.monotonic() - started)
        results.append(result)
        print(f"[{result.status}] {result.name}: {result.detail} ({result.duration_s:.2f}s)")

    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        args.json_report.write_text(
            json.dumps([asdict(result) for result in results], indent=2) + "\n",
            encoding="utf-8",
        )

    failures = [result for result in results if result.status == "FAIL"]
    skipped = [result for result in results if result.status == "SKIP"]
    print(
        f"\nSummary: {len(results) - len(failures) - len(skipped)} passed, "
        f"{len(skipped)} skipped, {len(failures)} failed."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
