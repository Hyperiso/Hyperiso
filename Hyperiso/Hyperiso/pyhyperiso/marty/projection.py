"""Configurable MARTY TreeLevel projection recipes and diagnostic scans.

The public HyperIso configuration keeps its existing per-coefficient fermion/order
maps.  Projection recipes are encoded into reserved keys of those maps so the
feature remains ABI-light and backwards compatible with the current pybind
configuration object.

A recipe is explicit and frozen: HyperIso never scans permutations during normal
predictions or fits.  The scanner helpers below run candidate profiles in
isolated subprocesses and export the selected profile as JSON.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations, product
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
from typing import Any, Iterable, Mapping, MutableMapping, Sequence

TREE_RECIPE_PREFIX = "__HYPERISO_MARTY_TREE_RECIPE__"
SUPPORTED_CURRENTS = frozenset({"VL", "VR", "V", "A"})
SUPPORTED_LAYOUTS = frozenset({"quark_first", "lepton_first"})
SEMILEPTONIC_VECTOR_COEFFICIENTS = frozenset({"C9", "C10", "CP9", "CP10"})
PROFILE_ENV = "HYPERISO_MARTY_PROJECTION_PROFILE"
SCAN_RESULT_ENV = "HYPERISO_MARTY_SCAN_RESULT"


def _normalise_order(order: Sequence[int], *, field_name: str) -> tuple[int, int, int, int]:
    values = tuple(int(v) for v in order)
    if len(values) != 4 or sorted(values) != [0, 1, 2, 3]:
        raise ValueError(f"{field_name} must be a permutation of [0,1,2,3]; received {list(values)}")
    return values  # type: ignore[return-value]


@dataclass(frozen=True)
class MartyProjectionTerm:
    """One weighted TreeLevel MARTY projector.

    ``fermion_order`` controls amplitude/Fierz matching (F); ``operator_order``
    controls ``dimension6Operator`` pairing (O).  ``layout='lepton_first'`` swaps
    the two Dirac-current arguments at projection time.
    """

    name: str
    weight: float
    left_current: str
    right_current: str
    fermion_order: Sequence[int]
    operator_order: Sequence[int]
    layout: str = "quark_first"

    def __post_init__(self) -> None:
        if not str(self.name).strip():
            raise ValueError("Projection term name must not be empty")
        if not math.isfinite(float(self.weight)):
            raise ValueError("Projection term weight must be finite")
        if self.left_current not in SUPPORTED_CURRENTS or self.right_current not in SUPPORTED_CURRENTS:
            raise ValueError(
                f"Recipe ABI v1 supports currents {sorted(SUPPORTED_CURRENTS)}; "
                f"received ({self.left_current}, {self.right_current})"
            )
        if self.layout not in SUPPORTED_LAYOUTS:
            raise ValueError(f"layout must be one of {sorted(SUPPORTED_LAYOUTS)}")
        object.__setattr__(self, "fermion_order", _normalise_order(self.fermion_order, field_name="fermion_order"))
        object.__setattr__(self, "operator_order", _normalise_order(self.operator_order, field_name="operator_order"))

    def reserved_key(self, coefficient: str, index: int) -> str:
        # Use a stable numeric prefix so lexicographic C++ sorting preserves term order.
        return "|".join(
            (
                TREE_RECIPE_PREFIX,
                coefficient,
                f"{index:04d}_{self.name}",
                format(float(self.weight), ".17g"),
                self.left_current,
                self.right_current,
                self.layout,
            )
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "weight": float(self.weight),
            "left_current": self.left_current,
            "right_current": self.right_current,
            "fermion_order": list(self.fermion_order),
            "operator_order": list(self.operator_order),
            "layout": self.layout,
        }

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> "MartyProjectionTerm":
        return cls(
            name=str(data["name"]),
            weight=float(data["weight"]),
            left_current=str(data["left_current"]),
            right_current=str(data["right_current"]),
            fermion_order=data["fermion_order"],
            operator_order=data["operator_order"],
            layout=str(data.get("layout", "quark_first")),
        )


@dataclass(frozen=True)
class MartyProjectionRecipe:
    coefficient: str
    terms: Sequence[MartyProjectionTerm] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        coefficient = str(self.coefficient).strip()
        if coefficient not in SEMILEPTONIC_VECTOR_COEFFICIENTS:
            raise ValueError(
                "Recipe ABI v1 is enabled for C9/C10/CP9/CP10; "
                f"received {coefficient!r}"
            )
        if not self.terms:
            raise ValueError("A projection recipe must contain at least one term")
        object.__setattr__(self, "coefficient", coefficient)
        object.__setattr__(self, "terms", tuple(self.terms))

    def to_json(self) -> dict[str, Any]:
        return {"terms": [term.to_json() for term in self.terms]}

    @classmethod
    def from_json(cls, coefficient: str, data: Mapping[str, Any]) -> "MartyProjectionRecipe":
        return cls(
            coefficient=coefficient,
            terms=tuple(MartyProjectionTerm.from_json(term) for term in data["terms"]),
        )


@dataclass(frozen=True)
class MartyProjectionProfile:
    tree: Mapping[str, MartyProjectionRecipe]
    version: int = 1

    def __post_init__(self) -> None:
        if self.version != 1:
            raise ValueError(f"Unsupported MARTY projection profile version {self.version}")
        for coefficient, recipe in self.tree.items():
            if coefficient != recipe.coefficient:
                raise ValueError(f"Profile key {coefficient!r} does not match recipe {recipe.coefficient!r}")

    def to_json(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "tree": {coefficient: recipe.to_json() for coefficient, recipe in self.tree.items()},
        }

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_json(), indent=2, sort_keys=True) + "\n")
        return path

    @classmethod
    def load(cls, path: str | Path) -> "MartyProjectionProfile":
        data = json.loads(Path(path).read_text())
        version = int(data.get("version", 1))
        tree = {
            coefficient: MartyProjectionRecipe.from_json(coefficient, recipe)
            for coefficient, recipe in data.get("tree", {}).items()
        }
        return cls(tree=tree, version=version)


def _remove_recipe_keys(mapping: MutableMapping[str, Sequence[int]], coefficient: str) -> None:
    prefix = f"{TREE_RECIPE_PREFIX}|{coefficient}|"
    for key in [key for key in mapping if str(key).startswith(prefix)]:
        del mapping[key]


def apply_tree_projection_profile(config: Any, profile: MartyProjectionProfile) -> Any:
    """Merge a frozen profile into an existing ``HyperisoConfig``.

    The recipe is authoritative term-by-term.  Coarse per-coefficient TreeLevel
    F/O entries for the same coefficient are removed; all unrelated coefficients
    and every OneLoop override are preserved.
    """

    fermion_orders: dict[str, Sequence[int]] = dict(getattr(config, "mty_tree_fermion_orders", {}))
    operator_orders: dict[str, Sequence[int]] = dict(getattr(config, "mty_tree_operator_orders", {}))

    for coefficient, recipe in profile.tree.items():
        _remove_recipe_keys(fermion_orders, coefficient)
        _remove_recipe_keys(operator_orders, coefficient)
        fermion_orders.pop(coefficient, None)
        operator_orders.pop(coefficient, None)
        for index, term in enumerate(recipe.terms):
            key = term.reserved_key(coefficient, index)
            fermion_orders[key] = list(term.fermion_order)
            operator_orders[key] = list(term.operator_order)

    config.mty_tree_fermion_orders = fermion_orders
    config.mty_tree_operator_orders = operator_orders
    return config


def apply_projection_profile_from_env(config: Any, *, env: Mapping[str, str] | None = None) -> Any:
    env = os.environ if env is None else env
    path = env.get(PROFILE_ENV)
    if not path:
        return config
    return apply_tree_projection_profile(config, MartyProjectionProfile.load(path))


def direct_semileptonic_profile(
    fermion_order: Sequence[int],
    operator_order: Sequence[int],
    *,
    layout: str = "quark_first",
) -> MartyProjectionProfile:
    """Direct V/A basis profile useful for neutral-vector exchange models."""

    f = _normalise_order(fermion_order, field_name="fermion_order")
    o = _normalise_order(operator_order, field_name="operator_order")
    currents = {
        "C9": ("VL", "V"),
        "C10": ("VL", "A"),
        "CP9": ("VR", "V"),
        "CP10": ("VR", "A"),
    }
    return MartyProjectionProfile(
        tree={
            coefficient: MartyProjectionRecipe(
                coefficient,
                (MartyProjectionTerm("direct", 1.0, left, right, f, o, layout),),
            )
            for coefficient, (left, right) in currents.items()
        }
    )


def lq_semileptonic_chiral_profile(
    fermion_order: Sequence[int] = (1, 0, 2, 3),
) -> MartyProjectionProfile:
    """Validated LQ TreeLevel LL/LR/RL/RR decomposition."""

    f = _normalise_order(fermion_order, field_name="fermion_order")
    o_ll = (0, 3, 1, 2)
    o_other = (1, 2, 0, 3)
    return MartyProjectionProfile(
        tree={
            "C9": MartyProjectionRecipe("C9", (
                MartyProjectionTerm("LL", +0.5, "VL", "VL", f, o_ll),
                MartyProjectionTerm("LR", +0.5, "VL", "VR", f, o_other),
            )),
            "C10": MartyProjectionRecipe("C10", (
                MartyProjectionTerm("LL", -0.5, "VL", "VL", f, o_ll),
                MartyProjectionTerm("LR", +0.5, "VL", "VR", f, o_other),
            )),
            "CP9": MartyProjectionRecipe("CP9", (
                MartyProjectionTerm("RL", +0.5, "VR", "VL", f, o_other),
                MartyProjectionTerm("RR", +0.5, "VR", "VR", f, o_other),
            )),
            "CP10": MartyProjectionRecipe("CP10", (
                MartyProjectionTerm("RL", -0.5, "VR", "VL", f, o_other),
                MartyProjectionTerm("RR", +0.5, "VR", "VR", f, o_other),
            )),
        }
    )


def write_projection_scan_result(values: Mapping[str, complex], path: str | Path | None = None) -> Path:
    """Write a scanner-compatible JSON result from an isolated runner process."""

    target = Path(path or os.environ.get(SCAN_RESULT_ENV, "projection_scan_result.json"))
    payload = {name: [float(complex(value).real), float(complex(value).imag)] for name, value in values.items()}
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return target


@dataclass(frozen=True)
class ProjectionScanCandidate:
    profile: MartyProjectionProfile
    value: complex | None
    score: float | None
    stdout: str
    stderr: str
    error: str | None = None


class SubprocessProjectionScanner:
    """Explore projection candidates without mutating the production process.

    ``runner`` is a command that initializes HyperIso/MARTY for one benchmark.
    The runner should call ``apply_projection_profile_from_env(config)`` before
    ``HyperisoMaster.init`` and ``write_projection_scan_result`` after computing
    the requested Wilson coefficient.
    """

    def __init__(self, runner: Sequence[str] | str, *, cwd: str | Path | None = None, timeout: float | None = None):
        self.runner = shlex.split(runner) if isinstance(runner, str) else [str(x) for x in runner]
        self.cwd = None if cwd is None else Path(cwd)
        self.timeout = timeout

    def evaluate(
        self,
        profile: MartyProjectionProfile,
        coefficient: str,
        *,
        expected: complex | None = None,
        strict: bool = False,
    ) -> ProjectionScanCandidate:
        with tempfile.TemporaryDirectory(prefix="pyhyperiso-marty-projection-") as tmp:
            tmpdir = Path(tmp)
            profile_path = profile.save(tmpdir / "profile.json")
            result_path = tmpdir / "result.json"
            env = os.environ.copy()
            env[PROFILE_ENV] = str(profile_path)
            env[SCAN_RESULT_ENV] = str(result_path)
            try:
                proc = subprocess.run(
                    self.runner,
                    cwd=self.cwd,
                    env=env,
                    text=True,
                    capture_output=True,
                    timeout=self.timeout,
                    check=False,
                )
            except (OSError, subprocess.TimeoutExpired) as exc:
                if strict:
                    raise
                return ProjectionScanCandidate(profile, None, None, "", str(exc), str(exc))

            error = None
            value: complex | None = None
            if proc.returncode != 0:
                error = f"runner exit code {proc.returncode}"
            elif not result_path.is_file():
                error = f"runner did not write {SCAN_RESULT_ENV}={result_path}"
            else:
                try:
                    payload = json.loads(result_path.read_text())
                    raw = payload[coefficient]
                    value = complex(float(raw[0]), float(raw[1]))
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    error = f"invalid scanner result: {exc}"

            if error is not None:
                if strict:
                    raise RuntimeError(
                        f"Projection runner failed: {error}\n"
                        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
                    )
                return ProjectionScanCandidate(
                    profile, None, None, proc.stdout, proc.stderr, error
                )

            assert value is not None
            score = None
            if expected is not None:
                scale = max(abs(expected), 1.0e-30)
                score = abs(value - expected) / scale
            return ProjectionScanCandidate(
                profile, value, score, proc.stdout, proc.stderr, None
            )

    def scan_single_term(
        self,
        *,
        coefficient: str,
        fermion_orders: Iterable[Sequence[int]],
        operator_orders: Iterable[Sequence[int]],
        current_pairs: Iterable[tuple[str, str]],
        layouts: Iterable[str] = ("quark_first", "lepton_first"),
        expected: complex | None = None,
        max_candidates: int = 10000,
    ) -> list[ProjectionScanCandidate]:
        fs = [_normalise_order(order, field_name="fermion_order") for order in fermion_orders]
        os_ = [_normalise_order(order, field_name="operator_order") for order in operator_orders]
        pairs = list(current_pairs)
        layouts_ = list(layouts)
        total = len(fs) * len(os_) * len(pairs) * len(layouts_)
        if total > max_candidates:
            raise ValueError(
                f"Projection scan would run {total} candidates; raise max_candidates explicitly "
                "or restrict F/O/current/layout candidates."
            )

        results: list[ProjectionScanCandidate] = []
        for index, (f, o, pair, layout) in enumerate(product(fs, os_, pairs, layouts_)):
            term = MartyProjectionTerm(
                name=f"scan_{index:05d}",
                weight=1.0,
                left_current=pair[0],
                right_current=pair[1],
                fermion_order=f,
                operator_order=o,
                layout=layout,
            )
            profile = MartyProjectionProfile(
                tree={coefficient: MartyProjectionRecipe(coefficient, (term,))}
            )
            results.append(self.evaluate(profile, coefficient, expected=expected))

        if expected is not None:
            results.sort(key=lambda item: float("inf") if item.score is None else item.score)
        else:
            results.sort(
                key=lambda item: -1.0 if item.value is None else abs(item.value),
                reverse=True,
            )
        return results


def all_fermion_orders() -> tuple[tuple[int, int, int, int], ...]:
    return tuple(permutations(range(4)))


__all__ = [
    "MartyProjectionTerm",
    "MartyProjectionRecipe",
    "MartyProjectionProfile",
    "ProjectionScanCandidate",
    "SubprocessProjectionScanner",
    "apply_tree_projection_profile",
    "apply_projection_profile_from_env",
    "direct_semileptonic_profile",
    "lq_semileptonic_chiral_profile",
    "write_projection_scan_result",
    "all_fermion_orders",
    "PROFILE_ENV",
    "SCAN_RESULT_ENV",
    "TREE_RECIPE_PREFIX",
]
