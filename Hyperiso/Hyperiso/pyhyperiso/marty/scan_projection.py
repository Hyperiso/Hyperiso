"""Command-line diagnostic for MARTY TreeLevel projection recipes.

The command never runs inside a production HyperIso process.  It repeatedly
launches a user-supplied isolated runner, ranks the candidate recipes, and can
freeze a validated candidate to JSON.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .projection import (
    SEMILEPTONIC_VECTOR_COEFFICIENTS,
    SUPPORTED_CURRENTS,
    SUPPORTED_LAYOUTS,
    MartyProjectionProfile,
    SubprocessProjectionScanner,
    all_fermion_orders,
)


def _order(text: str) -> tuple[int, int, int, int]:
    try:
        values = tuple(int(value.strip()) for value in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if len(values) != 4 or sorted(values) != [0, 1, 2, 3]:
        raise argparse.ArgumentTypeError("order must be a permutation such as 0,2,1,3")
    return values  # type: ignore[return-value]


def _current_pair(text: str) -> tuple[str, str]:
    fields = tuple(part.strip() for part in text.split(","))
    if len(fields) != 2 or any(field not in SUPPORTED_CURRENTS for field in fields):
        raise argparse.ArgumentTypeError(
            "current pair must contain two of " + ", ".join(sorted(SUPPORTED_CURRENTS))
        )
    return fields  # type: ignore[return-value]


def _complex_target(real: float | None, imag: float | None) -> complex | None:
    if real is None and imag is None:
        return None
    return complex(0.0 if real is None else real, 0.0 if imag is None else imag)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Scan isolated MARTY TreeLevel projection recipes."
    )
    parser.add_argument("--runner", required=True, help="isolated runner command")
    parser.add_argument(
        "--coefficient",
        required=True,
        choices=sorted(SEMILEPTONIC_VECTOR_COEFFICIENTS),
    )
    parser.add_argument("--cwd", type=Path, help="working directory for the runner")
    parser.add_argument("--timeout", type=float, default=None, help="timeout per candidate in seconds")

    parser.add_argument("--fermion-order", action="append", type=_order, default=[])
    parser.add_argument("--operator-order", action="append", type=_order, default=[])
    parser.add_argument("--all-fermion-orders", action="store_true")
    parser.add_argument("--all-operator-orders", action="store_true")

    parser.add_argument("--current-pair", action="append", type=_current_pair, default=[])
    parser.add_argument("--all-current-pairs", action="store_true")
    parser.add_argument(
        "--layout",
        action="append",
        choices=sorted(SUPPORTED_LAYOUTS),
        default=[],
    )
    parser.add_argument("--all-layouts", action="store_true")

    parser.add_argument("--expected-real", type=float)
    parser.add_argument("--expected-imag", type=float)
    parser.add_argument("--max-candidates", type=int, default=10000)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--export-best", type=Path)
    parser.add_argument(
        "--allow-magnitude-ranking",
        action="store_true",
        help=(
            "allow --export-best without an expected value. This is diagnostic only: "
            "largest/non-zero is not a physics validation."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.all_fermion_orders:
        fermion_orders = all_fermion_orders()
    elif args.fermion_order:
        fermion_orders = tuple(args.fermion_order)
    else:
        raise SystemExit("Specify --fermion-order at least once, or --all-fermion-orders")

    if args.all_operator_orders:
        operator_orders = all_fermion_orders()
    elif args.operator_order:
        operator_orders = tuple(args.operator_order)
    else:
        raise SystemExit("Specify --operator-order at least once, or --all-operator-orders")

    if args.all_current_pairs:
        current_pairs = tuple(
            (left, right)
            for left in sorted(SUPPORTED_CURRENTS)
            for right in sorted(SUPPORTED_CURRENTS)
        )
    elif args.current_pair:
        current_pairs = tuple(args.current_pair)
    else:
        raise SystemExit("Specify --current-pair at least once, or --all-current-pairs")

    if args.all_layouts:
        layouts = tuple(sorted(SUPPORTED_LAYOUTS))
    elif args.layout:
        layouts = tuple(args.layout)
    else:
        layouts = ("quark_first",)

    expected = _complex_target(args.expected_real, args.expected_imag)
    scanner = SubprocessProjectionScanner(
        args.runner,
        cwd=args.cwd,
        timeout=args.timeout,
    )
    results = scanner.scan_single_term(
        coefficient=args.coefficient,
        fermion_orders=fermion_orders,
        operator_orders=operator_orders,
        current_pairs=current_pairs,
        layouts=layouts,
        expected=expected,
        max_candidates=args.max_candidates,
    )

    successful = [result for result in results if result.value is not None]
    failed = len(results) - len(successful)
    print(
        f"MARTY projection scan: {len(results)} candidates, "
        f"{len(successful)} successful, {failed} failed"
    )
    if expected is not None:
        print(f"reference: {expected.real:+.12e}{expected.imag:+.12e}j")
    else:
        print("ranking: |C| only (diagnostic; non-zero does not establish correctness)")

    for rank, result in enumerate(successful[: max(args.top, 0)], start=1):
        recipe = result.profile.tree[args.coefficient]
        term = recipe.terms[0]
        score = "n/a" if result.score is None else f"{result.score:.6e}"
        value = result.value
        assert value is not None
        print(
            f"{rank:3d}  F={list(term.fermion_order)}  O={list(term.operator_order)}  "
            f"J=({term.left_current},{term.right_current})  layout={term.layout:12s}  "
            f"C={value.real:+.12e}{value.imag:+.12e}j  score={score}"
        )

    if args.export_best is not None:
        if not successful:
            raise SystemExit("Cannot export: no successful projection candidate")
        if expected is None and not args.allow_magnitude_ranking:
            raise SystemExit(
                "Refusing to export a largest/non-zero candidate as validated. "
                "Provide --expected-real/--expected-imag or explicitly pass "
                "--allow-magnitude-ranking for diagnostic use."
            )
        best: MartyProjectionProfile = successful[0].profile
        best.save(args.export_best)
        print(f"exported: {args.export_best}")

    return 0 if successful else 2


if __name__ == "__main__":
    sys.exit(main())
