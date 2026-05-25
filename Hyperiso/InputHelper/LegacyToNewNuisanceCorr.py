#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from itertools import product, combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

NewTarget = Tuple[str, str]
NewTargets = List[NewTarget]
CanonicalTarget = Tuple[str, str]
CanonicalPair = Tuple[CanonicalTarget, CanonicalTarget]


@dataclass
class LegacyCorrelation:
    line_no: int
    raw_name_1: str
    raw_name_2: str
    stat_correlation: float
    syst_correlation: float
    raw_line: str


@dataclass
class ConvertedCorrelation:
    line_no: int
    raw_name_1: str
    raw_name_2: str
    block_1: str
    id_1: str
    block_2: str
    id_2: str
    stat_correlation: float
    syst_correlation: float


@dataclass
class SkippedCorrelation:
    line_no: int
    raw_name_1: str
    raw_name_2: str
    reason: str
    details: Dict[str, Any]


@dataclass
class DuplicateCorrelation:
    line_no: int
    raw_name_1: str
    raw_name_2: str
    block_1: str
    id_1: str
    block_2: str
    id_2: str
    reason: str
    previous: Dict[str, Any]
    new: Dict[str, Any]


_FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _looks_like_float(token: str) -> bool:
    return bool(_FLOAT_RE.match(token))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert legacy nuisance correlations to the new JSON format")
    parser.add_argument("input", help="Legacy nuisance-correlation text file")
    parser.add_argument("output", help="Output JSON file")
    parser.add_argument(
        "--maps",
        required=True,
        help="Path to nuisance_map.py (or compatible mapping module)",
    )
    parser.add_argument(
        "--report-json",
        default=None,
        help="Optional JSON report path",
    )
    parser.add_argument(
        "--merge-with",
        default=None,
        help="Optional existing JSON file to merge into before writing output",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level for output files (default: 2)",
    )
    parser.add_argument(
        "--fail-on-conflict",
        action="store_true",
        help="Exit non-zero if the same output pair is produced with different values",
    )
    parser.add_argument(
        "--fail-on-skipped",
        action="store_true",
        help="Exit non-zero if at least one legacy line could not be converted",
    )
    return parser.parse_args()


def load_module(module_path: Union[str, Path]):
    module_path = Path(module_path)
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import mapping module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_legacy_line(line: str, line_no: int) -> Optional[LegacyCorrelation]:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None

    parts = stripped.split()
    if len(parts) < 3:
        raise ValueError(f"line {line_no}: expected at least 3 columns, got {len(parts)}")

    name1, name2 = parts[0], parts[1]
    if not _looks_like_float(parts[2]):
        raise ValueError(f"line {line_no}: cannot parse stat correlation: {parts[2]!r}")

    stat_corr = float(parts[2])
    syst_corr = 0.0
    if len(parts) >= 4:
        if not _looks_like_float(parts[3]):
            raise ValueError(f"line {line_no}: cannot parse syst correlation: {parts[3]!r}")
        syst_corr = float(parts[3])

    return LegacyCorrelation(
        line_no=line_no,
        raw_name_1=name1,
        raw_name_2=name2,
        stat_correlation=stat_corr,
        syst_correlation=syst_corr,
        raw_line=line.rstrip("\n"),
    )


def iter_legacy_correlations(path: Union[str, Path]) -> Iterable[LegacyCorrelation]:
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            entry = parse_legacy_line(line, line_no)
            if entry is not None:
                yield entry


def resolve_target(module: Any, legacy_name: str) -> Optional[Union[NewTarget, NewTargets]]:
    resolver = getattr(module, "get_new_nuisance_target", None)
    if callable(resolver):
        target = resolver(legacy_name)
        if target is not None:
            return target

    direct = getattr(module, "LEGACY_TO_NEW_NUISANCE_MAP", {}) or {}
    if legacy_name in direct:
        return direct[legacy_name]

    multi = getattr(module, "LEGACY_TO_NEW_MULTI_MAP", {}) or {}
    if legacy_name in multi:
        return multi[legacy_name]

    shared = getattr(module, "COMPRESSED_OR_SHARED_MAP", {}) or {}
    if legacy_name in shared:
        return shared[legacy_name]

    return None


def is_explicitly_unmatched(module: Any, legacy_name: str) -> bool:
    predicate = getattr(module, "is_unmatched_legacy_nuisance", None)
    if callable(predicate):
        try:
            return bool(predicate(legacy_name))
        except Exception:
            pass

    unmatched = getattr(module, "UNMATCHED_OR_NOT_EXPOSED_IN_NEW", []) or []
    return legacy_name in set(unmatched)


def normalize_targets(target: Union[NewTarget, Sequence[NewTarget]]) -> NewTargets:
    if isinstance(target, tuple) and len(target) == 2 and all(isinstance(x, str) for x in target):
        return [target]

    if isinstance(target, Sequence):
        out: NewTargets = []
        for item in target:
            if not (isinstance(item, tuple) and len(item) == 2 and all(isinstance(x, str) for x in item)):
                raise TypeError(f"Invalid target entry: {item!r}")
            out.append(item)
        return out

    raise TypeError(f"Unsupported target type: {type(target)!r}")


def canonical_target(block: str, nuisance_id: str) -> CanonicalTarget:
    return (str(block), str(nuisance_id))


def canonical_pair(t1: CanonicalTarget, t2: CanonicalTarget) -> CanonicalPair:
    return tuple(sorted((t1, t2)))  # type: ignore[return-value]


def record_from_pair(pair: CanonicalPair, stat_corr: float, syst_corr: float) -> Dict[str, Any]:
    return {
        "block_1": pair[0][0],
        "id_1": pair[0][1],
        "block_2": pair[1][0],
        "id_2": pair[1][1],
        "stat_correlation": stat_corr,
        "syst_correlation": syst_corr,
    }


def records_equal(lhs: Dict[str, Any], rhs: Dict[str, Any]) -> bool:
    keys = set(lhs) | set(rhs)
    for key in keys:
        lv = lhs.get(key)
        rv = rhs.get(key)
        if isinstance(lv, (int, float)) and isinstance(rv, (int, float)):
            if not math.isclose(float(lv), float(rv), rel_tol=0.0, abs_tol=1e-15):
                return False
        else:
            if lv != rv:
                return False
    return True


def load_existing_output(path: Optional[Union[str, Path]]) -> Dict[CanonicalPair, Dict[str, Any]]:
    if path is None:
        return {}

    with Path(path).open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"{path}: existing JSON is not an object")
    corr = data.get("correlations", [])
    if not isinstance(corr, list):
        raise ValueError(f"{path}: top-level 'correlations' must be a list")

    out: Dict[CanonicalPair, Dict[str, Any]] = {}
    for idx, item in enumerate(corr):
        if not isinstance(item, dict):
            raise ValueError(f"{path}: correlations[{idx}] must be an object")
        required = {"block_1", "id_1", "block_2", "id_2", "stat_correlation", "syst_correlation"}
        missing = required - set(item)
        if missing:
            raise ValueError(f"{path}: correlations[{idx}] missing keys: {sorted(missing)}")

        t1 = canonical_target(item["block_1"], item["id_1"])
        t2 = canonical_target(item["block_2"], item["id_2"])
        pair = canonical_pair(t1, t2)
        out[pair] = record_from_pair(pair, float(item["stat_correlation"]), float(item["syst_correlation"]))

    return out


def make_report(
    *,
    parsed_lines: int,
    expanded_pairs: int,
    written_pairs: int,
    multi_target_self_pairs_added: int,
    skipped: List[SkippedCorrelation],
    identical_duplicates: List[DuplicateCorrelation],
    conflicting_duplicates: List[DuplicateCorrelation],
) -> Dict[str, Any]:
    return {
        "summary": {
            "legacy_lines_total": parsed_lines,
            "expanded_pairs": expanded_pairs,
            "written_pairs": written_pairs,
            "multi_target_self_pairs_added": multi_target_self_pairs_added,
            "skipped": len(skipped),
            "conflicting_duplicates": len(conflicting_duplicates),
            "identical_duplicates": len(identical_duplicates),
        },
        "skipped_entries": [asdict(item) for item in skipped],
        "conflicting_duplicates": [asdict(item) for item in conflicting_duplicates],
        "identical_duplicates": [asdict(item) for item in identical_duplicates],
    }


def inject_multi_target_self_correlations(
    module: Any,
    existing: Dict[CanonicalPair, Dict[str, Any]],
    conflicting_duplicates: List[DuplicateCorrelation],
) -> int:
    """
    If one legacy nuisance is split into several new nuisance entries, those new
    entries represent the same underlying random variable. Their mutual
    correlation must therefore be +1 unless an existing/explicit entry says
    something else, in which case we report a conflict.

    Example:
        old a0_HO_gpp_LbLll -> new Lb_L:1_3_0 and Lb_L:1_6_0
        => add corr(Lb_L:1_3_0, Lb_L:1_6_0) = +1
    """
    added = 0
    multi_map = getattr(module, "LEGACY_TO_NEW_MULTI_MAP", {}) or {}

    for legacy_name, raw_targets in multi_map.items():
        targets = normalize_targets(raw_targets)

        for (block1, id1), (block2, id2) in combinations(targets, 2):
            t1 = canonical_target(block1, id1)
            t2 = canonical_target(block2, id2)
            pair = canonical_pair(t1, t2)

            if pair[0] == pair[1]:
                continue

            record = record_from_pair(pair, 1.0, 0.0)
            previous = existing.get(pair)

            if previous is None:
                existing[pair] = record
                added += 1
                continue

            if not records_equal(previous, record):
                conflicting_duplicates.append(
                    DuplicateCorrelation(
                        line_no=0,
                        raw_name_1=legacy_name,
                        raw_name_2=legacy_name,
                        block_1=pair[0][0],
                        id_1=pair[0][1],
                        block_2=pair[1][0],
                        id_2=pair[1][1],
                        reason="multi_target_self_correlation_conflict",
                        previous=previous,
                        new=record,
                    )
                )

    return added


def convert(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    maps_path: Union[str, Path],
    report_path: Optional[Union[str, Path]] = None,
    merge_with: Optional[Union[str, Path]] = None,
    indent: int = 2,
) -> Tuple[int, Dict[str, Any]]:
    module = load_module(maps_path)

    existing = load_existing_output(merge_with)
    skipped: List[SkippedCorrelation] = []
    identical_duplicates: List[DuplicateCorrelation] = []
    conflicting_duplicates: List[DuplicateCorrelation] = []

    multi_target_self_pairs_added = inject_multi_target_self_correlations(
        module,
        existing,
        conflicting_duplicates,
    )

    parsed_lines = 0
    expanded_pairs = 0
    written_pairs = len(existing)

    for entry in iter_legacy_correlations(input_path):
        parsed_lines += 1

        target1 = resolve_target(module, entry.raw_name_1)
        target2 = resolve_target(module, entry.raw_name_2)

        if target1 is None:
            reason = "explicitly_unmatched_1" if is_explicitly_unmatched(module, entry.raw_name_1) else "unknown_nuisance_1"
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_name_1=entry.raw_name_1,
                    raw_name_2=entry.raw_name_2,
                    reason=reason,
                    details={
                        "stat_correlation": entry.stat_correlation,
                        "syst_correlation": entry.syst_correlation,
                    },
                )
            )
            continue

        if target2 is None:
            reason = "explicitly_unmatched_2" if is_explicitly_unmatched(module, entry.raw_name_2) else "unknown_nuisance_2"
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_name_1=entry.raw_name_1,
                    raw_name_2=entry.raw_name_2,
                    reason=reason,
                    details={
                        "stat_correlation": entry.stat_correlation,
                        "syst_correlation": entry.syst_correlation,
                    },
                )
            )
            continue

        try:
            targets1 = normalize_targets(target1)
            targets2 = normalize_targets(target2)
        except Exception as exc:
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_name_1=entry.raw_name_1,
                    raw_name_2=entry.raw_name_2,
                    reason="invalid_target",
                    details={"error": str(exc), "target1": repr(target1), "target2": repr(target2)},
                )
            )
            continue

        seen_pairs_this_line: set[CanonicalPair] = set()
        for (block1, id1), (block2, id2) in product(targets1, targets2):
            t1 = canonical_target(block1, id1)
            t2 = canonical_target(block2, id2)
            pair = canonical_pair(t1, t2)

            if pair[0] == pair[1]:
                skipped.append(
                    SkippedCorrelation(
                        line_no=entry.line_no,
                        raw_name_1=entry.raw_name_1,
                        raw_name_2=entry.raw_name_2,
                        reason="collapsed_to_same_target",
                        details={
                            "target": {"block": pair[0][0], "id": pair[0][1]},
                            "stat_correlation": entry.stat_correlation,
                            "syst_correlation": entry.syst_correlation,
                        },
                    )
                )
                continue

            if pair in seen_pairs_this_line:
                continue
            seen_pairs_this_line.add(pair)
            expanded_pairs += 1

            record = record_from_pair(pair, entry.stat_correlation, entry.syst_correlation)
            previous = existing.get(pair)
            if previous is None:
                existing[pair] = record
                written_pairs += 1
                continue

            dup = DuplicateCorrelation(
                line_no=entry.line_no,
                raw_name_1=entry.raw_name_1,
                raw_name_2=entry.raw_name_2,
                block_1=pair[0][0],
                id_1=pair[0][1],
                block_2=pair[1][0],
                id_2=pair[1][1],
                reason="identical_duplicate" if records_equal(previous, record) else "conflicting_duplicate",
                previous=previous,
                new=record,
            )
            if records_equal(previous, record):
                identical_duplicates.append(dup)
            else:
                conflicting_duplicates.append(dup)

    correlations_out = [existing[pair] for pair in sorted(existing.keys())]
    output_obj = {"correlations": correlations_out}

    output_path = Path(output_path)
    output_path.write_text(json.dumps(output_obj, indent=indent) + "\n", encoding="utf-8")

    report = make_report(
        parsed_lines=parsed_lines,
        expanded_pairs=expanded_pairs,
        written_pairs=written_pairs,
        multi_target_self_pairs_added=multi_target_self_pairs_added,
        skipped=skipped,
        identical_duplicates=identical_duplicates,
        conflicting_duplicates=conflicting_duplicates,
    )

    if report_path is not None:
        Path(report_path).write_text(json.dumps(report, indent=indent) + "\n", encoding="utf-8")

    summary = report["summary"]
    print(f"[OK] Wrote {output_path}")
    print(f"  total legacy lines     : {summary['legacy_lines_total']}")
    print(f"  expanded pairs         : {summary['expanded_pairs']}")
    print(f"  written pairs          : {summary['written_pairs']}")
    print(f"  multi self pairs added : {summary['multi_target_self_pairs_added']}")
    print(f"  skipped                : {summary['skipped']}")
    print(f"  conflicting duplicates : {summary['conflicting_duplicates']}")
    print(f"  identical duplicates   : {summary['identical_duplicates']}")

    if skipped:
        print("\nSkipped entries:")
        for item in skipped[:30]:
            print(f"  line {item.line_no}: {item.raw_name_1} <-> {item.raw_name_2} -> {item.reason}")
        if len(skipped) > 30:
            print(f"  ... and {len(skipped) - 30} more")

    if conflicting_duplicates:
        print("\nConflicting duplicates:")
        for item in conflicting_duplicates[:20]:
            print(
                f"  line {item.line_no}: {item.raw_name_1} <-> {item.raw_name_2} -> "
                f"{item.block_1}/{item.id_1} <-> {item.block_2}/{item.id_2}"
            )
        if len(conflicting_duplicates) > 20:
            print(f"  ... and {len(conflicting_duplicates) - 20} more")

    return (0 if not conflicting_duplicates else 2), report


def main() -> int:
    args = parse_args()
    rc, report = convert(
        input_path=args.input,
        output_path=args.output,
        maps_path=args.maps,
        report_path=args.report_json,
        merge_with=args.merge_with,
        indent=args.indent,
    )

    if args.fail_on_conflict and report["summary"]["conflicting_duplicates"]:
        rc = max(rc, 2)
    if args.fail_on_skipped and report["summary"]["skipped"]:
        rc = max(rc, 3)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())