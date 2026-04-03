#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, List, Optional, Sequence, Tuple, Union

NewTarget = Tuple[str, str]
NewTargets = List[NewTarget]


@dataclass
class LegacyEntry:
    line_no: int
    raw_name: str
    central_value: float
    stat_error: float
    nuisance_kind: Optional[int]
    raw_line: str


@dataclass
class SkippedEntry:
    line_no: int
    raw_name: str
    reason: str
    details: Dict[str, Any]


@dataclass
class DuplicateEntry:
    line_no: int
    raw_name: str
    target_block: str
    target_id: str
    reason: str
    previous: Dict[str, Any]
    new: Dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Legacy nuisance input text file")
    parser.add_argument("output", help="Output JSON file in the new nuisance format")
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
        "--merge-existing",
        default=None,
        help=(
            "Optional existing new-format JSON to merge into before writing. "
            "If omitted, the output is built from scratch."
        ),
    )
    parser.add_argument(
        "--merge-output-if-exists",
        action="store_true",
        help=(
            "Merge into the output JSON itself if it already exists. "
            "Ignored when --merge-existing is provided."
        ),
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level for output files (default: 2)",
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


_FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def _looks_like_float(token: str) -> bool:
    return bool(_FLOAT_RE.match(token))


def parse_legacy_line(line: str, line_no: int) -> Optional[LegacyEntry]:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None

    parts = stripped.split()
    if len(parts) < 3:
        raise ValueError(f"line {line_no}: expected at least 3 columns, got {len(parts)}")

    name = parts[0]
    if not _looks_like_float(parts[1]) or not _looks_like_float(parts[2]):
        raise ValueError(
            f"line {line_no}: could not parse central/stat values from columns 2/3: {parts[1:3]}"
        )

    central_value = float(parts[1])
    stat_error = float(parts[2])
    nuisance_kind: Optional[int] = None
    if len(parts) >= 4:
        try:
            nuisance_kind = int(float(parts[3]))
        except ValueError:
            nuisance_kind = None

    return LegacyEntry(
        line_no=line_no,
        raw_name=name,
        central_value=central_value,
        stat_error=stat_error,
        nuisance_kind=nuisance_kind,
        raw_line=line.rstrip("\n"),
    )


def iter_legacy_entries(path: Union[str, Path]) -> Iterable[LegacyEntry]:
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


def build_payload(entry: LegacyEntry) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "central_value": entry.central_value,
        "stat_error": entry.stat_error,
        "syst_error": 0,
    }
    return payload


def payloads_equal(lhs: Dict[str, Any], rhs: Dict[str, Any]) -> bool:
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


def load_existing_output(path: Union[str, Path]) -> Dict[str, Dict[str, Any]]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Existing JSON file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise TypeError(f"Existing JSON root must be an object: {path}")

    normalized: Dict[str, Dict[str, Any]] = {}
    for block, entries in data.items():
        if not isinstance(block, str):
            raise TypeError(f"Block name must be a string in {path}: {block!r}")
        if not isinstance(entries, dict):
            raise TypeError(f"Block {block!r} must contain an object in {path}")

        normalized_block: Dict[str, Any] = {}
        for nuisance_id, payload in entries.items():
            if not isinstance(nuisance_id, str):
                raise TypeError(
                    f"Nuisance id must be a string in block {block!r} from {path}: {nuisance_id!r}"
                )
            if not isinstance(payload, dict):
                raise TypeError(
                    f"Payload for {block}/{nuisance_id} must be an object in {path}"
                )
            normalized_block[nuisance_id] = dict(payload)
        normalized[block] = normalized_block

    return normalized


def choose_merge_source(
    *,
    explicit_merge_path: Optional[Union[str, Path]],
    merge_output_if_exists: bool,
    output_path: Union[str, Path],
) -> Optional[Path]:
    if explicit_merge_path is not None:
        return Path(explicit_merge_path)

    output_path = Path(output_path)
    if merge_output_if_exists and output_path.exists():
        return output_path

    return None


def make_report(
    *,
    parsed_lines: int,
    converted_lines: int,
    written_entries: int,
    preloaded_entries: int,
    skipped: List[SkippedEntry],
    identical_duplicates: List[DuplicateEntry],
    conflicting_duplicates: List[DuplicateEntry],
    output_data: Dict[str, Dict[str, Any]],
    merge_source: Optional[str],
) -> Dict[str, Any]:
    return {
        "summary": {
            "legacy_lines_total": parsed_lines,
            "converted_lines": converted_lines,
            "written_entries": written_entries,
            "preloaded_entries": preloaded_entries,
            "final_entries": sum(len(v) for v in output_data.values()),
            "skipped": len(skipped),
            "conflicting_duplicates": len(conflicting_duplicates),
            "identical_duplicates": len(identical_duplicates),
            "blocks_written": sorted(output_data.keys()),
            "merge_source": merge_source,
        },
        "skipped_entries": [asdict(item) for item in skipped],
        "conflicting_duplicates": [asdict(item) for item in conflicting_duplicates],
        "identical_duplicates": [asdict(item) for item in identical_duplicates],
    }


def convert(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    maps_path: Union[str, Path],
    report_path: Optional[Union[str, Path]] = None,
    merge_existing_path: Optional[Union[str, Path]] = None,
    merge_output_if_exists: bool = False,
    indent: int = 2,
) -> int:
    module = load_module(maps_path)

    merge_source_path = choose_merge_source(
        explicit_merge_path=merge_existing_path,
        merge_output_if_exists=merge_output_if_exists,
        output_path=output_path,
    )

    output_data: DefaultDict[str, Dict[str, Any]] = defaultdict(dict)
    preloaded_entries = 0

    if merge_source_path is not None:
        existing_data = load_existing_output(merge_source_path)
        for block, entries in existing_data.items():
            output_data[block].update(entries)
            preloaded_entries += len(entries)

    skipped: List[SkippedEntry] = []
    identical_duplicates: List[DuplicateEntry] = []
    conflicting_duplicates: List[DuplicateEntry] = []

    parsed_lines = 0
    converted_lines = 0
    written_entries = 0

    for entry in iter_legacy_entries(input_path):
        parsed_lines += 1

        target = resolve_target(module, entry.raw_name)
        if target is None:
            reason = "explicitly_unmatched" if is_explicitly_unmatched(module, entry.raw_name) else "unknown_nuisance"
            skipped.append(
                SkippedEntry(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    reason=reason,
                    details={
                        "central_value": entry.central_value,
                        "stat_error": entry.stat_error,
                        "nuisance_kind": entry.nuisance_kind,
                    },
                )
            )
            continue

        try:
            targets = normalize_targets(target)
        except Exception as exc:
            skipped.append(
                SkippedEntry(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    reason="invalid_target",
                    details={"error": str(exc), "target": repr(target)},
                )
            )
            continue

        payload = build_payload(entry)
        converted_lines += 1

        for block, nuisance_id in targets:
            existing = output_data[block].get(nuisance_id)
            if existing is None:
                output_data[block][nuisance_id] = payload
                written_entries += 1
                continue

            dup = DuplicateEntry(
                line_no=entry.line_no,
                raw_name=entry.raw_name,
                target_block=block,
                target_id=nuisance_id,
                reason="identical_duplicate" if payloads_equal(existing, payload) else "conflicting_duplicate",
                previous=existing,
                new=payload,
            )
            if payloads_equal(existing, payload):
                identical_duplicates.append(dup)
            else:
                conflicting_duplicates.append(dup)

    sorted_output_data = {
        block: dict(sorted(entries.items(), key=lambda item: item[0]))
        for block, entries in sorted(output_data.items(), key=lambda item: item[0])
    }

    output_path = Path(output_path)
    output_path.write_text(json.dumps(sorted_output_data, indent=indent) + "\n", encoding="utf-8")

    report = make_report(
        parsed_lines=parsed_lines,
        converted_lines=converted_lines,
        written_entries=written_entries,
        preloaded_entries=preloaded_entries,
        skipped=skipped,
        identical_duplicates=identical_duplicates,
        conflicting_duplicates=conflicting_duplicates,
        output_data=sorted_output_data,
        merge_source=str(merge_source_path) if merge_source_path is not None else None,
    )

    if report_path is not None:
        Path(report_path).write_text(json.dumps(report, indent=indent) + "\n", encoding="utf-8")

    print(f"[OK] Wrote {output_path}")
    print(f"  total legacy lines     : {parsed_lines}")
    print(f"  converted lines        : {converted_lines}")
    print(f"  preloaded entries      : {preloaded_entries}")
    print(f"  written entries        : {written_entries}")
    print(f"  final entries          : {sum(len(v) for v in sorted_output_data.values())}")
    print(f"  skipped                : {len(skipped)}")
    print(f"  conflicting duplicates : {len(conflicting_duplicates)}")
    print(f"  identical duplicates   : {len(identical_duplicates)}")
    print(f"  merge source           : {merge_source_path if merge_source_path is not None else '(none)'}")
    print(f"  blocks written         : {', '.join(sorted(sorted_output_data.keys())) if sorted_output_data else '(none)'}")

    if skipped:
        print("\nSkipped entries:")
        for item in skipped:
            print(f"  line {item.line_no}: {item.raw_name} -> {item.reason}")

    if conflicting_duplicates:
        print("\nConflicting duplicates:")
        for item in conflicting_duplicates:
            print(
                f"  line {item.line_no}: {item.raw_name} -> {item.target_block}/{item.target_id}"
            )

    return 0 if not conflicting_duplicates else 2


def main() -> int:
    args = parse_args()
    return convert(
        input_path=args.input,
        output_path=args.output,
        maps_path=args.maps,
        report_path=args.report_json,
        merge_existing_path=args.merge_existing,
        merge_output_if_exists=args.merge_output_if_exists,
        indent=args.indent,
    )


if __name__ == "__main__":
    raise SystemExit(main())
