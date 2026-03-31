#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

NUM_RE = re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')


def is_number_token(s: str) -> bool:
    return bool(NUM_RE.fullmatch(s.strip()))


def norm_zero(x: float) -> float:
    return 0.0 if x == 0.0 else x


def trim_decimal_string(s: str) -> str:
    if "." not in s:
        return s
    s = s.rstrip("0").rstrip(".")
    return s if s else "0"


@dataclass(frozen=True)
class EncodedBin:
    int_part: int
    frac_part: int
    frac_digits: int


def encode_bin_gev(x: float) -> EncodedBin:
    if not math.isfinite(x):
        raise ValueError(f"Bin value must be finite, got {x!r}")

    x = norm_zero(x)
    k_max_frac_digits = 12
    k_scale = 1.0e12

    xr = round(float(x) * k_scale) / k_scale
    neg = math.copysign(1.0, xr) < 0.0
    ax = abs(xr)

    s = f"{ax:.{k_max_frac_digits}f}"
    s = trim_decimal_string(s)

    if "." in s:
        int_str, frac_str = s.split(".", 1)
    else:
        int_str, frac_str = s, ""

    if not int_str:
        int_str = "0"

    int_part = int(int_str)
    frac_part = int(frac_str) if frac_str else 0
    frac_digits = len(frac_str)

    if neg:
        if int_part != 0:
            return EncodedBin(-int_part, frac_part, frac_digits)
        if frac_part != 0:
            return EncodedBin(0, -frac_part, frac_digits)

    return EncodedBin(int_part, frac_part, frac_digits)


def insert_bin_parts(unbinned_parts: Sequence[int], low: float, high: float) -> List[int]:
    low_e = encode_bin_gev(low)
    high_e = encode_bin_gev(high)

    bin_parts = [
        low_e.int_part,
        low_e.frac_part,
        low_e.frac_digits,
        high_e.int_part,
        high_e.frac_part,
        high_e.frac_digits,
    ]
    parts = list(unbinned_parts)
    return parts[:2] + bin_parts + parts[2:]


def parts_to_key(parts: Sequence[int]) -> str:
    return "_".join(str(x) for x in parts)


def canonical_pair(key1: str, key2: str) -> Tuple[str, str]:
    return tuple(sorted((key1, key2)))


def canonicalize_observable_name(s: str) -> str:
    s = s.strip()

    replacements = [
        ("dGamma/dq2", "dGamma_dq2"),
        ("dBR/dq2", "dBR_dq2"),
        ("K0*", "Kstar0"),
        ("K*0", "Kstar0"),
        ("K*", "Kstar"),
        ("D0*", "Dstar0"),
        ("D*0", "Dstar0"),
        ("D*", "Dstar"),
        ("Lambda_b", "Lambdab"),
        ("P'_", "Pprime_"),
        ("P'", "Pprime"),
        ("ATReCP", "A_T_RE_CPV"),
        ("ATImCP", "A_T_IM_CPV"),
        ("ATRe", "A_T_RE"),
        ("ATIm", "A_T_IM"),
        ("AFBCP", "A_FB_CPV"),
        ("|epsilon_K|", "absepsilonk"),
    ]

    for old, new in replacements:
        s = s.replace(old, new)

    s = re.sub(r'\bP([4568])prime\b', r'Pprime_\1', s)
    s = re.sub(r'\bP([4568])prime_', r'Pprime_\1_', s)
    s = s.lower()
    s = re.sub(r'[^a-z0-9]+', '', s)
    return s


@dataclass
class ParsedLegacyName:
    raw_name: str
    base_name: str
    section: str
    has_bin: bool
    low: float
    high: float


@dataclass
class ResolvedObservable:
    raw_name: str
    base_name: str
    section: str
    has_bin: bool
    low: float
    high: float
    enum_name: str
    canonical_name: Optional[str]
    key: str
    parts: List[int]


@dataclass
class LegacyCorrelation:
    line_no: int
    raw_obs1: str
    raw_obs2: str
    stat_correlation: float


@dataclass
class ConvertedCorrelation:
    line_no: int
    raw_obs1: str
    raw_obs2: str
    section: str
    key1: str
    key2: str
    pair_key: Tuple[str, str]
    stat_correlation: float
    syst_correlation: Optional[float]
    obs1: Dict[str, Any]
    obs2: Dict[str, Any]


@dataclass
class SkippedCorrelation:
    line_no: int
    raw_obs1: str
    raw_obs2: str
    reason: str
    details: Dict[str, Any]


def load_python_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load python module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_alias_map(
    observable_mapping: Dict[str, str],
    observable_flha_mapping: Dict[str, Sequence[int]],
    manual_aliases: Dict[str, str],
) -> Dict[str, str]:
    alias_to_enum: Dict[str, str] = {}

    for enum_name, canon_name in observable_mapping.items():
        if enum_name not in observable_flha_mapping:
            continue
        alias_to_enum[canonicalize_observable_name(canon_name)] = enum_name

    for legacy_name, enum_name in manual_aliases.items():
        if enum_name not in observable_flha_mapping:
            continue
        alias_to_enum[canonicalize_observable_name(legacy_name)] = enum_name

    return alias_to_enum


def parse_legacy_correlations(path: Path) -> List[LegacyCorrelation]:
    out: List[LegacyCorrelation] = []

    with path.open("r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue

            cols = stripped.split()
            if len(cols) < 3:
                raise ValueError(
                    f"{path}:{line_no}: expected at least 3 whitespace-separated columns, got: {raw.rstrip()}"
                )

            raw_obs1 = cols[0]
            raw_obs2 = cols[1]
            try:
                stat_corr = float(cols[2])
            except Exception as exc:
                raise ValueError(f"{path}:{line_no}: cannot parse correlation value in {raw.rstrip()}") from exc

            out.append(
                LegacyCorrelation(
                    line_no=line_no,
                    raw_obs1=raw_obs1,
                    raw_obs2=raw_obs2,
                    stat_correlation=stat_corr,
                )
            )

    return out


def candidate_parses(raw_name: str) -> List[ParsedLegacyName]:
    tokens = raw_name.split("_")
    candidates: List[ParsedLegacyName] = []

    def add(base_tokens: Sequence[str], section: str, has_bin: bool, low: float, high: float) -> None:
        base_name = "_".join(base_tokens).strip()
        if not base_name:
            return
        cand = ParsedLegacyName(
            raw_name=raw_name,
            base_name=base_name,
            section=section or "DEFAULT",
            has_bin=has_bin,
            low=low,
            high=high,
        )
        if cand not in candidates:
            candidates.append(cand)

    if len(tokens) >= 4 and (not is_number_token(tokens[-1])) and is_number_token(tokens[-2]) and is_number_token(tokens[-3]):
        add(tokens[:-3], tokens[-1], True, float(tokens[-3]), float(tokens[-2]))

    if len(tokens) >= 3 and is_number_token(tokens[-1]) and is_number_token(tokens[-2]):
        add(tokens[:-2], "DEFAULT", True, float(tokens[-2]), float(tokens[-1]))

    if len(tokens) >= 2 and (not is_number_token(tokens[-1])):
        add(tokens[:-1], tokens[-1], False, 0.0, 0.0)

    add(tokens, "DEFAULT", False, 0.0, 0.0)
    return candidates


def resolve_legacy_observable(
    raw_name: str,
    alias_map: Dict[str, str],
    observable_mapping: Dict[str, str],
    observable_flha_mapping: Dict[str, Sequence[int]],
) -> Tuple[Optional[ResolvedObservable], Dict[str, Any]]:
    attempts = []

    for cand in candidate_parses(raw_name):
        norm_base = canonicalize_observable_name(cand.base_name)
        enum_name = alias_map.get(norm_base)

        attempts.append(
            {
                "base_name": cand.base_name,
                "section": cand.section,
                "has_bin": cand.has_bin,
                "low": cand.low,
                "high": cand.high,
                "normalized_base": norm_base,
                "enum_name": enum_name,
            }
        )

        if enum_name is None:
            continue
        if enum_name not in observable_flha_mapping:
            return None, {"reason": "missing_flha_mapping", "enum_name": enum_name, "attempts": attempts}

        parts = insert_bin_parts(observable_flha_mapping[enum_name], cand.low, cand.high)
        return (
            ResolvedObservable(
                raw_name=raw_name,
                base_name=cand.base_name,
                section=cand.section,
                has_bin=cand.has_bin,
                low=cand.low,
                high=cand.high,
                enum_name=enum_name,
                canonical_name=observable_mapping.get(enum_name),
                key=parts_to_key(parts),
                parts=parts,
            ),
            {"attempts": attempts},
        )

    return None, {"reason": "unknown_legacy_name", "attempts": attempts}


def ensure_output_root(merge_with: Optional[Path]) -> Dict[str, Any]:
    if merge_with is None:
        return {"correlations": {}}

    with merge_with.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError(f"{merge_with}: existing JSON is not an object")
    if "correlations" not in data:
        data["correlations"] = {}
    if not isinstance(data["correlations"], dict):
        raise ValueError(f"{merge_with}: top-level correlations is not an object")
    return data


def canonicalize_existing_correlation_lists(data: Dict[str, Any]) -> Dict[str, Dict[Tuple[str, str], Dict[str, Any]]]:
    """
    Convert existing correlations lists into section -> canonical_pair -> record.
    Also normalizes pair ordering in memory.
    """
    out: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]] = {}
    corr = data["correlations"]

    for section, items in corr.items():
        if not isinstance(items, list):
            raise ValueError(f"Existing correlations[{section!r}] must be a list")

        sec_dict: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for idx, item in enumerate(items):
            if not isinstance(item, dict):
                raise ValueError(f"Existing correlations[{section!r}][{idx}] must be an object")
            if "id_1" not in item or "id_2" not in item or "stat_correlation" not in item:
                raise ValueError(f"Existing correlations[{section!r}][{idx}] missing id_1/id_2/stat_correlation")

            key1 = str(item["id_1"])
            key2 = str(item["id_2"])
            pair = canonical_pair(key1, key2)
            sec_dict[pair] = {
                **item,
                "id_1": pair[0],
                "id_2": pair[1],
            }
        out[section] = sec_dict

    return out


def make_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert legacy correlation text format into the new JSON correlations format.")
    p.add_argument("legacy_file", type=Path, help="Old flat correlation text file")
    p.add_argument("output_json", type=Path, help="Output JSON file")
    p.add_argument(
        "--maps",
        type=Path,
        default=Path("observable_maps_fixed.py"),
        help="Python module containing OBSERVABLE_MAPPING / OBSERVABLE_FLHA_MAPPING / optional MANUAL_LEGACY_ALIASES",
    )
    p.add_argument(
        "--merge-with",
        type=Path,
        default=None,
        help="Optional existing JSON file to merge into before writing output",
    )
    p.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional JSON report with converted/skipped/conflicting entries",
    )
    p.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Indent level for output JSON (default: 2)",
    )
    p.add_argument(
        "--write-syst-zero",
        action="store_true",
        help="Also write syst_correlation = 0.0 for all converted entries",
    )
    p.add_argument(
        "--fail-on-conflict",
        action="store_true",
        help="Exit with non-zero status if the same output pair is produced with different values",
    )
    p.add_argument(
        "--fail-on-skipped",
        action="store_true",
        help="Exit with non-zero status if at least one legacy line could not be converted",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_cli().parse_args(argv)

    maps_mod = load_python_module(args.maps)
    observable_mapping = getattr(maps_mod, "OBSERVABLE_MAPPING")
    observable_flha_mapping = getattr(maps_mod, "OBSERVABLE_FLHA_MAPPING")
    manual_aliases = getattr(maps_mod, "MANUAL_LEGACY_ALIASES", {})

    if not isinstance(observable_mapping, dict) or not isinstance(observable_flha_mapping, dict):
        raise TypeError("maps module must define OBSERVABLE_MAPPING and OBSERVABLE_FLHA_MAPPING as dict")

    alias_map = build_alias_map(
        observable_mapping=observable_mapping,
        observable_flha_mapping=observable_flha_mapping,
        manual_aliases=manual_aliases,
    )

    legacy_corrs = parse_legacy_correlations(args.legacy_file)

    out = ensure_output_root(args.merge_with)
    existing = canonicalize_existing_correlation_lists(out)

    converted: List[ConvertedCorrelation] = []
    skipped: List[SkippedCorrelation] = []
    conflicts: List[Dict[str, Any]] = []
    duplicate_identical: List[Dict[str, Any]] = []

    for entry in legacy_corrs:
        obs1, meta1 = resolve_legacy_observable(
            entry.raw_obs1,
            alias_map=alias_map,
            observable_mapping=observable_mapping,
            observable_flha_mapping=observable_flha_mapping,
        )
        obs2, meta2 = resolve_legacy_observable(
            entry.raw_obs2,
            alias_map=alias_map,
            observable_mapping=observable_mapping,
            observable_flha_mapping=observable_flha_mapping,
        )

        if obs1 is None:
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    reason="unknown_obs1",
                    details=meta1,
                )
            )
            continue

        if obs2 is None:
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    reason="unknown_obs2",
                    details=meta2,
                )
            )
            continue

        if obs1.section != obs2.section:
            skipped.append(
                SkippedCorrelation(
                    line_no=entry.line_no,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    reason="inconsistent_legacy_sections",
                    details={"obs1": asdict(obs1), "obs2": asdict(obs2)},
                )
            )
            continue

        section = obs1.section
        pair = canonical_pair(obs1.key, obs2.key)
        record = {
            "id_1": pair[0],
            "id_2": pair[1],
            "stat_correlation": entry.stat_correlation,
        }
        if args.write_syst_zero:
            record["syst_correlation"] = 0.0

        sec_dict = existing.setdefault(section, {})
        prev = sec_dict.get(pair)
        if prev is not None:
            same = (
                float(prev.get("stat_correlation")) == float(record["stat_correlation"])
                and (
                    ("syst_correlation" not in prev and "syst_correlation" not in record)
                    or float(prev.get("syst_correlation", 0.0)) == float(record.get("syst_correlation", 0.0))
                )
            )
            if same:
                duplicate_identical.append(
                    {
                        "line_no": entry.line_no,
                        "raw_obs1": entry.raw_obs1,
                        "raw_obs2": entry.raw_obs2,
                        "section": section,
                        "pair_key": pair,
                    }
                )
            else:
                conflicts.append(
                    {
                        "line_no": entry.line_no,
                        "raw_obs1": entry.raw_obs1,
                        "raw_obs2": entry.raw_obs2,
                        "section": section,
                        "pair_key": pair,
                        "previous": prev,
                        "new": record,
                    }
                )

        sec_dict[pair] = record

        converted.append(
            ConvertedCorrelation(
                line_no=entry.line_no,
                raw_obs1=entry.raw_obs1,
                raw_obs2=entry.raw_obs2,
                section=section,
                key1=pair[0],
                key2=pair[1],
                pair_key=pair,
                stat_correlation=entry.stat_correlation,
                syst_correlation=(0.0 if args.write_syst_zero else None),
                obs1=asdict(obs1),
                obs2=asdict(obs2),
            )
        )

    # Convert canonicalized dicts back to JSON lists
    out["correlations"] = {}
    for section in sorted(existing.keys()):
        records = list(existing[section].values())
        records.sort(key=lambda r: (str(r["id_1"]), str(r["id_2"])))
        out["correlations"][section] = records

    with args.output_json.open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=args.indent, sort_keys=True)
        fh.write("\n")

    summary = {
        "legacy_lines_total": len(legacy_corrs),
        "converted_count": len(converted),
        "skipped_count": len(skipped),
        "conflict_count": len(conflicts),
        "duplicate_identical_count": len(duplicate_identical),
        "sections_written": sorted(out["correlations"].keys()),
    }

    print(f"[OK] Wrote {args.output_json}")
    print(f"  total legacy lines       : {summary['legacy_lines_total']}")
    print(f"  converted                : {summary['converted_count']}")
    print(f"  skipped                  : {summary['skipped_count']}")
    print(f"  conflicting duplicates   : {summary['conflict_count']}")
    print(f"  identical duplicates     : {summary['duplicate_identical_count']}")
    print(f"  sections written         : {', '.join(summary['sections_written']) if summary['sections_written'] else '(none)'}")

    if skipped:
        print("\nSkipped entries:")
        for s in skipped[:30]:
            print(f"  line {s.line_no}: {s.raw_obs1} <-> {s.raw_obs2} -> {s.reason}")
        if len(skipped) > 30:
            print(f"  ... and {len(skipped) - 30} more")

    if conflicts:
        print("\nConflicting duplicates:")
        for c in conflicts[:20]:
            print(f"  line {c['line_no']}: {c['raw_obs1']} <-> {c['raw_obs2']} -> {c['section']} / {c['pair_key']}")
        if len(conflicts) > 20:
            print(f"  ... and {len(conflicts) - 20} more")

    if args.report_json is not None:
        report = {
            "summary": summary,
            "skipped": [asdict(x) for x in skipped],
            "converted": [asdict(x) for x in converted],
            "conflicts": conflicts,
            "duplicate_identical": duplicate_identical,
        }
        with args.report_json.open("w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"[OK] Wrote report {args.report_json}")

    rc = 0
    if args.fail_on_conflict and conflicts:
        rc = 2
    if args.fail_on_skipped and skipped:
        rc = max(rc, 3)

    return rc


if __name__ == "__main__":
    sys.exit(main())
