from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

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


def decode_bin_gev(int_part: int, frac_part: int, frac_digits: int) -> float:
    if frac_digits < 0:
        raise ValueError("decode_bin_gev: negative frac_digits")

    neg = (int_part < 0) or (frac_part < 0)
    abs_int = abs(int_part)
    abs_frac = abs(frac_part)

    out = float(abs_int)
    if frac_digits > 0:
        out += abs_frac / (10 ** frac_digits)

    if neg:
        out = -out

    return norm_zero(out)


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


def parse_key_parts(key: str) -> List[int]:
    try:
        return [int(x) for x in key.split("_")]
    except Exception as exc:
        raise ValueError(f"Invalid FLHA key {key!r}") from exc


def unbinned_parts_from_binned(full_parts: Sequence[int]) -> Optional[List[int]]:
    if len(full_parts) < 8:
        return None
    return list(full_parts[:2]) + list(full_parts[8:])


def decode_bins_from_binned(full_parts: Sequence[int]) -> Optional[Tuple[float, float]]:
    if len(full_parts) < 8:
        return None
    low = decode_bin_gev(full_parts[2], full_parts[3], full_parts[4])
    high = decode_bin_gev(full_parts[5], full_parts[6], full_parts[7])
    return (low, high)


def canonicalize_observable_name(s: str) -> str:
    """
    Normalize old and canonical mapping names to a comparable token stream.

    This function tries to absorb common legacy shorthand:
    - K* / Kstar, D* / Dstar
    - P4prime / P'_4 style variants
    - ATReCP / ATImCP / AFBCP spellings
    - punctuation / separators / case
    """
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
class LegacyEntry:
    line_no: int
    raw_name: str
    base_name: str
    section: str
    has_bin: bool
    low: float
    high: float
    central: float
    stat: float
    syst: float
    enum_name: Optional[str] = None
    canonical_name: Optional[str] = None
    expected_key: Optional[str] = None


@dataclass
class NewEntry:
    section: str
    key: str
    parts: List[int]
    unbinned_parts: Optional[List[int]]
    bins: Optional[Tuple[float, float]]
    central: float
    stat: float
    syst: float


@dataclass
class DiffRecord:
    line_no: int
    raw_name: str
    base_name: str
    section: str
    expected_key: Optional[str]
    reason: str
    legacy: Optional[Dict[str, float]]
    new: Optional[Dict[str, float]]
    details: Optional[Any] = None


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
) -> Tuple[Dict[str, str], List[Tuple[str, List[str]]]]:
    alias_to_enum: Dict[str, str] = {}
    collisions: Dict[str, List[str]] = {}

    for enum_name, canon_name in observable_mapping.items():
        if enum_name not in observable_flha_mapping:
            continue
        alias = canonicalize_observable_name(canon_name)
        prev = alias_to_enum.get(alias)
        if prev is not None and prev != enum_name:
            collisions.setdefault(alias, [prev]).append(enum_name)
        else:
            alias_to_enum[alias] = enum_name

    for legacy_name, enum_name in manual_aliases.items():
        alias = canonicalize_observable_name(legacy_name)
        prev = alias_to_enum.get(alias)
        if prev is not None and prev != enum_name:
            collisions.setdefault(alias, [prev]).append(enum_name)
        alias_to_enum[alias] = enum_name

    collision_list = []
    for alias, enums in collisions.items():
        uniq = []
        for e in enums:
            if e not in uniq:
                uniq.append(e)
        collision_list.append((alias, uniq))

    return alias_to_enum, collision_list


def parse_legacy_file(path: Path, valid_sections: Iterable[str]) -> List[LegacyEntry]:
    valid_sections = set(valid_sections)
    out: List[LegacyEntry] = []

    with path.open("r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue

            cols = stripped.split()
            if len(cols) < 4:
                raise ValueError(
                    f"{path}:{line_no}: expected at least 4 whitespace-separated columns, got: {raw.rstrip()}"
                )

            name = cols[0]
            try:
                central = float(cols[1])
                stat = float(cols[2])
                syst = float(cols[3])
            except Exception as exc:
                raise ValueError(f"{path}:{line_no}: cannot parse numeric columns in {raw.rstrip()}") from exc

            section = "DEFAULT"
            core = name
            tokens = core.split("_")
            if tokens and tokens[-1] in valid_sections and tokens[-1] != "DEFAULT":
                section = tokens[-1]
                core = "_".join(tokens[:-1])
                tokens = core.split("_")

            has_bin = len(tokens) >= 3 and is_number_token(tokens[-1]) and is_number_token(tokens[-2])
            if has_bin:
                low = float(tokens[-2])
                high = float(tokens[-1])
                base_name = "_".join(tokens[:-2])
            else:
                low = 0.0
                high = 0.0
                base_name = core

            out.append(
                LegacyEntry(
                    line_no=line_no,
                    raw_name=name,
                    base_name=base_name,
                    section=section,
                    has_bin=has_bin,
                    low=low,
                    high=high,
                    central=central,
                    stat=stat,
                    syst=syst,
                )
            )

    return out


def parse_new_json(path: Path) -> Tuple[Dict[str, Dict[str, NewEntry]], Dict[str, Dict[Tuple[int, ...], List[NewEntry]]]]:
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    fobs = data.get("FOBS")
    if not isinstance(fobs, dict):
        raise ValueError(f"{path}: top-level key 'FOBS' missing or not an object")

    by_section: Dict[str, Dict[str, NewEntry]] = {}
    by_unbinned: Dict[str, Dict[Tuple[int, ...], List[NewEntry]]] = {}

    for section, payload in fobs.items():
        if not isinstance(payload, dict):
            raise ValueError(f"{path}: FOBS[{section!r}] is not an object")

        sec_entries: Dict[str, NewEntry] = {}
        sec_unbinned: Dict[Tuple[int, ...], List[NewEntry]] = {}

        for key, value in payload.items():
            if not isinstance(value, dict):
                raise ValueError(f"{path}: FOBS[{section!r}][{key!r}] is not an object")

            parts = parse_key_parts(key)
            unbinned = unbinned_parts_from_binned(parts)
            bins = decode_bins_from_binned(parts)

            entry = NewEntry(
                section=section,
                key=key,
                parts=parts,
                unbinned_parts=unbinned,
                bins=bins,
                central=float(value["central_value"]),
                stat=float(value["stat_error"]),
                syst=float(value["syst_error"]),
            )
            sec_entries[key] = entry
            if unbinned is not None:
                sec_unbinned.setdefault(tuple(unbinned), []).append(entry)

        by_section[section] = sec_entries
        by_unbinned[section] = sec_unbinned

    return by_section, by_unbinned


def compare_float(a: float, b: float, rel_tol: float, abs_tol: float) -> bool:
    return math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol)


def compare_entry(old: LegacyEntry, new: NewEntry, rel_tol: float, abs_tol: float) -> List[str]:
    diffs = []
    if not compare_float(old.central, new.central, rel_tol, abs_tol):
        diffs.append("central_value")
    if not compare_float(old.stat, new.stat, rel_tol, abs_tol):
        diffs.append("stat_error")
    if not compare_float(old.syst, new.syst, rel_tol, abs_tol):
        diffs.append("syst_error")
    return diffs


def fmt(x: float) -> str:
    return f"{x:.16g}"


def make_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Compare legacy observable inputs against the new JSON FOBS format.")
    p.add_argument("legacy_file", type=Path, help="Old flat text input file")
    p.add_argument("new_json", type=Path, help="New JSON file")
    p.add_argument(
        "--maps",
        type=Path,
        default=Path("observable_maps.py"),
        help="Python module containing OBSERVABLE_MAPPING / OBSERVABLE_FLHA_MAPPING / optional MANUAL_LEGACY_ALIASES",
    )
    p.add_argument("--rel-tol", type=float, default=0.0, help="Relative tolerance for numerical comparison")
    p.add_argument("--abs-tol", type=float, default=0.0, help="Absolute tolerance for numerical comparison")
    p.add_argument("--json-report", type=Path, default=None, help="Optional JSON output report")
    p.add_argument(
        "--show-matched",
        action="store_true",
        help="Also print the fully matched entries",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_cli().parse_args(argv)

    maps_mod = load_python_module(args.maps)
    observable_mapping = getattr(maps_mod, "OBSERVABLE_MAPPING")
    observable_flha_mapping = getattr(maps_mod, "OBSERVABLE_FLHA_MAPPING")
    manual_aliases = getattr(maps_mod, "MANUAL_LEGACY_ALIASES", {})

    if not isinstance(observable_mapping, dict) or not isinstance(observable_flha_mapping, dict):
        raise TypeError("observable_maps.py must define OBSERVABLE_MAPPING and OBSERVABLE_FLHA_MAPPING as dict")

    alias_map, alias_collisions = build_alias_map(
        observable_mapping=observable_mapping,
        observable_flha_mapping=observable_flha_mapping,
        manual_aliases=manual_aliases,
    )

    new_by_section, new_by_unbinned = parse_new_json(args.new_json)
    legacy_entries = parse_legacy_file(args.legacy_file, valid_sections=new_by_section.keys())

    matched_keys = set()
    unknown_legacy: List[DiffRecord] = []
    missing_in_new: List[DiffRecord] = []
    value_mismatches: List[DiffRecord] = []
    fully_matched: List[LegacyEntry] = []
    missing_flha_mapping: List[DiffRecord] = []

    for entry in legacy_entries:
        norm_base = canonicalize_observable_name(entry.base_name)
        enum_name = alias_map.get(norm_base)

        if enum_name is None:
            unknown_legacy.append(
                DiffRecord(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    base_name=entry.base_name,
                    section=entry.section,
                    expected_key=None,
                    reason="unknown_legacy_name",
                    legacy={"central_value": entry.central, "stat_error": entry.stat, "syst_error": entry.syst},
                    new=None,
                    details={"normalized_base_name": norm_base},
                )
            )
            continue

        if enum_name not in observable_flha_mapping:
            missing_flha_mapping.append(
                DiffRecord(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    base_name=entry.base_name,
                    section=entry.section,
                    expected_key=None,
                    reason="missing_flha_mapping",
                    legacy={"central_value": entry.central, "stat_error": entry.stat, "syst_error": entry.syst},
                    new=None,
                    details={"enum_name": enum_name},
                )
            )
            continue

        entry.enum_name = enum_name
        entry.canonical_name = observable_mapping.get(enum_name)
        expected_parts = insert_bin_parts(observable_flha_mapping[enum_name], entry.low, entry.high)
        expected_key = parts_to_key(expected_parts)
        entry.expected_key = expected_key

        sec_entries = new_by_section.get(entry.section)
        if sec_entries is None:
            missing_in_new.append(
                DiffRecord(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    base_name=entry.base_name,
                    section=entry.section,
                    expected_key=expected_key,
                    reason="missing_section_in_new_json",
                    legacy={"central_value": entry.central, "stat_error": entry.stat, "syst_error": entry.syst},
                    new=None,
                    details={"enum_name": enum_name},
                )
            )
            continue

        new_entry = sec_entries.get(expected_key)
        if new_entry is None:
            same_obs_candidates = new_by_unbinned.get(entry.section, {}).get(tuple(observable_flha_mapping[enum_name]), [])
            details = {
                "enum_name": enum_name,
                "canonical_name": entry.canonical_name,
            }
            if same_obs_candidates:
                details["same_observable_other_bins"] = [
                    {
                        "key": cand.key,
                        "bins": list(cand.bins) if cand.bins is not None else None,
                        "central_value": cand.central,
                        "stat_error": cand.stat,
                        "syst_error": cand.syst,
                    }
                    for cand in same_obs_candidates
                ]
                reason = "missing_exact_bin_in_new_json"
            else:
                reason = "missing_in_new_json"

            missing_in_new.append(
                DiffRecord(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    base_name=entry.base_name,
                    section=entry.section,
                    expected_key=expected_key,
                    reason=reason,
                    legacy={"central_value": entry.central, "stat_error": entry.stat, "syst_error": entry.syst},
                    new=None,
                    details=details,
                )
            )
            continue

        matched_keys.add((entry.section, expected_key))
        diffs = compare_entry(entry, new_entry, rel_tol=args.rel_tol, abs_tol=args.abs_tol)
        if diffs:
            value_mismatches.append(
                DiffRecord(
                    line_no=entry.line_no,
                    raw_name=entry.raw_name,
                    base_name=entry.base_name,
                    section=entry.section,
                    expected_key=expected_key,
                    reason="value_mismatch",
                    legacy={"central_value": entry.central, "stat_error": entry.stat, "syst_error": entry.syst},
                    new={"central_value": new_entry.central, "stat_error": new_entry.stat, "syst_error": new_entry.syst},
                    details={"fields": diffs, "enum_name": enum_name, "canonical_name": entry.canonical_name},
                )
            )
        else:
            fully_matched.append(entry)

    extra_in_new = []
    for section, entries in new_by_section.items():
        for key, new_entry in entries.items():
            if (section, key) in matched_keys:
                continue
            extra_in_new.append(
                {
                    "section": section,
                    "key": key,
                    "bins": list(new_entry.bins) if new_entry.bins is not None else None,
                    "central_value": new_entry.central,
                    "stat_error": new_entry.stat,
                    "syst_error": new_entry.syst,
                    "unbinned_parts": new_entry.unbinned_parts,
                }
            )

    summary = {
        "legacy_entries_total": len(legacy_entries),
        "fully_matched": len(fully_matched),
        "value_mismatches": len(value_mismatches),
        "missing_in_new": len(missing_in_new),
        "unknown_legacy_names": len(unknown_legacy),
        "missing_flha_mapping": len(missing_flha_mapping),
        "extra_in_new": len(extra_in_new),
        "alias_collisions": len(alias_collisions),
    }

    print("\n=== Summary ===")
    for k, v in summary.items():
        print(f"{k}: {v}")

    if alias_collisions:
        print("\n=== Alias normalization collisions ===")
        for alias, enums in alias_collisions:
            print(f"- {alias}: {', '.join(enums)}")

    if unknown_legacy:
        print("\n=== Unknown legacy observable names ===")
        for rec in unknown_legacy:
            print(f"- line {rec.line_no}: {rec.raw_name}  (base={rec.base_name}, normalized={rec.details['normalized_base_name']})")

    if missing_flha_mapping:
        print("\n=== Legacy observables resolved but missing FLHA mapping ===")
        for rec in missing_flha_mapping:
            print(f"- line {rec.line_no}: {rec.raw_name} -> {rec.details['enum_name']}")

    if missing_in_new:
        print("\n=== Missing in new JSON / bin mismatches ===")
        for rec in missing_in_new:
            msg = f"- line {rec.line_no}: {rec.raw_name} -> expected {rec.expected_key} in section {rec.section}"
            if rec.reason == "missing_exact_bin_in_new_json":
                msg += "  [same observable found with different bins]"
            elif rec.reason == "missing_section_in_new_json":
                msg += "  [section absent from JSON]"
            print(msg)
            if rec.details and "same_observable_other_bins" in rec.details:
                for cand in rec.details["same_observable_other_bins"]:
                    print(f"    candidate: {cand['key']} bins={cand['bins']}")

    if value_mismatches:
        print("\n=== Numerical mismatches ===")
        for rec in value_mismatches:
            fields = ", ".join(rec.details["fields"])
            print(f"- line {rec.line_no}: {rec.raw_name} [{fields}]")
            print(
                f"    legacy: central={fmt(rec.legacy['central_value'])} stat={fmt(rec.legacy['stat_error'])} syst={fmt(rec.legacy['syst_error'])}"
            )
            print(
                f"    new   : central={fmt(rec.new['central_value'])} stat={fmt(rec.new['stat_error'])} syst={fmt(rec.new['syst_error'])}"
            )
            print(f"    key   : {rec.expected_key} (section={rec.section})")

    if args.show_matched and fully_matched:
        print("\n=== Fully matched entries ===")
        for entry in fully_matched:
            print(f"- line {entry.line_no}: {entry.raw_name} -> {entry.expected_key} (section={entry.section})")

    if extra_in_new:
        print("\n=== Entries present in new JSON but not matched from legacy ===")
        for item in extra_in_new:
            print(f"- {item['section']} :: {item['key']}  bins={item['bins']}")

    if args.json_report is not None:
        report = {
            "summary": summary,
            "alias_collisions": [{"normalized_alias": alias, "enum_names": enums} for alias, enums in alias_collisions],
            "unknown_legacy_names": [asdict(x) for x in unknown_legacy],
            "missing_flha_mapping": [asdict(x) for x in missing_flha_mapping],
            "missing_in_new": [asdict(x) for x in missing_in_new],
            "value_mismatches": [asdict(x) for x in value_mismatches],
            "extra_in_new": extra_in_new,
            "matched": [
                {
                    "line_no": x.line_no,
                    "raw_name": x.raw_name,
                    "base_name": x.base_name,
                    "section": x.section,
                    "low": x.low,
                    "high": x.high,
                    "expected_key": x.expected_key,
                    "enum_name": x.enum_name,
                    "canonical_name": x.canonical_name,
                }
                for x in fully_matched
            ],
        }
        with args.json_report.open("w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        print(f"\nJSON report written to: {args.json_report}")

    has_problems = any(
        [
            alias_collisions,
            unknown_legacy,
            missing_flha_mapping,
            missing_in_new,
            value_mismatches,
        ]
    )
    return 1 if has_problems else 0


if __name__ == "__main__":
    sys.exit(main())
