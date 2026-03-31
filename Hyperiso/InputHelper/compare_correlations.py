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
class NewCorrelation:
    section: str
    pair_key: Tuple[str, str]
    key1: str
    key2: str
    stat_correlation: float
    syst_correlation: Optional[float]
    raw_record: Dict[str, Any]


@dataclass
class DiffRecord:
    line_no: Optional[int]
    section: str
    raw_obs1: Optional[str]
    raw_obs2: Optional[str]
    expected_pair: Optional[Tuple[str, str]]
    reason: str
    legacy_stat: Optional[float]
    new_stat: Optional[float]
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


def parse_new_correlations(path: Path) -> Tuple[Dict[str, Dict[Tuple[str, str], NewCorrelation]], List[Dict[str, Any]]]:
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")

    correlations = data.get("correlations")
    if correlations is None:
        raise ValueError(f"{path}: missing top-level key 'correlations'")
    if not isinstance(correlations, dict):
        raise ValueError(f"{path}: 'correlations' must be an object mapping section -> list")

    out: Dict[str, Dict[Tuple[str, str], NewCorrelation]] = {}
    conflicts: List[Dict[str, Any]] = []

    for section, items in correlations.items():
        if not isinstance(items, list):
            raise ValueError(f"{path}: correlations[{section!r}] must be a list")

        sec_dict: Dict[Tuple[str, str], NewCorrelation] = {}

        for idx, item in enumerate(items):
            if not isinstance(item, dict):
                raise ValueError(f"{path}: correlations[{section!r}][{idx}] must be an object")

            try:
                key1 = str(item["id_1"])
                key2 = str(item["id_2"])
                stat = float(item["stat_correlation"])
            except KeyError as exc:
                raise ValueError(f"{path}: correlations[{section!r}][{idx}] missing key {exc}") from exc
            except Exception as exc:
                raise ValueError(f"{path}: correlations[{section!r}][{idx}] has invalid fields") from exc

            syst = item.get("syst_correlation")
            syst = None if syst is None else float(syst)

            pair = canonical_pair(key1, key2)
            rec = NewCorrelation(
                section=section,
                pair_key=pair,
                key1=key1,
                key2=key2,
                stat_correlation=stat,
                syst_correlation=syst,
                raw_record=item,
            )

            prev = sec_dict.get(pair)
            if prev is not None:
                same = (
                    prev.stat_correlation == rec.stat_correlation
                    and prev.syst_correlation == rec.syst_correlation
                )
                if not same:
                    conflicts.append(
                        {
                            "section": section,
                            "pair_key": pair,
                            "previous": prev.raw_record,
                            "new": rec.raw_record,
                        }
                    )
            sec_dict[pair] = rec

        out[section] = sec_dict

    return out, conflicts


def is_close(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    return math.isclose(a, b, abs_tol=abs_tol, rel_tol=rel_tol)


def make_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Compare legacy correlation inputs against the new JSON correlations format.")
    p.add_argument("legacy_file", type=Path, help="Old flat text correlation file")
    p.add_argument("new_json", type=Path, help="New JSON correlation file")
    p.add_argument(
        "--maps",
        type=Path,
        default=Path("observable_maps_fixed.py"),
        help="Python module containing OBSERVABLE_MAPPING / OBSERVABLE_FLHA_MAPPING / optional MANUAL_LEGACY_ALIASES",
    )
    p.add_argument("--abs-tol", type=float, default=1e-12, help="Absolute tolerance for correlation comparison")
    p.add_argument("--rel-tol", type=float, default=1e-12, help="Relative tolerance for correlation comparison")
    p.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional JSON report file",
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
    new_corrs, new_conflicts = parse_new_correlations(args.new_json)

    diffs: List[DiffRecord] = []
    matched_new_pairs: set[Tuple[str, Tuple[str, str]]] = set()
    matched_old_count = 0

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
            diffs.append(
                DiffRecord(
                    line_no=entry.line_no,
                    section="?",
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    expected_pair=None,
                    reason="unknown_obs1",
                    legacy_stat=entry.stat_correlation,
                    new_stat=None,
                    details=meta1,
                )
            )
            continue

        if obs2 is None:
            diffs.append(
                DiffRecord(
                    line_no=entry.line_no,
                    section=obs1.section,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    expected_pair=None,
                    reason="unknown_obs2",
                    legacy_stat=entry.stat_correlation,
                    new_stat=None,
                    details=meta2,
                )
            )
            continue

        if obs1.section != obs2.section:
            diffs.append(
                DiffRecord(
                    line_no=entry.line_no,
                    section=f"{obs1.section} vs {obs2.section}",
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    expected_pair=None,
                    reason="inconsistent_legacy_sections",
                    legacy_stat=entry.stat_correlation,
                    new_stat=None,
                    details={
                        "obs1": asdict(obs1),
                        "obs2": asdict(obs2),
                    },
                )
            )
            continue

        section = obs1.section
        pair = canonical_pair(obs1.key, obs2.key)
        new_rec = new_corrs.get(section, {}).get(pair)

        if new_rec is None:
            diffs.append(
                DiffRecord(
                    line_no=entry.line_no,
                    section=section,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    expected_pair=pair,
                    reason="missing_in_new_json",
                    legacy_stat=entry.stat_correlation,
                    new_stat=None,
                    details={
                        "obs1": asdict(obs1),
                        "obs2": asdict(obs2),
                    },
                )
            )
            continue

        matched_new_pairs.add((section, pair))
        matched_old_count += 1

        if not is_close(entry.stat_correlation, new_rec.stat_correlation, args.abs_tol, args.rel_tol):
            diffs.append(
                DiffRecord(
                    line_no=entry.line_no,
                    section=section,
                    raw_obs1=entry.raw_obs1,
                    raw_obs2=entry.raw_obs2,
                    expected_pair=pair,
                    reason="stat_correlation_mismatch",
                    legacy_stat=entry.stat_correlation,
                    new_stat=new_rec.stat_correlation,
                    details={
                        "obs1": asdict(obs1),
                        "obs2": asdict(obs2),
                        "new_record": new_rec.raw_record,
                    },
                )
            )

    unmatched_new: List[Dict[str, Any]] = []
    for section, sec_dict in new_corrs.items():
        for pair, rec in sec_dict.items():
            if (section, pair) not in matched_new_pairs:
                unmatched_new.append(
                    {
                        "section": section,
                        "pair_key": pair,
                        "record": rec.raw_record,
                    }
                )

    summary = {
        "legacy_pairs_total": len(legacy_corrs),
        "matched_old_pairs": matched_old_count,
        "diff_count": len(diffs),
        "new_unmatched_count": len(unmatched_new),
        "new_conflict_count": len(new_conflicts),
        "sections_in_new_json": sorted(new_corrs.keys()),
    }

    print(f"[OK] Compared legacy correlations against {args.new_json}")
    print(f"  total legacy pairs      : {summary['legacy_pairs_total']}")
    print(f"  matched old pairs       : {summary['matched_old_pairs']}")
    print(f"  diffs / missing / errs  : {summary['diff_count']}")
    print(f"  unmatched new pairs     : {summary['new_unmatched_count']}")
    print(f"  conflicts inside new    : {summary['new_conflict_count']}")
    print(f"  sections in new JSON    : {', '.join(summary['sections_in_new_json']) if summary['sections_in_new_json'] else '(none)'}")

    if diffs:
        print("\nDifferences:")
        for d in diffs[:40]:
            print(
                f"  line {d.line_no if d.line_no is not None else '?'} | {d.reason} | "
                f"{d.raw_obs1 or '?'} <-> {d.raw_obs2 or '?'} | section={d.section}"
            )
        if len(diffs) > 40:
            print(f"  ... and {len(diffs) - 40} more")

    if unmatched_new:
        print("\nPairs present only in new JSON:")
        for x in unmatched_new[:20]:
            print(f"  section={x['section']} pair={x['pair_key']}")
        if len(unmatched_new) > 20:
            print(f"  ... and {len(unmatched_new) - 20} more")

    if new_conflicts:
        print("\nConflicting duplicate pairs inside new JSON:")
        for x in new_conflicts[:20]:
            print(f"  section={x['section']} pair={x['pair_key']}")
        if len(new_conflicts) > 20:
            print(f"  ... and {len(new_conflicts) - 20} more")

    if args.report_json is not None:
        report = {
            "summary": summary,
            "diffs": [asdict(x) for x in diffs],
            "unmatched_new": unmatched_new,
            "new_conflicts": new_conflicts,
        }
        with args.report_json.open("w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"[OK] Wrote report {args.report_json}")

    return 1 if diffs or unmatched_new or new_conflicts else 0


if __name__ == "__main__":
    sys.exit(main())
