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
class ParsedLegacyName:
    raw_name: str
    base_name: str
    section: str
    has_bin: bool
    low: float
    high: float


@dataclass
class LegacyLine:
    line_no: int
    raw_name: str
    central: float
    stat: float
    syst: float


@dataclass
class ResolvedLine:
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
    enum_name: str
    canonical_name: Optional[str]
    key: str
    parts: List[int]


@dataclass
class SkippedLine:
    line_no: int
    raw_name: str
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
        if enum_name not in observable_flha_mapping:
            continue
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


def parse_legacy_lines(path: Path) -> List[LegacyLine]:
    out: List[LegacyLine] = []

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

            out.append(LegacyLine(line_no=line_no, raw_name=name, central=central, stat=stat, syst=syst))

    return out


def candidate_parses(raw_name: str) -> List[ParsedLegacyName]:
    """
    Build plausible interpretations of the legacy observable name.

    We do NOT assume that the last token is always a section.
    Instead we try several structured interpretations and later keep
    the first one whose base name resolves through the mapper.
    """
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

    # Most specific first: section + bin
    if len(tokens) >= 4 and (not is_number_token(tokens[-1])) and is_number_token(tokens[-2]) and is_number_token(tokens[-3]):
        add(tokens[:-3], tokens[-1], True, float(tokens[-3]), float(tokens[-2]))

    # Bin only
    if len(tokens) >= 3 and is_number_token(tokens[-1]) and is_number_token(tokens[-2]):
        add(tokens[:-2], "DEFAULT", True, float(tokens[-2]), float(tokens[-1]))

    # Section only
    if len(tokens) >= 2 and (not is_number_token(tokens[-1])):
        add(tokens[:-1], tokens[-1], False, 0.0, 0.0)

    # Plain name
    add(tokens, "DEFAULT", False, 0.0, 0.0)

    return candidates


def resolve_legacy_name(
    raw_name: str,
    alias_map: Dict[str, str],
    observable_mapping: Dict[str, str],
    observable_flha_mapping: Dict[str, Sequence[int]],
) -> Tuple[Optional[ParsedLegacyName], Optional[str], Optional[str], Dict[str, Any]]:
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
            return cand, enum_name, observable_mapping.get(enum_name), {
                "attempts": attempts,
                "reason": "missing_flha_mapping",
            }
        return cand, enum_name, observable_mapping.get(enum_name), {"attempts": attempts}

    return None, None, None, {"attempts": attempts, "reason": "unknown_legacy_name"}


def ensure_output_root(merge_with: Optional[Path]) -> Dict[str, Any]:
    if merge_with is None:
        return {"FOBS": {}}

    with merge_with.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError(f"{merge_with}: existing JSON is not an object")
    if "FOBS" not in data:
        data["FOBS"] = {}
    if not isinstance(data["FOBS"], dict):
        raise ValueError(f"{merge_with}: top-level FOBS is not an object")
    return data


def make_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert legacy flat observable input format into the new JSON FOBS format.")
    p.add_argument("legacy_file", type=Path, help="Old flat text input file")
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
        "--fail-on-conflict",
        action="store_true",
        help="Exit with non-zero status if the same output key is produced with different values",
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

    alias_map, alias_collisions = build_alias_map(
        observable_mapping=observable_mapping,
        observable_flha_mapping=observable_flha_mapping,
        manual_aliases=manual_aliases,
    )

    legacy_lines = parse_legacy_lines(args.legacy_file)
    out = ensure_output_root(args.merge_with)
    fobs = out["FOBS"]

    converted: List[ResolvedLine] = []
    skipped: List[SkippedLine] = []
    conflicts: List[Dict[str, Any]] = []
    duplicate_identical: List[Dict[str, Any]] = []

    for line in legacy_lines:
        parsed, enum_name, canonical_name, meta = resolve_legacy_name(
            raw_name=line.raw_name,
            alias_map=alias_map,
            observable_mapping=observable_mapping,
            observable_flha_mapping=observable_flha_mapping,
        )

        if parsed is None or enum_name is None:
            skipped.append(
                SkippedLine(
                    line_no=line.line_no,
                    raw_name=line.raw_name,
                    reason=meta.get("reason", "unknown_legacy_name"),
                    details=meta,
                )
            )
            continue

        if enum_name not in observable_flha_mapping:
            skipped.append(
                SkippedLine(
                    line_no=line.line_no,
                    raw_name=line.raw_name,
                    reason="missing_flha_mapping",
                    details={"enum_name": enum_name, **meta},
                )
            )
            continue

        parts = insert_bin_parts(observable_flha_mapping[enum_name], parsed.low, parsed.high)
        key = parts_to_key(parts)

        section_obj = fobs.setdefault(parsed.section, {})
        new_payload = {
            "central_value": line.central,
            "stat_error": line.stat,
            "syst_error": line.syst,
        }

        previous = section_obj.get(key)
        if previous is not None:
            same = (
                previous.get("central_value") == new_payload["central_value"]
                and previous.get("stat_error") == new_payload["stat_error"]
                and previous.get("syst_error") == new_payload["syst_error"]
            )
            if same:
                duplicate_identical.append(
                    {
                        "line_no": line.line_no,
                        "raw_name": line.raw_name,
                        "section": parsed.section,
                        "key": key,
                    }
                )
            else:
                conflicts.append(
                    {
                        "line_no": line.line_no,
                        "raw_name": line.raw_name,
                        "section": parsed.section,
                        "key": key,
                        "previous": previous,
                        "new": new_payload,
                        "enum_name": enum_name,
                    }
                )

        section_obj[key] = new_payload

        converted.append(
            ResolvedLine(
                line_no=line.line_no,
                raw_name=line.raw_name,
                base_name=parsed.base_name,
                section=parsed.section,
                has_bin=parsed.has_bin,
                low=parsed.low,
                high=parsed.high,
                central=line.central,
                stat=line.stat,
                syst=line.syst,
                enum_name=enum_name,
                canonical_name=canonical_name,
                key=key,
                parts=parts,
            )
        )

    with args.output_json.open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=args.indent, sort_keys=True)
        fh.write("\n")

    summary = {
        "legacy_lines_total": len(legacy_lines),
        "converted_count": len(converted),
        "skipped_count": len(skipped),
        "conflict_count": len(conflicts),
        "duplicate_identical_count": len(duplicate_identical),
        "alias_collision_count": len(alias_collisions),
        "sections_written": sorted(fobs.keys()),
    }

    print(f"[OK] Wrote {args.output_json}")
    print(f"  total legacy lines       : {summary['legacy_lines_total']}")
    print(f"  converted                : {summary['converted_count']}")
    print(f"  skipped                  : {summary['skipped_count']}")
    print(f"  conflicting duplicates   : {summary['conflict_count']}")
    print(f"  identical duplicates     : {summary['duplicate_identical_count']}")
    print(f"  alias collisions         : {summary['alias_collision_count']}")
    print(f"  sections written         : {', '.join(summary['sections_written']) if summary['sections_written'] else '(none)'}")

    if skipped:
        print("\nSkipped entries:")
        for s in skipped[:30]:
            print(f"  line {s.line_no}: {s.raw_name} -> {s.reason}")
        if len(skipped) > 30:
            print(f"  ... and {len(skipped) - 30} more")

    if conflicts:
        print("\nConflicting duplicates:")
        for c in conflicts[:20]:
            print(f"  line {c['line_no']}: {c['raw_name']} -> {c['section']} / {c['key']}")
        if len(conflicts) > 20:
            print(f"  ... and {len(conflicts) - 20} more")

    if alias_collisions:
        print("\nAlias collisions in mapper normalization:")
        for alias, enums in alias_collisions[:20]:
            print(f"  {alias}: {', '.join(enums)}")
        if len(alias_collisions) > 20:
            print(f"  ... and {len(alias_collisions) - 20} more")

    if args.report_json is not None:
        report = {
            "summary": summary,
            "skipped": [asdict(x) for x in skipped],
            "converted": [asdict(x) for x in converted],
            "conflicts": conflicts,
            "duplicate_identical": duplicate_identical,
            "alias_collisions": [{"alias": a, "enum_names": e} for a, e in alias_collisions],
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