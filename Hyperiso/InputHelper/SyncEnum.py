#!/usr/bin/env python3
import argparse
import re
import shutil
from pathlib import Path
from typing import List, Optional, Tuple


def strip_cpp_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//.*?$", "", text, flags=re.MULTILINE)
    return text


def extract_enum_members(header_text: str, enum_name: str) -> List[str]:
    pattern = re.compile(
        rf"enum\s+(?:class\s+)?{re.escape(enum_name)}\s*(?:\:\s*[^{{]+)?\s*{{(?P<body>.*?)}}\s*;",
        flags=re.DOTALL,
    )
    match = pattern.search(header_text)
    if not match:
        raise ValueError(f"Enum '{enum_name}' introuvable dans le header.")

    body = strip_cpp_comments(match.group("body"))
    members: List[str] = []

    for raw_item in body.split(","):
        item = raw_item.strip()
        if not item:
            continue

        # Gère NAME, NAME = 3, NAME = expr
        member_match = re.match(r"^([A-Za-z_]\w*)\b", item)
        if member_match:
            members.append(member_match.group(1))

    if not members:
        raise ValueError(f"Aucun membre extrait pour l'enum '{enum_name}'.")

    return members


def find_pybind_block(text: str, enum_name: str) -> Tuple[int, int, str, str]:
    start_pattern = re.compile(
        rf'^(?P<indent>[ \t]*)py::enum_<\s*{re.escape(enum_name)}\s*>\s*\(.*?"{re.escape(enum_name)}".*$',
        flags=re.MULTILINE,
    )
    match = start_pattern.search(text)
    if not match:
        raise ValueError(f"Bloc pybind11 pour '{enum_name}' introuvable.")

    indent = match.group("indent")
    start = match.start()

    # On prend tout jusqu'au premier ';' après le début du bloc
    semicolon_idx = text.find(";", match.end())
    if semicolon_idx == -1:
        raise ValueError(f"Fin de bloc pybind11 pour '{enum_name}' introuvable (pas de ';').")

    end = semicolon_idx + 1
    block = text[start:end]
    return start, end, indent, block


def generate_pybind_block(enum_name: str, members: List[str], indent: str, export_values: bool) -> str:
    lines = [f'{indent}py::enum_<{enum_name}>(m, "{enum_name}")']
    for member in members:
        lines.append(f'{indent}    .value("{member}", {enum_name}::{member})')

    if export_values:
        lines.append(f"{indent}    .export_values();")
    else:
        lines[-1] += ";"

    return "\n".join(lines)


def update_pybind_file(binding_path: Path, enum_name: str, members: List[str]) -> bool:
    text = binding_path.read_text(encoding="utf-8")
    start, end, indent, old_block = find_pybind_block(text, enum_name)
    export_values = ".export_values()" in old_block
    new_block = generate_pybind_block(enum_name, members, indent, export_values)

    if old_block == new_block:
        return False

    updated = text[:start] + new_block + text[end:]
    binding_path.write_text(updated, encoding="utf-8")
    return True


def find_python_class_block(text: str, enum_name: str) -> Tuple[int, int, str, str, str]:
    class_pattern = re.compile(
        rf"^(?P<indent>[ \t]*)class\s+{re.escape(enum_name)}\s*\((?P<bases>[^)]*)\)\s*:\s*\n",
        flags=re.MULTILINE,
    )
    match = class_pattern.search(text)
    if not match:
        raise ValueError(f"Classe Python '{enum_name}' introuvable.")

    indent = match.group("indent")
    bases = match.group("bases").strip()
    start = match.start()
    body_start = match.end()

    # Cherche la prochaine classe de même niveau
    next_class_pattern = re.compile(
        rf"^(?:{re.escape(indent)})class\s+\w+\b.*:\s*$",
        flags=re.MULTILINE,
    )
    next_match = next_class_pattern.search(text, body_start)
    end = next_match.start() if next_match else len(text)

    block = text[start:end].rstrip("\n")
    return start, end, indent, bases, block


def extract_existing_docstring(class_block: str, indent: str) -> Optional[str]:
    parts = class_block.split("\n", 1)
    if len(parts) < 2:
        return None
    body_text = parts[1]

    pattern = re.compile(
        "^" + re.escape(indent) + "    (?P<quote>'''|\"\"\")(?P<doc>.*?)(?P=quote)",
        flags=re.DOTALL,
    )
    match = pattern.search(body_text)
    if not match:
        return None

    quote = match.group("quote")
    doc = match.group("doc")
    return f"{indent}    {quote}{doc}{quote}"


def generate_python_class(
    enum_name: str,
    members: List[str],
    indent: str,
    bases: str,
    docstring: Optional[str],
) -> str:
    cpp_enum_name = f"_Cpp{enum_name}"
    lines = [f"{indent}class {enum_name}({bases}):"]

    if docstring:
        lines.append(docstring)
        lines.append("")

    for member in members:
        lines.append(f"{indent}    {member} = {cpp_enum_name}.{member}")

    return "\n".join(lines)


def update_python_wrapper(python_path: Path, enum_name: str, members: List[str]) -> bool:
    text = python_path.read_text(encoding="utf-8")
    start, end, indent, bases, old_block = find_python_class_block(text, enum_name)
    docstring = extract_existing_docstring(old_block, indent)
    new_block = generate_python_class(enum_name, members, indent, bases, docstring)

    if old_block == new_block:
        return False

    tail = text[end:]
    if tail:
        tail = "\n\n" + tail.lstrip("\n")

    updated = text[:start] + new_block + tail
    python_path.write_text(updated, encoding="utf-8")
    return True


def backup_file(path: Path) -> None:
    backup = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Synchronise un enum C++ avec son binding pybind11 et son wrapper Python."
    )
    parser.add_argument("--enum-name", required=True, help="Nom de l'enum (ex: Observables)")
    parser.add_argument("--header", required=True, help="Chemin vers le header C++ source")
    parser.add_argument("--binding", required=True, help="Chemin vers le fichier pybind11")
    parser.add_argument("--python-wrapper", required=True, help="Chemin vers le wrapper Python")
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="N'écrit pas de fichiers .bak avant modification",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Affiche les membres détectés sans modifier les fichiers",
    )

    args = parser.parse_args()

    header_path = Path(args.header)
    binding_path = Path(args.binding)
    python_path = Path(args.python_wrapper)

    header_text = header_path.read_text(encoding="utf-8")
    members = extract_enum_members(header_text, args.enum_name)

    if args.dry_run:
        print(f"Enum: {args.enum_name}")
        print("Membres détectés:")
        for member in members:
            print(f"  - {member}")
        return

    if not args.no_backup:
        backup_file(binding_path)
        backup_file(python_path)

    binding_changed = update_pybind_file(binding_path, args.enum_name, members)
    python_changed = update_python_wrapper(python_path, args.enum_name, members)

    print(f"[binding] {'modifié' if binding_changed else 'déjà à jour'}: {binding_path}")
    print(f"[python ] {'modifié' if python_changed else 'déjà à jour'}: {python_path}")


if __name__ == "__main__":
    main()