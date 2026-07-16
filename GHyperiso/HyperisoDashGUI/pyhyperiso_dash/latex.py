from __future__ import annotations

from typing import Any

from dash import dcc, html

from pyhyperiso_dash.domain import code_to_display
from pyhyperiso_dash.latex_data.decay_name_to_latex_map import DECAY_ENUM_TO_LATEX_MAP
from pyhyperiso_dash.latex_data.new_to_latex_nuisance_map import (
    NEW_TO_LATEX_NUISANCE_MAP,
)
from pyhyperiso_dash.latex_data.observable_flha_to_latex_map import (
    OBSERVABLE_ENUM_TO_LATEX_MAP,
)
from pyhyperiso_dash.latex_data.wilson_lha_to_latex_map import (
    WILSON_LHA_ID_TO_LATEX_MAP,
    WILSON_LHA_ID_TO_NAME_MAP,
    WILSON_NAME_TO_LATEX_MAP,
)


def _strip_math(label: str | None) -> str | None:
    if label is None:
        return None
    text = str(label).strip()
    if len(text) >= 2 and text[0] == "$" and text[-1] == "$":
        return text[1:-1]
    return text


def math(label: str | None, fallback: str = "") -> str:
    """Return a dollar-delimited LaTeX label."""
    body = _strip_math(label) or fallback
    if not body:
        return ""
    return f"${body}$"


def text_latex(label: str | None, fallback: str = "") -> str:
    """Return a readable LaTeX string without dollar delimiters.

    Dash DataTable and some dropdown renderers do not run MathJax inside cell
    content.  Returning the raw LaTeX body avoids showing confusing ``$...$``
    delimiters while still presenting the physics notation.
    """
    return _strip_math(label) or str(fallback or "")


_SUBSCRIPT_TRANSLATION = str.maketrans(
    {
        "0": "₀",
        "1": "₁",
        "2": "₂",
        "3": "₃",
        "4": "₄",
        "5": "₅",
        "6": "₆",
        "7": "₇",
        "8": "₈",
        "9": "₉",
        "+": "₊",
        "-": "₋",
        "=": "₌",
        "(": "₍",
        ")": "₎",
        "a": "ₐ",
        "e": "ₑ",
        "h": "ₕ",
        "i": "ᵢ",
        "j": "ⱼ",
        "k": "ₖ",
        "l": "ₗ",
        "m": "ₘ",
        "n": "ₙ",
        "o": "ₒ",
        "p": "ₚ",
        "r": "ᵣ",
        "s": "ₛ",
        "t": "ₜ",
        "u": "ᵤ",
        "v": "ᵥ",
        "x": "ₓ",
        "S": "ₛ",
        "P": "ₚ",
        "L": "ₗ",
        "T": "ₜ",
        "V": "ᵥ",
    }
)

_SUPERSCRIPT_TRANSLATION = str.maketrans(
    {
        "0": "⁰",
        "1": "¹",
        "2": "²",
        "3": "³",
        "4": "⁴",
        "5": "⁵",
        "6": "⁶",
        "7": "⁷",
        "8": "⁸",
        "9": "⁹",
        "+": "⁺",
        "-": "⁻",
        "=": "⁼",
        "(": "⁽",
        ")": "⁾",
        "a": "ᵃ",
        "b": "ᵇ",
        "c": "ᶜ",
        "d": "ᵈ",
        "e": "ᵉ",
        "f": "ᶠ",
        "g": "ᵍ",
        "h": "ʰ",
        "i": "ⁱ",
        "j": "ʲ",
        "k": "ᵏ",
        "l": "ˡ",
        "m": "ᵐ",
        "n": "ⁿ",
        "o": "ᵒ",
        "p": "ᵖ",
        "r": "ʳ",
        "s": "ˢ",
        "t": "ᵗ",
        "u": "ᵘ",
        "v": "ᵛ",
        "w": "ʷ",
        "x": "ˣ",
        "y": "ʸ",
        "z": "ᶻ",
        "B": "ᴮ",
        "D": "ᴰ",
        "K": "ᴷ",
        "P": "ᴾ",
        "S": "ˢ",
        "T": "ᵀ",
        "V": "ⱽ",
    }
)


def compact_math_text(label: str | None, fallback: str = "") -> str:
    """Return a compact Unicode rendering of the simple Wilson LaTeX map.

    This deliberately returns a plain string rather than a Dash component.
    Large multi-select Dropdowns use react-virtualized-select; component labels
    can make its legacy memoizer hit ``undefined.join`` on initial rendering.
    The Wilson map only needs subscripts, superscripts, primes and a tilde, all
    of which have a stable text representation here.
    """
    text = _strip_math(label) or str(fallback or "")
    text = text.replace(r"\tilde{C}", "C̃").replace(r"\prime", "′")

    import re

    def sub(match):
        return match.group(1).replace("_", "").translate(_SUBSCRIPT_TRANSLATION)

    def sup(match):
        return match.group(1).replace("_", "").translate(_SUPERSCRIPT_TRANSLATION)

    # The supplied Wilson map uses flat {...} groups, so these substitutions
    # are sufficient and avoid embedding React/MathJax nodes in the selector.
    text = re.sub(r"_\{([^{}]*)\}", sub, text)
    text = re.sub(r"\^\{([^{}]*)\}", sup, text)
    return text.replace("{", "").replace("}", "")


def dropdown_label(latex: str | None, raw: str, *, include_raw: bool = True):
    """Return a dropdown label with rendered MathJax and a raw code.

    Dash allows component labels in dropdown options.  For actual MathJax
    rendering we use ``dcc.Markdown(mathjax=True)`` instead of a plain text
    span.  The important stability points are: never include ``None`` in the
    children list, and always provide a plain string ``search`` field in the
    option dictionaries that use this component label.
    """
    raw = str(raw)
    if not latex:
        return raw

    children = [
        dcc.Markdown(
            math(latex, raw),
            mathjax=True,
            className="latex-inline-math",
        )
    ]
    if include_raw:
        children.append(html.Span(raw, className="option-raw"))

    return html.Span(
        className="latex-option",
        title=raw,
        children=children,
    )


def plain_label(latex: str | None, raw: str | None = None) -> str:
    """Return a Plotly/Markdown label, keeping MathJax delimiters."""
    return math(latex, str(raw or "")) if latex else str(raw or "")


def table_label(latex: str | None, raw: str | None = None) -> str:
    """Return a MathJax-friendly label for tables and rendered cells.

    Dash DataTable itself does not run MathJax automatically, but the app ships
    a small MutationObserver asset that asks MathJax to typeset any ``$...$``
    text after Dash updates the DOM.  Keeping dollar delimiters here gives the
    desired visual rendering while companion raw columns preserve traceability.
    """
    return math(latex, str(raw or "")) if latex else str(raw or "")


def observable_latex(obs_name: str | None) -> str | None:
    if not obs_name:
        return None
    return OBSERVABLE_ENUM_TO_LATEX_MAP.get(str(obs_name))


def observable_label(obs_name: str | None) -> str:
    return plain_label(observable_latex(obs_name), str(obs_name or ""))


def observable_table_label(obs_name: str | None) -> str:
    return table_label(observable_latex(obs_name), str(obs_name or ""))


def observable_option(obs_name: str) -> dict:
    latex = observable_latex(obs_name)
    return {
        "label": dropdown_label(latex, obs_name),
        "value": obs_name,
        "title": obs_name,
        "search": f"{obs_name} {text_latex(latex, '')}",
    }


def decay_latex(decay_name: str | None) -> str | None:
    if not decay_name:
        return None
    return DECAY_ENUM_TO_LATEX_MAP.get(str(decay_name))


def decay_label(decay_name: str | None) -> str:
    return plain_label(decay_latex(decay_name), str(decay_name or ""))


def decay_table_label(decay_name: str | None) -> str:
    return table_label(decay_latex(decay_name), str(decay_name or ""))


def decay_option(decay_name: str) -> dict:
    latex = decay_latex(decay_name)
    return {
        "label": dropdown_label(latex, decay_name),
        "value": decay_name,
        "title": decay_name,
        "search": f"{decay_name} {text_latex(latex, '')}",
    }


def _int_or_none(value: Any) -> int | None:
    try:
        text = str(value).strip()
        if text.startswith("+"):
            text = text[1:]
        if text and text.lstrip("-").isdigit():
            return int(text)
    except Exception:
        pass
    return None


def _code_parts(code: Any) -> list[int]:
    text = code_to_display(code).replace(",", "_").replace(" ", "_")
    out: list[int] = []
    for part in text.split("_"):
        part = part.strip()
        if not part:
            continue
        val = _int_or_none(part)
        if val is not None:
            out.append(val)
    return out


def wilson_parameter_latex(block: str | None, code: Any) -> str | None:
    """Return a Wilson coefficient label from either raw name or LHA ids.

    Wilson parameters can appear in several forms depending on the block used by
    HyperISO: sometimes the code is already the full pair ``block_id_coeff_id``;
    sometimes the block carries the numeric Wilson block id and the code carries
    only the coefficient id.  This helper tries both forms before falling back.
    """
    block_s = str(block or "")
    code_s = code_to_display(code)

    # Direct raw-name lookup, useful for WCoeff enum names.
    for key in (code_s, block_s):
        if key in WILSON_NAME_TO_LATEX_MAP:
            return WILSON_NAME_TO_LATEX_MAP[key]

    block_i = _int_or_none(block_s)
    parts = _code_parts(code_s)

    candidates: list[tuple[int, int]] = []
    if block_i is not None:
        for p in parts:
            candidates.append((block_i, p))
    if len(parts) >= 2:
        candidates.append((parts[0], parts[1]))
        candidates.append((parts[-2], parts[-1]))
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                candidates.append((parts[i], parts[j]))

    seen: set[tuple[int, int]] = set()
    for key in candidates:
        if key in seen:
            continue
        seen.add(key)
        label = WILSON_LHA_ID_TO_LATEX_MAP.get(key)
        if label:
            return label
        raw_name = WILSON_LHA_ID_TO_NAME_MAP.get(key)
        if raw_name:
            return WILSON_NAME_TO_LATEX_MAP.get(raw_name)

    return None


def parameter_latex(
    block: str | None, code: Any, param_type: str | None = None
) -> str | None:
    if not block or code in (None, ""):
        return None
    block_s = str(block)
    code_s = code_to_display(code)

    # For ParameterType.WILSON, use the Wilson-specific map first.
    if str(param_type or "").upper() == "WILSON":
        return wilson_parameter_latex(block_s, code_s) or NEW_TO_LATEX_NUISANCE_MAP.get(
            (block_s, code_s)
        )

    # Nuisance/ordinary parameter map first, with Wilson as a permissive fallback.
    return NEW_TO_LATEX_NUISANCE_MAP.get((block_s, code_s)) or wilson_parameter_latex(
        block_s, code_s
    )


def parameter_label(
    block: str | None, code: Any, raw: str | None = None, param_type: str | None = None
) -> str:
    code_s = code_to_display(code)
    return plain_label(
        parameter_latex(block, code_s, param_type), raw or f"{block}:{code_s}"
    )


def parameter_table_label(
    block: str | None, code: Any, raw: str | None = None, param_type: str | None = None
) -> str:
    code_s = code_to_display(code)
    return table_label(
        parameter_latex(block, code_s, param_type), raw or f"{block}:{code_s}"
    )


_LATEX_SYMBOLS = {
    r"\alpha": "α",
    r"\beta": "β",
    r"\gamma": "γ",
    r"\delta": "δ",
    r"\Delta": "Δ",
    r"\epsilon": "ε",
    r"\lambda": "λ",
    r"\mu": "μ",
    r"\nu": "ν",
    r"\rho": "ρ",
    r"\sigma": "σ",
    r"\tau": "τ",
    r"\phi": "φ",
    r"\Phi": "Φ",
    r"\pi": "π",
    r"\bar": "",
    r"\parallel": "∥",
    r"\perp": "⊥",
    r"\to": "→",
    r"\rightarrow": "→",
    r"\ell": "ℓ",
    r"\pm": "±",
    r"\times": "×",
    r"\cdot": "·",
    r"\mathrm": "",
    r"\text": "",
    r"\operatorname": "",
    r"\left": "",
    r"\right": "",
    r"\,": " ",
    r"\;": " ",
    r"\!": "",
}


def _latex_token_to_text(token: str) -> str:
    """Return a compact text fallback for a LaTeX token.

    This is intentionally conservative and used only inside Dash dropdown
    options.  React-Select can crash when MathJax mutates clickable option DOM
    nodes, so parameter-code options use static HTML ``sub``/``sup`` nodes
    instead of ``dcc.Markdown(mathjax=True)``.
    """
    text = _strip_math(token) or ""
    for src, dst in _LATEX_SYMBOLS.items():
        text = text.replace(src, dst)
    text = text.replace("{", "").replace("}", "")
    text = text.replace("\\", "")
    return " ".join(text.split())


def _read_braced_group(text: str, start: int) -> tuple[str, int]:
    """Read a LaTeX group starting at ``start`` and return body/new index."""
    if start >= len(text):
        return "", start
    if text[start] != "{":
        return text[start], start + 1
    depth = 0
    out = []
    i = start
    while i < len(text):
        ch = text[i]
        if ch == "{":
            depth += 1
            if depth > 1:
                out.append(ch)
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return "".join(out), i + 1
            out.append(ch)
        else:
            out.append(ch)
        i += 1
    return "".join(out), i


def _latex_static_nodes(label: str) -> list[Any]:
    """Convert a useful subset of LaTeX to static HTML nodes.

    The goal is not to be a full TeX renderer; it gives stable visual labels in
    dropdowns without invoking MathJax inside React-Select menus.  It supports
    the nuisance/parameter labels used most often, especially forms such as
    ``f_{B_s}``, ``m_b`` and ``C_{9}^{\\prime}``.
    """
    text = _strip_math(label) or ""
    nodes: list[Any] = []
    buf: list[str] = []

    def flush():
        if buf:
            nodes.append(_latex_token_to_text("".join(buf)))
            buf.clear()

    i = 0
    while i < len(text):
        ch = text[i]
        if ch in "_^":
            flush()
            body, j = _read_braced_group(text, i + 1)
            content = _latex_token_to_text(body)
            if content:
                cls = "latex-html-sub" if ch == "_" else "latex-html-sup"
                tag = html.Sub if ch == "_" else html.Sup
                nodes.append(tag(content, className=cls))
            i = j
            continue
        if ch == "\\":
            matched = False
            for cmd in sorted(_LATEX_SYMBOLS, key=len, reverse=True):
                if text.startswith(cmd, i):
                    buf.append(_LATEX_SYMBOLS[cmd])
                    i += len(cmd)
                    matched = True
                    break
            if matched:
                continue
            # Unknown command: keep its name without the backslash.
            j = i + 1
            while j < len(text) and text[j].isalpha():
                j += 1
            buf.append(text[i + 1 : j])
            i = j
            continue
        if ch in "{}":
            i += 1
            continue
        buf.append(ch)
        i += 1
    flush()
    return [n for n in nodes if not (isinstance(n, str) and n == "")]


def parameter_dropdown_label(latex: str | None, raw_code: str, raw_full: str):
    """Return a stable parameter-code dropdown label.

    If a LaTeX name is known, it is shown first using static HTML nodes
    (``sub``/``sup``) rather than MathJax.  This avoids the intermittent
    React-Select crash ``Cannot read properties of undefined (reading 'props')``
    observed when MathJax/Markdown components are embedded inside clickable
    dropdown options.  The raw code remains visible on the right.
    """
    raw_code = str(raw_code)
    raw_full = str(raw_full)
    children = []
    if latex:
        rendered = _latex_static_nodes(latex)
        if rendered:
            children.append(
                html.Span(
                    rendered,
                    className="parameter-code-latex parameter-code-latex-static",
                )
            )
    children.append(html.Span(raw_code, className="parameter-code-raw"))
    return html.Span(
        className="latex-option parameter-code-option",
        title=raw_full,
        children=children,
    )


def parameter_option(block: str, code: Any, param_type: str | None = None) -> dict:
    code_s = code_to_display(code)
    raw = f"{block}:{code_s}"
    latex = parameter_latex(block, code_s, param_type)
    return {
        "label": parameter_dropdown_label(latex, code_s, raw),
        "value": code_s,
        "title": raw,
        "search": f"{raw} {code_s} {text_latex(latex, '')}",
    }


def wilson_latex(coeff_name: str | None) -> str | None:
    if not coeff_name:
        return None
    return WILSON_NAME_TO_LATEX_MAP.get(str(coeff_name))


def wilson_base_label(coeff_name: str | None) -> str:
    return plain_label(wilson_latex(coeff_name), str(coeff_name or ""))


def wilson_table_label(coeff_name: str | None) -> str:
    return table_label(wilson_latex(coeff_name), str(coeff_name or ""))


def wilson_option(coeff_name: str) -> dict:
    latex = wilson_latex(coeff_name)
    return {
        # Static sub/sup nodes keep the physics notation readable without
        # letting MathJax mutate React-Select's clickable option tree.
        "label": parameter_dropdown_label(latex, coeff_name, coeff_name),
        "value": coeff_name,
        "title": coeff_name,
        "search": f"{coeff_name} {text_latex(latex, '')}",
    }


def _order_index(order_name: str | None) -> str:
    return {"LO": "0", "NLO": "1", "NNLO": "2", "NONE": "0"}.get(
        str(order_name or "LO"), "0"
    )


def wilson_request_latex(
    method: str, coeff_name: str, order_name: str, contribution: str
) -> str:
    """Return the displayed Wilson coefficient with order/full/contribution marks."""
    base = _strip_math(wilson_latex(coeff_name)) or str(coeff_name)
    n = _order_index(order_name)
    full = str(method).upper() in {"FM", "FR"}
    cont = str(contribution or "TOTAL").upper()

    bracket = r"\left[" if full else r"\left("
    close = r"\right]" if full else r"\right)"
    order_part = rf"^{{\leq {n}}}" if full else rf"^{{({n})}}"

    if cont == "BSM":
        body = rf"\delta {bracket}{base}{close}{order_part}"
    elif cont == "SM":
        body = rf"{bracket}{base}{close}{order_part}_{{\rm SM}}"
    else:
        body = rf"{bracket}{base}{close}{order_part}_{{\rm TOT}}"

    return f"${body}$"


def strip_for_title(label: str) -> str:
    return label or ""


def wilson_request_table_latex(
    method: str, coeff_name: str, order_name: str, contribution: str
) -> str:
    """Return a MathJax-friendly Wilson request label for tables."""
    return wilson_request_latex(method, coeff_name, order_name, contribution)
