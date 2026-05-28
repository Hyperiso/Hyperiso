from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Iterable, Sequence


def enum_by_name(enum_cls: type, name: str):
    if name is None:
        raise ValueError(f"Missing value for {enum_cls.__name__}")
    try:
        return enum_cls[str(name)]
    except KeyError as exc:
        valid = ", ".join(x.name for x in enum_cls)
        raise ValueError(f"Unknown {enum_cls.__name__}: {name}. Valid values: {valid}") from exc


def as_float(value: Any, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except Exception:
        return default


def as_int(value: Any, default: int | None = None) -> int | None:
    if value is None or value == "":
        return default
    try:
        return int(value)
    except Exception:
        return default


def _call_or_value(obj: Any, name: str) -> Any:
    """Return ``obj.name`` or ``obj.name()`` when it is callable."""
    attr = getattr(obj, name)
    return attr() if callable(attr) else attr


def _primitive_float(value: Any) -> float | None:
    """Convert obvious numeric/scalar objects without recursive unwrapping."""
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, complex):
        if abs(value.imag) > 0.0:
            raise TypeError("complex value is not purely real")
        return float(value.real)

    # pyhyperiso scalar_t is not always registered as a Python float. The
    # wrapper exposes to_double(), real() and imag(); prefer those over float().
    for name in ("to_double", "real"):
        if hasattr(value, name):
            try:
                attr = getattr(value, name)
                out = attr() if callable(attr) else attr
                if out is value:
                    continue
                if isinstance(out, (int, float)):
                    return float(out)
            except Exception:
                pass

    try:
        return float(value)
    except Exception:
        return None


def scalar_components(value: Any, _seen: set[int] | None = None) -> tuple[float, float]:
    """Return real and imaginary parts of a scalar-like C++ value.

    Hyperiso uses ``scalar_t`` for Wilson coefficients and several parameter or
    observable quantities. Depending on the pybind11 binding, that object may be
    exposed as a native number, a Python complex, a wrapper object with
    ``value``/``_cpp_obj``, or a bound C++ object with real/imag accessors. This
    helper keeps the Dash layer independent from the exact binding shape.

    Args:
        value: Python number, complex, wrapper, or bound C++ ``scalar_t``.
        _seen: Internal recursion guard.

    Returns:
        ``(real, imag)`` as Python floats.

    Raises:
        TypeError: If the value cannot be interpreted as a scalar.
    """
    if _seen is None:
        _seen = set()
    oid = id(value)
    if oid in _seen:
        raise TypeError(f"Recursive scalar wrapper while converting {value!r}")
    _seen.add(oid)

    if isinstance(value, bool):
        return float(value), 0.0
    if isinstance(value, (int, float)):
        return float(value), 0.0
    if isinstance(value, complex):
        return float(value.real), float(value.imag)

    # Fast path for pyhyperiso's Python Scalar wrapper and raw scalar_t binding.
    # User-side Scalar.py exposes real(), imag() and to_double(); raw pybind
    # scalar_t usually exposes the same methods.
    if hasattr(value, "real") or hasattr(value, "to_double"):
        try:
            real_obj = _call_or_value(value, "real") if hasattr(value, "real") else _call_or_value(value, "to_double")
            imag_obj = _call_or_value(value, "imag") if hasattr(value, "imag") else 0.0
            real = _primitive_float(real_obj)
            imag = _primitive_float(imag_obj)
            if real is not None and imag is not None:
                return real, imag
        except Exception:
            pass

    try:
        c = complex(value)
        return float(c.real), float(c.imag)
    except Exception:
        pass

    for real_name, imag_name in [
        ("real", "imag"),
        ("real", "imaginary"),
        ("re", "im"),
        ("Re", "Im"),
        ("get_real", "get_imag"),
        ("real_part", "imag_part"),
    ]:
        if hasattr(value, real_name):
            try:
                real_obj = _call_or_value(value, real_name)
                imag_obj = _call_or_value(value, imag_name) if hasattr(value, imag_name) else 0.0
                real = _primitive_float(real_obj)
                imag = _primitive_float(imag_obj)
                if real is not None and imag is not None:
                    return real, imag
                return (
                    scalar_components(real_obj, _seen)[0],
                    scalar_components(imag_obj, _seen)[0] if imag_obj is not None else 0.0,
                )
            except Exception:
                pass

    for attr_name in ("value", "val", "scalar", "_cpp_obj"):
        if hasattr(value, attr_name):
            try:
                inner = _call_or_value(value, attr_name)
                if inner is not value:
                    return scalar_components(inner, _seen)
            except Exception:
                pass

    text = str(value)
    # Fallbacks for common pybind/repr shapes, e.g. ``(1.2,-0.3)``,
    # ``1.2+0.3i`` or ``scalar_t(real=1.2, imag=-0.3)``. Do not parse
    # default Python object reprs because their memory addresses contain digits.
    if "object at 0x" in text or (text.startswith("<") and text.endswith(">")):
        raise TypeError(f"Cannot convert opaque scalar object {type(value)!r}; expose real/imag accessors in the binding.")
    nums = re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    if nums:
        if len(nums) >= 2 and any(token in text.lower() for token in ("j", "i", "imag", ",")):
            return float(nums[0]), float(nums[1])
        return float(nums[0]), 0.0

    raise TypeError(f"Cannot convert {value!r} to scalar components")


def scalar_to_float(value: Any, component: str = "real") -> float:
    """Convert a scalar-like value to a float component.

    Args:
        value: Python scalar, complex, wrapper, or bound C++ ``scalar_t``.
        component: Component to extract. Supported values are ``"real"``,
            ``"imag"``/``"imaginary"``, ``"abs"``/``"magnitude"`` and
            ``"phase"``.

    Returns:
        Requested component as a Python float.
    """
    real, imag = scalar_components(value)
    key = str(component or "real").lower()
    if key in {"real", "re"}:
        return real
    if key in {"imag", "imaginary", "im"}:
        return imag
    if key in {"abs", "absolute", "magnitude", "mod", "modulus"}:
        return math.hypot(real, imag)
    if key in {"phase", "arg", "angle"}:
        return math.atan2(imag, real)
    raise ValueError(f"Unknown scalar component: {component}")


def scalar_to_display(value: Any) -> str:
    """Return a compact display string for a scalar-like value."""
    real, imag = scalar_components(value)
    if abs(imag) < 1e-15:
        return f"{real:.12g}"
    sign = "+" if imag >= 0 else "-"
    return f"{real:.12g} {sign} {abs(imag):.12g}i"


def code_to_display(code: Any) -> str:
    if hasattr(code, "to_string"):
        try:
            return str(code.to_string())
        except Exception:
            pass
    if hasattr(code, "get_parts"):
        try:
            return "_".join(str(x) for x in code.get_parts())
        except Exception:
            pass
    return str(code)


def key_to_code(key: Any) -> Any:
    if hasattr(key, "to_string") or hasattr(key, "get_parts"):
        return key
    return key


def parse_code(value: Any) -> Any:
    if value is None or value == "":
        return 0
    if isinstance(value, (int, float)) and float(value).is_integer():
        return int(value)
    text = str(value).strip()
    if "," in text:
        return [int(x.strip()) for x in text.split(",") if x.strip()]
    if "_" in text:
        return text
    try:
        return int(text)
    except Exception:
        return text


def linspace(min_value: Any, max_value: Any, n_points: Any) -> list[float]:
    lo = float(min_value)
    hi = float(max_value)
    n = max(2, int(n_points))
    if n == 1:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def bins_from_step(low: Any, high: Any, step: Any) -> list[tuple[float, float]]:
    lo = float(low)
    hi = float(high)
    st = float(step)
    if st <= 0:
        raise ValueError("Bin step must be strictly positive")
    bins = []
    cur = lo
    guard = 0
    while cur < hi - 1e-15 and guard < 10000:
        nxt = min(hi, cur + st)
        bins.append((float(cur), float(nxt)))
        cur = nxt
        guard += 1
    return bins


def maybe_path(value: str | None) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    return Path(str(value).strip()).expanduser()


def split_csv(text: str | None) -> list[str]:
    if not text:
        return []
    return [x.strip() for x in str(text).replace("\n", ",").split(",") if x.strip()]
