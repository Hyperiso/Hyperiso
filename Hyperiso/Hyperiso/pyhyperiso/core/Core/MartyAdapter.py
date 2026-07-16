"""Python accessors for MARTY-related configuration."""

from __future__ import annotations

from enum import Enum
from typing import Optional

from pyhyperiso.phyperiso.pyhyperiso.core import MartyAdapter as _CppMartyAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import MartyPath as _CppMartyPath


class MartyPath(Enum):
    """MARTY path keys exposed by the C++ runtime.

    Attributes:
        MODEL_FILE: Path to the configured MARTY model file.
        TEMPLATE_DIR: Read-only MARTY template directory.
        PARAM_MAPPING_DIR: Read-only directory containing MARTY/Hyperiso mappings.
        SM_MAPPING_FILE: Read-only SM mapping JSON.
        BSM_MAPPING_FILE: Optional user-provided BSM mapping JSON.
        MARTY_TEMP_DIR: Writable MARTY generated-code/cache directory.
    """

    MODEL_FILE = _CppMartyPath.MODEL_FILE
    TEMPLATE_DIR = _CppMartyPath.TEMPLATE_DIR
    PARAM_MAPPING_DIR = _CppMartyPath.PARAM_MAPPING_DIR
    SM_MAPPING_FILE = _CppMartyPath.SM_MAPPING_FILE
    BSM_MAPPING_FILE = _CppMartyPath.BSM_MAPPING_FILE
    MARTY_TEMP_DIR = _CppMartyPath.MARTY_TEMP_DIR


class MartyAdapter:
    """Inspect MARTY paths and flags used by Hyperiso.

    This adapter is read-only from Python. It is typically useful after
    ``HyperisoMaster.init(...)`` when debugging which MARTY model, mappings,
    templates and writable cache directory were selected.
    """

    def __init__(self) -> None:
        """Create an adapter bound to the current C++ MARTY configuration."""
        self._cpp_obj = _CppMartyAdapter()

    def get_path(self, path_kind: MartyPath) -> str:
        """Return a required MARTY path by kind."""
        return str(self._cpp_obj.get_path(path_kind.value))

    def get_optional_path(self, path_kind: MartyPath) -> Optional[str]:
        """Return an optional MARTY path, or ``None`` when it is not configured."""
        value = self._cpp_obj.get_optional_path(path_kind.value)
        return None if value is None else str(value)

    def check_flag(self, flag) -> bool:
        """Return the value of a MARTY-related C++ flag."""
        return bool(self._cpp_obj.check_flag(flag))

    def get_model_name(self) -> str:
        """Return the configured MARTY model name."""
        return str(self._cpp_obj.get_marty_model_name())

    def __repr__(self) -> str:
        """Return a compact representation including the model name."""
        return f"<PyMartyAdapter model='{self.get_model_name()}'>"


__all__ = ["MartyAdapter", "MartyPath"]
