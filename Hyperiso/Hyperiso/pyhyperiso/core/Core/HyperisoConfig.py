"""Configuration objects used to initialize Hyperiso."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

from pyhyperiso.phyperiso.pyhyperiso.core import ExternalFlag as _CppExternalFlag
from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoConfig as _CppHyperisoConfig
from pyhyperiso.core.Common.GeneralEnum import MartyOrderPolicy, Model

PathLike = Union[str, Path]


class ExternalFlag(Enum):
    """External input flags understood by the C++ initialization layer.

    Attributes:
        IS_LHA_SPECTRUM: Whether the input file is a spectrum-style LHA file.
        HAS_WILSON_INPUT: Whether Wilson coefficients are provided externally.
        HAS_TH_OBSERVABLE_INPUT: Whether theory-observable inputs are provided
            externally.
        HYP_AS_SM_MARTY: Whether Hyperiso should expose SM-like inputs to the
            MARTY backend.
    """

    IS_LHA_SPECTRUM = _CppExternalFlag.IS_LHA_SPECTRUM
    HAS_WILSON_INPUT = _CppExternalFlag.HAS_WILSON_INPUT
    HAS_TH_OBSERVABLE_INPUT = _CppExternalFlag.HAS_TH_OBSERVABLE_INPUT
    HYP_AS_SM_MARTY = _CppExternalFlag.HYP_AS_SM_MARTY


@dataclass
class HyperisoConfig:
    """Python configuration for ``HyperisoMaster``.

    Attributes:
        flags: Mapping of external flags to booleans. Missing flags are not
            explicitly set in C++; the default mapping mirrors the common
            Python workflow.
        model: Physics model selected for the session.
        mty_model_name: Optional MARTY model name used by the C++ backend.
        mty_model_path: Optional path to the MARTY model directory or file.
        mty_bsm_mapping_path: Optional user-provided BSM MARTY/Hyperiso mapping JSON.
        mty_order_policy: Explicit MARTY order policy for BSM Wilson coefficients.
        mty_tree_fermion_order: Optional common TreeLevel external-fermion order.

    Examples:
        >>> from pathlib import Path
        >>> cfg = HyperisoConfig(
        ...     model=Model.SM,
        ...     mty_model_name="MSSM_UFO",
        ...     mty_model_path=Path("/path/to/marty/model"),
        ...     mty_bsm_mapping_path=Path("/path/to/zprime_mapping.json"),
        ... )
        >>> cpp_cfg = cfg.to_cpp()
    """

    flags: Dict[ExternalFlag, bool] = field(
        default_factory=lambda: {
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: False,
        }
    )
    model: Model = Model.SM
    mty_model_name: Optional[str] = None
    mty_model_path: Optional[PathLike] = None
    mty_bsm_mapping_path: Optional[PathLike] = None
    mty_order_policy: MartyOrderPolicy = MartyOrderPolicy.AUTO
    mty_tree_fermion_order: Optional[Sequence[int]] = None

    def to_cpp(self) -> _CppHyperisoConfig:
        """Convert this Python config into the bound C++ config.

        Returns:
            C++ ``HyperisoConfig`` instance ready to pass to pybind methods.
        """
        cpp = _CppHyperisoConfig()
        cpp.flags = {flag.value: val for flag, val in self.flags.items()}
        cpp.model = self.model.value

        if self.mty_model_name is not None:
            cpp.mty_model_name = self.mty_model_name

        if self.mty_model_path is not None:
            cpp.mty_model_path = str(self.mty_model_path)

        if self.mty_bsm_mapping_path is not None:
            cpp.mty_bsm_mapping_path = str(self.mty_bsm_mapping_path)

        cpp.mty_order_policy = self.mty_order_policy.value
        if self.mty_tree_fermion_order is not None:
            order = [int(value) for value in self.mty_tree_fermion_order]
            if sorted(order) != list(range(len(order))):
                raise ValueError(
                    "mty_tree_fermion_order must be a permutation of 0..N-1; "
                    f"received {order}"
                )
            cpp.mty_tree_fermion_order = order
        return cpp

    def __repr__(self) -> str:
        """Return a compact representation suitable for logs."""
        return (
            "HyperisoConfig("
            f"model={self.model}, "
            f"flags={self.flags}, "
            f"mty_model_name={self.mty_model_name!r}, "
            f"mty_model_path={self.mty_model_path!r}, "
            f"mty_bsm_mapping_path={self.mty_bsm_mapping_path!r}, "
            f"mty_order_policy={self.mty_order_policy}, "
            f"mty_tree_fermion_order={self.mty_tree_fermion_order}"
            ")"
        )


__all__ = ["ExternalFlag", "HyperisoConfig"]
