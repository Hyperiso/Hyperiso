"""Configuration objects used to initialize Hyperiso."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Union

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
        mty_tree_fermion_orders: Per-coefficient TreeLevel fermion-order overrides.
        mty_one_loop_fermion_orders: Per-coefficient OneLoop fermion-order overrides.

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
    mty_tree_fermion_orders: Mapping[str, Sequence[int]] = field(default_factory=dict)
    mty_one_loop_fermion_orders: Mapping[str, Sequence[int]] = field(default_factory=dict)

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
        cpp.mty_tree_fermion_orders = self._normalise_fermion_orders(
            self.mty_tree_fermion_orders,
            field_name="mty_tree_fermion_orders",
        )
        cpp.mty_one_loop_fermion_orders = self._normalise_fermion_orders(
            self.mty_one_loop_fermion_orders,
            field_name="mty_one_loop_fermion_orders",
        )
        return cpp

    @staticmethod
    def _normalise_fermion_orders(
        orders: Mapping[str, Sequence[int]], *, field_name: str
    ) -> Dict[str, list[int]]:
        """Validate and copy a coefficient-to-permutation mapping."""
        normalised: Dict[str, list[int]] = {}
        for raw_coefficient, raw_order in orders.items():
            coefficient = str(raw_coefficient).strip()
            if not coefficient:
                raise ValueError(f"{field_name} contains an empty coefficient name")
            order = [int(value) for value in raw_order]
            if not order:
                raise ValueError(
                    f"{field_name}[{coefficient!r}] must not be empty"
                )
            if sorted(order) != list(range(len(order))):
                raise ValueError(
                    f"{field_name}[{coefficient!r}] must be a permutation of "
                    f"0..N-1; received {order}"
                )
            normalised[coefficient] = order
        return normalised

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
            f"mty_tree_fermion_orders={dict(self.mty_tree_fermion_orders)!r}, "
            f"mty_one_loop_fermion_orders={dict(self.mty_one_loop_fermion_orders)!r}"
            ")"
        )


__all__ = ["ExternalFlag", "HyperisoConfig"]
