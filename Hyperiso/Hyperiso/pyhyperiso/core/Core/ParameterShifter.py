"""Mutable parameter shifter for nuisance/scan workflows."""

from __future__ import annotations

from typing import Union

from pyhyperiso.phyperiso.pyhyperiso.core import ParameterShifter as _CppParameterShifter
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.Parameter import ParameterMode
from pyhyperiso.core.Math.scalar import Scalar, _to_scalar


class ParameterShifter:
    """Python wrapper around the C++ ``ParameterShifter``.

    This class exposes the same Python surface as ``ParameterSetter`` but is
    backed by the C++ shifter object. It is useful for workflows that treat
    parameter changes as shifts rather than direct value replacement in the C++
    parameter layer.
    """

    def __init__(self) -> None:
        """Create a shifter bound to the current C++ parameter storage."""
        self._cpp_obj = _CppParameterShifter()

    def mutate(self, pid: ParamId, value: Union[float, complex, Scalar]) -> None:
        """Apply a new shifted value to a parameter.

        Args:
            pid: Parameter to shift.
            value: Shifted value. Numeric inputs are converted to ``Scalar``.
        """
        self._cpp_obj.mutate(pid._cpp_obj, _to_scalar(value)._cpp_obj)

    def change_mode(self, pid: ParamId, mode: ParameterMode) -> None:
        """Change a parameter update mode.

        Args:
            pid: Parameter to update.
            mode: New mode, such as ``ParameterMode.FIXED`` or
                ``ParameterMode.SHIFTABLE``.
        """
        self._cpp_obj.change_mode(pid._cpp_obj, mode.value)


__all__ = ["ParameterShifter"]

        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
    from pyhyperiso.core.Common.GeneralEnum import Model, ParameterType
    
    print("🔧 Initializing PyHyperisoMaster with custom PyHyperisoConfig...")

    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/camilia.flha"

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    param_shifter = ParameterShifter()
    param_prov = ParameterProvider(ParameterType.SM)
    
    print("param before mutate : ", param_prov.get_by_block("MASS", 24))
    param_shifter.mutate(ParamId(ParameterType.SM, "MASS", 24), 82)
    print("param before mutate : ", param_prov.get_by_block("MASS", 24))
    