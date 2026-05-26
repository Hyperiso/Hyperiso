"""Mutable parameter setter for central values and modes."""

from __future__ import annotations

from typing import Union

from pyhyperiso.phyperiso.pyhyperiso.core import ParameterSetter as _CppParameterSetter
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.Parameter import ParameterMode
from pyhyperiso.core.Math.scalar import Scalar, _to_scalar


class ParameterSetter:
    """Python wrapper around the C++ ``ParameterSetter``.

    ``ParameterSetter`` mutates stored parameter values directly. For workflows
    that require recomputing observables after mutation, remember to call the
    relevant reload/update method on the observable or statistic interface.

    Examples:
        >>> setter = ParameterSetter()
        >>> setter.mutate(ParamId(ParameterType.SM, "MASS", 24), 80.379)
    """

    def __init__(self) -> None:
        """Create a setter bound to the current C++ parameter storage."""
        self._cpp_obj = _CppParameterSetter()

    def mutate(self, pid: ParamId, value: Union[float, complex, Scalar]) -> None:
        """Set a parameter value.

        Args:
            pid: Parameter to mutate.
            value: New value. Numeric inputs are converted to ``Scalar``.
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


__all__ = ["ParameterSetter"]
        
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
            # ExternalFlag.USE_MARTY: False
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
    
    param_shifter = ParameterSetter()
    param_prov = ParameterProvider(ParameterType.SM)
    
    print("param before mutate : ", param_prov.get_by_block("MASS", 24))
    param_shifter.mutate(ParamId(ParameterType.SM, "MASS", 24), 82)
    print("param before mutate at 82: ", param_prov.get_by_block("MASS", 24))
    