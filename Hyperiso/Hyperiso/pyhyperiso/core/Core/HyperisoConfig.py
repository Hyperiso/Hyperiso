from dataclasses import dataclass, field
from typing import Dict, Optional
from pathlib import Path
from enum import Enum
from pyhyperiso.phyperiso.pyhyperiso.core import HyperisoConfig as _CppHyperisoConfig
from pyhyperiso.phyperiso.pyhyperiso.core import ExternalFlag as _CppExternalFlag
from pyhyperiso.core.Common.GeneralEnum import Model

class ExternalFlag(Enum):
    IS_LHA_SPECTRUM = _CppExternalFlag.IS_LHA_SPECTRUM
    HAS_WILSON_INPUT = _CppExternalFlag.HAS_WILSON_INPUT
    HAS_TH_OBSERVABLE_INPUT = _CppExternalFlag.HAS_TH_OBSERVABLE_INPUT
    HYP_AS_SM_MARTY = _CppExternalFlag.HYP_AS_SM_MARTY
    # USE_MARTY = _CppExternalFlag.USE_MARTY
    

@dataclass
class HyperisoConfig:
    """Python wrapper for the C++ Config struct."""
    flags: Dict[ExternalFlag, bool] = field(default_factory=lambda: {
        ExternalFlag.IS_LHA_SPECTRUM: False,
        ExternalFlag.HAS_WILSON_INPUT: False,
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        ExternalFlag.HYP_AS_SM_MARTY: True,
        # ExternalFlag.USE_MARTY: False,
    })
    model: Model = Model.SM
    mty_model_name: Optional[str] = None
    mty_model_path: Optional[Path] = None

    def to_cpp(self) -> _CppHyperisoConfig:
        """Converts this Python wrapper to a C++ Config object."""
        cpp = _CppHyperisoConfig()
        cpp.flags = {flag.value: val for flag, val in self.flags.items()}
        cpp.model = self.model.value

        if self.mty_model_name is not None:
            cpp.mty_model_name = self.mty_model_name

        if self.mty_model_path is not None:
            cpp.mty_model_path = str(self.mty_model_path)

        return cpp
    
    def __repr__(self):
        return f"HyperisoConfig(model={self.model}, flags={self.flags})"
