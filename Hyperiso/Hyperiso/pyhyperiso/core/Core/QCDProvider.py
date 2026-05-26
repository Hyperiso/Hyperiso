"""High-level Python access to QCD running quantities from the C++ core.

The :class:`QCDProvider` wrapper delegates to the C++ ``QCDProvider`` and
exposes the few operations needed by Python workflows: running coupling,
running quark masses, and backend QCD constants.
"""

from pyhyperiso.phyperiso.pyhyperiso.core import QCDProvider as _CppQCDProvider
from pyhyperiso.core.Common.Configs import AlphasConfig, MassConfig
from pyhyperiso.core.Core.QCDConstants import QCDConstants


class QCDProvider:
    """Compute QCD quantities through the bound C++ provider.

    The provider expects the global Hyperiso state to have been initialized
    before use, typically through ``HyperisoMaster.init(...)``. The exact
    numerical scheme and thresholds are managed by the C++ backend.

    Example:
        >>> qcd = QCDProvider()
        >>> alpha_s = qcd.get_alphas(AlphasConfig(42.0))
        >>> mt_run = qcd.get_qcd_masses(MassConfig(5.0, pdg_id=6))
        >>> constants = qcd.get_qcd_constants()
    """

    def __init__(self):
        """Create a Python wrapper around the C++ ``QCDProvider``."""
        self._cpp_obj = _CppQCDProvider()

    def get_alphas(self, alpha_config: AlphasConfig):
        """Compute the strong coupling for the requested configuration.

        Args:
            alpha_config: Python configuration object describing the scale and
                any backend-specific options required by the C++ provider.

        Returns:
            float: The value returned by ``QCDProvider::compute_alphas``.
        """
        return self._cpp_obj.compute_alphas(alpha_config.to_cpp())

    def get_qcd_masses(self, mass_config: MassConfig):
        """Compute a running QCD mass for the requested configuration.

        Args:
            mass_config: Python configuration object containing the scale, PDG
                identifier, and mass convention expected by the C++ backend.

        Returns:
            float: The value returned by ``QCDProvider::compute_mass``.
        """
        return self._cpp_obj.compute_mass(mass_config.to_cpp())

    def get_qcd_constants(self) -> QCDConstants:
        """Return QCD constants used internally by the backend.

        Returns:
            QCDConstants: Read-only Python wrapper around the C++ constants
            object.
        """
        cpp_constants = self._cpp_obj.get_constants()
        return QCDConstants(cpp_constants)


__all__ = ["QCDProvider"]
        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model
    
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

    print("🔧 HyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/camilia.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    qcd_params = QCDProvider()

    print("alphas_s at 42 GeV : ", qcd_params.get_alphas(AlphasConfig(42)))
    
    print("mass at 5 GeV : ", qcd_params.get_qcd_masses(MassConfig(5, pdg_id=6)))
    
    consts = qcd_params.get_qcd_constants()

    print(consts.Nc)      # 3
    print(consts.C_F)     # 4/3
    print(consts.beta[0]) # par ex., [31/3, 134/3, ...]