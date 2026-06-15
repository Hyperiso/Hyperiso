"""Python wrapper for parameter and observable correlations."""

from __future__ import annotations

from enum import Enum

from pyhyperiso.phyperiso.pyhyperiso.core import CorrelationProvider as _CppCorrelationProvider
# from pyhyperiso.phyperiso.pyhyperiso.core import CorrelationType as _CppCorrelationType
_CppCorrelationType = _CppCorrelationProvider.CorrelationType
from pyhyperiso.core.Common.GeneralEnum import Observables
from pyhyperiso.core.Common.ParamId import ParamId


class CorrelationType(Enum):
    """Type of uncertainty correlation to query.

    Attributes:
        STAT: Statistical-only correlation.
        SYST: Systematic-only correlation.
        COMBINED: Correlation after combining statistical and systematic
            components according to the C++ provider rules.
    """

    STAT = _CppCorrelationType.STAT
    SYST = _CppCorrelationType.SYST
    COMBINED = _CppCorrelationType.COMBINED


class CorrelationProvider:
    """Read correlations stored by the C++ core.

    The provider supports two correlation domains:

    * parameter correlations, addressed with ``ParamId``;
    * experimental observable correlations, addressed with ``Observables``.

    Examples:
        >>> corr = CorrelationProvider()
        >>> corr.correlation_from_observable(
        ...     Observables.BR_B_XS_GAMMA,
        ...     Observables.BR_B_XS_GAMMA,
        ...     CorrelationType.COMBINED,
        ... )
        1.0
    """

    def __init__(self) -> None:
        """Create a provider bound to the current C++ correlation storage."""
        self._cpp_obj = _CppCorrelationProvider()

    def correlation_from_paramid(
        self,
        pid_1: ParamId,
        pid_2: ParamId,
        corr_type: CorrelationType,
    ) -> float:
        """Return the correlation between two parameters.

        Args:
            pid_1: First parameter identifier.
            pid_2: Second parameter identifier.
            corr_type: Statistical, systematic or combined correlation type.

        Returns:
            Correlation coefficient, typically in ``[-1, 1]``.
        """
        return float(self._cpp_obj(pid_1._cpp_obj, pid_2._cpp_obj, corr_type.value))

    def correlation_from_observable(
        self,
        obs_1: Observables,
        obs_2: Observables,
        corr_type: CorrelationType,
        experiment : str = "DEFAULT",
    ) -> float:
        """Return the experimental correlation between two observables.

        Args:
            obs_1: First observable enum.
            obs_2: Second observable enum.
            corr_type: Statistical, systematic or combined correlation type.

        Returns:
            Correlation coefficient, typically in ``[-1, 1]``.
        """
        return float(self._cpp_obj(experiment, obs_1.value, obs_2.value, corr_type.value))


__all__ = ["CorrelationProvider", "CorrelationType"]
        
if __name__ == "__main__":
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pathlib import Path
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Core.ParamaterProvider import ParameterProvider
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

    print("🔧 PyHyperisoConfig content:")
    print(config)

    hyp = HyperisoMaster()
    lha_file_path = "lha/si_input.flha" 

    print("\n🚀 Calling init with config...")
    hyp.init(lha_file=lha_file_path, config=config)
    
    corr_provider = CorrelationProvider()
    
    print("correlation between two obs (same obs)", corr_provider.correlation_from_observable(Observables.BR_B_XS_GAMMA, Observables.BR_B_XS_GAMMA, CorrelationType.COMBINED))