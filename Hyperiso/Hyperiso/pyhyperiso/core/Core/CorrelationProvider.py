"""Python wrapper for parameter and observable correlations."""

from __future__ import annotations

from enum import Enum

from pyhyperiso.core.Common.GeneralEnum import Observables
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.phyperiso.pyhyperiso.core import CorrelationProvider as _CppCorrelationProvider

_CppCorrelationType = _CppCorrelationProvider.CorrelationType


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
        experiment: str = "DEFAULT",
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
