"""High-level Python access to QED/electroweak coupling quantities from the C++ core.

The :class:`QEDProvider` wrapper delegates to the C++ ``QEDProvider`` and
exposes the electromagnetic coupling computation backed by ``EWHelper``.
"""

from pyhyperiso.phyperiso.pyhyperiso.core import QEDProvider as _CppQEDProvider
from pyhyperiso.core.Common.Configs import AlphasConfig


class QEDProvider:
    r"""Compute QED/electroweak quantities through the bound C++ provider.

    The provider expects the global Hyperiso state to have been initialized
    before use, typically through ``HyperisoMaster.init(...)``. The numerical
    implementation is delegated to the C++ ``EWHelper`` backend.

    Note:
        At the C++ level, generic running through ``EWHelper::alpha_em``
        raises an explicit error because it is not implemented. The precomputed
        values of :math:`\alpha_{\mathrm{em}}` are stored in the ``EW`` block
        initialized by ``EWHelper::Init()``.

    Example:
        >>> qed = QEDProvider()
        >>> alpha_em = qed.get_alphaem(AlphasConfig(91.1876))
    """

    def __init__(self):
        """Create a Python wrapper around the C++ ``QEDProvider``."""
        self._cpp_obj = _CppQEDProvider()

    def get_alphaem(self, alpha_config: AlphasConfig):
        """Compute the electromagnetic coupling for the requested configuration.

        Args:
            alpha_config: Python configuration object describing the target
                scale and convertible to the C++ ``AlphasConfig``.

        Returns:
            float: The value returned by ``QEDProvider::compute_alphaem``.
        """
        return self._cpp_obj.compute_alphaem(alpha_config.to_cpp())

    def get_alpha_em(self, alpha_config: AlphasConfig):
        """Alias for :meth:`get_alphaem` using the conventional alpha_em spelling."""
        return self.get_alphaem(alpha_config)


__all__ = ["QEDProvider"]
