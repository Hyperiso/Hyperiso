"""Python wrapper for immutable QCD constants exposed by the C++ core.

This module provides a small read-only facade over the bound C++
``QCDConstants`` object. The object is usually obtained through
:class:`pyhyperiso.core.Core.QCDProvider.QCDProvider` and contains standard
color factors and perturbative coefficients used by the QCD running routines.

Example:
    >>> provider = QCDProvider()
    >>> constants = provider.get_qcd_constants()
    >>> constants.Nc
    3
"""

from pyhyperiso.phyperiso.pyhyperiso.core import QCDConstants as _CppQCDConstants


class QCDConstants:
    """Read-only view of QCD color factors and perturbative coefficients.

    The wrapped C++ object is not constructed directly from Python. It is
    returned by ``QCDProvider.get_qcd_constants()`` after the Hyperiso core has
    been initialized.

    Args:
        cpp_obj: Bound C++ ``QCDConstants`` instance.

    Attributes:
        Nc: Number of colors, usually ``3``.
        C_F: Fundamental Casimir color factor.
        C_A: Adjoint Casimir color factor.
        beta: Beta-function coefficients used by the QCD running backend.
        gamma: Anomalous-dimension coefficients used by the QCD running backend.
    """

    def __init__(self, cpp_obj: _CppQCDConstants):
        self._cpp_obj = cpp_obj

    @property
    def Nc(self):
        """Number of QCD colors.

        Returns:
            int | float: Value stored by the C++ backend.
        """
        return self._cpp_obj.Nc

    @property
    def C_F(self):
        """Fundamental Casimir color factor.

        Returns:
            float: The coefficient commonly denoted :math:`C_F`.
        """
        return self._cpp_obj.C_F

    @property
    def C_A(self):
        """Adjoint Casimir color factor.

        Returns:
            float: The coefficient commonly denoted :math:`C_A`.
        """
        return self._cpp_obj.C_A

    @property
    def beta(self):
        """QCD beta-function coefficients.

        Returns:
            Sequence[float]: Coefficients exposed by the C++ QCD provider.
        """
        return self._cpp_obj.beta

    @property
    def gamma(self):
        """QCD anomalous-dimension coefficients.

        Returns:
            Sequence[float]: Coefficients exposed by the C++ QCD provider.
        """
        return self._cpp_obj.gamma

    def __repr__(self):
        """Return a compact developer-oriented representation."""
        return (f"QCDConstants(Nc={self.Nc}, C_F={self.C_F}, C_A={self.C_A}, "
                f"beta={self.beta}, gamma={self.gamma})")


__all__ = ["QCDConstants"]
