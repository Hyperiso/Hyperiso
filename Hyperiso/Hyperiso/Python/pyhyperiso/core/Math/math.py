"""
Python wrapper for C++ bindings generated using Pybind11.
Each function corresponds to a mathematical or numerical utility available in the C++ backend.
"""
from pyhyperiso.phyperiso.pyhyperiso import math as mb
from pyhyperiso.core.Math.scalar import Scalar
# Math functions

def Li2(x: float) -> float:
    """Compute the dilogarithm function.

    Args:
        x: A real number.

    Returns:
        The value of Li2(x).
    """
    return mb.Li2(x)

def Li3(x: float) -> float:
    """Compute the trilogarithm function.

    Args:
        x: A real number.

    Returns:
        The value of Li3(x).
    """
    return mb.Li3(x)

def CLi2(x: Scalar) -> complex:
    """Compute the complex dilogarithm function.

    Args:
        x: A complex number.

    Returns:
        Complex dilogarithm value.
    """
    if type(x) == complex:
        x = Scalar.from_complex(x)

    return Scalar.from_cpp(mb.CLi2(x._cpp_obj))

def Cl2(x: float) -> float:
    """Compute the Clausen function Cl2.

    Args:
        x: A real number.

    Returns:
        Value of Cl2(x).
    """
    return mb.Cl2(x)

def H2(x: float, y: float) -> float:
    """Compute the H2(x, y) special function.

    Args:
        x: First input value.
        y: Second input value.

    Returns:
        Value of H2(x, y).
    """
    return mb.H2(x, y)

def B(m1: float, m2: float, Q: float) -> float:
    """Compute the B-function of two masses and a scale Q.

    Args:
        m1: First mass.
        m2: Second mass.
        Q: Renormalization scale.

    Returns:
        Value of B(m1, m2, Q).
    """
    return mb.B(m1, m2, Q)

def kron(x: int, y: int) -> int:
    """Kronecker delta function.

    Args:
        x: First integer.
        y: Second integer.

    Returns:
        1 if x == y, else 0.
    """
    return mb.kron(x, y)

def integrate(f, l: float, u: float, prec: float = 1e-6) -> float:
    """Integrate a real-valued function f from l to u.

    Args:
        f: Callable real function.
        l: Lower integration bound.
        u: Upper integration bound.
        prec: Desired precision.

    Returns:
        Integral value.
    """
    return mb.integrate(f, l, u, prec)

def c_integrate(f, l: float, u: float, prec: float = 1e-6) -> Scalar:
    """Integrate a complex-valued function f over [l, u].

    Args:
        f: Callable complex function.
        l: Lower bound.
        u: Upper bound.
        prec: Precision.

    Returns:
        Complex integral value.
    """
    return Scalar.from_cpp(mb.c_integrate(f, l, u, prec))

def psi(n: int) -> float:
    """Compute the digamma function ψ(n) for integer input.

    Args:
        n: Integer input (n > 0).

    Returns:
        Value of ψ(n).
    """
    return mb.psi(n)

# Constants
PI = mb.PI
E = mb.E
ZETA3 = mb.ZETA3
GAMMA = mb.GAMMA

# Wilson coefficients
def A0t(x: float) -> float:
    """Wilson coefficient A0t.

    Args:
        x: A real parameter.

    Returns:
        Value of A0t(x).
    """
    return mb.A0t(x)

def F0t(x: float) -> float:
    """Wilson coefficient F0t."""
    return mb.F0t(x)

def B0t(x: float) -> float:
    """Wilson coefficient B0t."""
    return mb.B0t(x)

def C0t(x: float) -> float:
    """Wilson coefficient C0t."""
    return mb.C0t(x)

def D0t(x: float) -> float:
    """Wilson coefficient D0t."""
    return mb.D0t(x)

def E0t(x: float) -> float:
    """Wilson coefficient E0t."""
    return mb.E0t(x)

def T(x: float) -> float:
    """Wilson coefficient T."""
    return mb.T(x)

def A1t(x: float, l: float) -> float:
    """Wilson coefficient A1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of A1t(x, l).
    """
    return mb.A1t(x, l)

def B1t(x: float, l: float) -> float:
    """Wilson coefficient B1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of B1t(x, l).
    """
    return mb.B1t(x, l)

def C1t(x: float, l: float) -> float:
    """Wilson coefficient C1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of C1t(x, l).
    """
    return mb.C1t(x, l)

def D1t(x: float, l: float) -> float:
    """Wilson coefficient D1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of D1t(x, l).
    """
    return mb.D1t(x, l)

def E1t(x: float, l: float) -> float:
    """Wilson coefficient E1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of E1t(x, l).
    """
    return mb.E1t(x, l)

def F1t(x: float, l: float) -> float:
    """Wilson coefficient F1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of F1t(x, l).
    """
    return mb.F1t(x, l)

def G1t(x: float, l: float) -> float:
    """Wilson coefficient G1t.

    Args:
        x: A real parameter.
        l: Logarithmic scale.

    Returns:
        Value of G1t(x, l).
    """
    return mb.G1t(x, l)

def C7t2mt(x: float) -> float:
    """Wilson coefficient C7t at 2m_t scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C7t2mt(x).
    """
    return mb.C7t2mt(x)

def C7c2MW(x: float) -> float:
    """Wilson coefficient C7c at M_W scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C7c2MW(x).
    """
    return mb.C7c2MW(x)

def C8t2mt(x: float) -> float:
    """Wilson coefficient C8t at 2m_t scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C8t2mt(x).
    """
    return mb.C8t2mt(x)

def C8c2MW(x: float) -> float:
    """Wilson coefficient C8c at M_W scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C8c2MW(x).
    """
    return mb.C8c2MW(x)

def F7_1(x: float) -> float:
    """Wilson coefficient F7_1.

    Args:
        x: A real parameter.

    Returns:
        Value of F7_1(x).
    """
    return mb.F7_1(x)

def F7_2(x: float) -> float:
    """Wilson coefficient F7_2.

    Args:
        x: A real parameter.

    Returns:
        Value of F7_2(x).
    """
    return mb.F7_2(x)

def F8_1(x: float) -> float:
    """Wilson coefficient F8_1.

    Args:
        x: A real parameter.

    Returns:
        Value of F8_1(x).
    """
    return mb.F8_1(x)

def F8_2(x: float) -> float:
    """Wilson coefficient F8_2.

    Args:
        x: A real parameter.

    Returns:
        Value of F8_2(x).
    """
    return mb.F8_2(x)

def G3H(x: float, lu: float) -> float:
    """Wilson coefficient G3H.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of G3H(x, lu).
    """
    return mb.G3H(x, lu)

def G4H(x: float, lu: float) -> float:
    """Wilson coefficient G4H.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of G4H(x, lu).
    """
    return mb.G4H(x, lu)

def G7H(x: float, lu: float, ld: float) -> float:
    """Wilson coefficient G7H.

    Args:
        x: A real parameter.
        lu: Upper log scale.
        ld: Lower log scale.

    Returns:
        Value of G7H(x, lu, ld).
    """
    return mb.G7H(x, lu, ld)

def G8H(x: float, lu: float, ld: float) -> float:
    """Wilson coefficient G8H.

    Args:
        x: A real parameter.
        lu: Upper log scale.
        ld: Lower log scale.

    Returns:
        Value of G8H(x, lu, ld).
    """
    return mb.G8H(x, lu, ld)

def EH(x: float, lu: float) -> float:
    """Wilson coefficient EH.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of EH(x, lu).
    """
    return mb.EH(x, lu)

def D9H0(x: float, lu: float) -> float:
    """Wilson coefficient D9H0.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of D9H0(x, lu).
    """
    return mb.D9H0(x, lu)

def D9H1(x: float, lu: float, L: float) -> float:
    """Wilson coefficient D9H1.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.
        L: Additional parameter.

    Returns:
        Value of D9H1(x, lu, L).
    """
    return mb.D9H1(x, lu, L)

def Delta3H(x: float, lu: float) -> float:
    """Wilson coefficient Delta3H.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of Delta3H(x, lu).
    """
    return mb.Delta3H(x, lu)

def Delta4H(x: float, lu: float) -> float:
    """Wilson coefficient Delta4H.

    Args:
        x: A real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of Delta4H(x, lu).
    """
    return mb.Delta4H(x, lu)

def Delta7H(x: float, lu: float, ld: float) -> float:
    """Wilson coefficient Delta7H.

    Args:
        x: A real parameter.
        lu: Upper log scale.
        ld: Lower log scale.

    Returns:
        Value of Delta7H(x, lu, ld).
    """
    return mb.Delta7H(x, lu, ld)

def Delta8H(x: float, lu: float, ld: float) -> float:
    """Wilson coefficient Delta8H.

    Args:
        x: A real parameter.
        lu: Upper log scale.
        ld: Lower log scale.

    Returns:
        Value of Delta8H(x, lu, ld).
    """
    return mb.Delta8H(x, lu, ld)

def C9llH0(x: float, y: float, lu: float) -> float:
    """Wilson coefficient C9llH0.

    Args:
        x: First real parameter.
        y: Second real parameter.
        lu: Logarithmic scale.

    Returns:
        Value of C9llH0(x, y, lu).
    """
    return mb.C9llH0(x, y, lu)

def C9llH1(x: float, y: float, lu: float, L: float) -> float:
    """Wilson coefficient C9llH1.

    Args:
        x: First real parameter.
        y: Second real parameter.
        lu: Logarithmic scale.
        L: Additional scale.

    Returns:
        Value of C9llH1(x, y, lu, L).
    """
    return mb.C9llH1(x, y, lu, L)

def C10Wt2mt(x: float) -> float:
    """Wilson coefficient C10Wt at 2m_t scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C10Wt2mt(x).
    """
    return mb.C10Wt2mt(x)

def C10Wc2MW(x: float) -> float:
    """Wilson coefficient C10Wc at M_W scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C10Wc2MW(x).
    """
    return mb.C10Wc2MW(x)

def C10Zt2mt(x: float) -> float:
    """Wilson coefficient C10Zt at 2m_t scale.

    Args:
        x: A real parameter.

    Returns:
        Value of C10Zt2mt(x).
    """
    return mb.C10Zt2mt(x)

def C10Z2tri(x: float) -> float:
    """Wilson coefficient C10Z from triangle diagrams.

    Args:
        x: A real parameter.

    Returns:
        Value of C10Z2tri(x).
    """
    return mb.C10Z2tri(x)

def F0SP(xt: float) -> float:
    """Wilson coefficient F0SP.

    Args:
        xt: A real parameter.

    Returns:
        Value of F0SP(xt).
    """
    return mb.F0SP(xt)

