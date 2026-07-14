from pyhyperiso.core.Math.Math import *
from pyhyperiso.core.Math.Scalar import Scalar
import math
import cmath


def test_kron():
    assert kron(1, 1) == 1
    assert kron(2, 3) == 0


def test_Li2():
    assert math.isclose(Li2(0), 0.0, abs_tol=1e-10)
    assert math.isclose(Li2(1), math.pi**2 / 6, rel_tol=1e-5)


def test_Li3():
    assert math.isclose(Li3(0), 0.0, abs_tol=1e-10)


def test_CLi2():
    z = complex(0.5, 0.5)
    assert isinstance(CLi2(z), Scalar)


def test_Cl2():
    assert math.isclose(Cl2(0), 0.0, abs_tol=1e-10)


def test_H2():
    assert isinstance(H2(1.0, 2.0), float)


def test_B():
    assert isinstance(B(1.0, 2.0, 1.0), float)


def test_integrate():
    result = integrate(lambda x: x**2, 0, 1)
    assert math.isclose(result, 1 / 3, rel_tol=1e-4)


def test_c_integrate():
    result = c_integrate(lambda x: Scalar(x, 0)._cpp_obj, 0, 1)
    assert cmath.isclose(Scalar.from_cpp(result), Scalar(0.5, 0), abs_tol=1e-4)


def test_psi():
    assert math.isclose(psi(1), -0.5772156649, rel_tol=1e-5)


def test_A0t():
    assert isinstance(A0t(1.0), float)


def test_B0t():
    assert isinstance(B0t(1.0), float)


def test_C0t():
    assert isinstance(C0t(1.0), float)


def test_D0t():
    assert isinstance(D0t(1.0), float)


def test_E0t():
    assert isinstance(E0t(1.0), float)


def test_F0t():
    assert isinstance(F0t(1.0), float)


def test_T():
    assert isinstance(T(1.0), float)


def test_A1t():
    assert isinstance(A1t(1.0, 0.1), float)


def test_B1t():
    assert isinstance(B1t(1.0, 0.1), float)


def test_C1t():
    assert isinstance(C1t(1.0, 0.1), float)


def test_D1t():
    assert isinstance(D1t(1.0, 0.1), float)


def test_E1t():
    assert isinstance(E1t(1.0, 0.1), float)


def test_F1t():
    assert isinstance(F1t(1.0, 0.1), float)


def test_G1t():
    assert isinstance(G1t(1.0, 0.1), float)


def test_C7C8():
    assert isinstance(C7t2mt(1.0), float)
    assert isinstance(C7c2MW(1.0), float)
    assert isinstance(C8t2mt(1.0), float)
    assert isinstance(C8c2MW(1.0), float)


def test_F78():
    assert isinstance(F7_1(1.0), float)
    assert isinstance(F7_2(1.0), float)
    assert isinstance(F8_1(1.0), float)
    assert isinstance(F8_2(1.0), float)


def test_GH_EH():
    assert isinstance(G3H(1.0, 0.1), float)
    assert isinstance(G4H(1.0, 0.1), float)
    assert isinstance(G7H(1.0, 0.1, 0.05), float)
    assert isinstance(G8H(1.0, 0.1, 0.05), float)
    assert isinstance(EH(1.0, 0.1), float)


def test_Deltas():
    assert isinstance(D9H0(1.0, 0.1), float)
    assert isinstance(D9H1(1.0, 0.1, 0.5), float)
    assert isinstance(Delta3H(1.0, 0.1), float)
    assert isinstance(Delta4H(1.0, 0.1), float)
    assert isinstance(Delta7H(1.0, 0.1, 0.05), float)
    assert isinstance(Delta8H(1.0, 0.1, 0.05), float)


def test_C9ll_C10():
    assert isinstance(C9llH0(1.0, 2.0, 0.1), float)
    assert isinstance(C9llH1(1.0, 2.0, 0.1, 0.5), float)
    assert isinstance(C10Wt2mt(1.0), float)
    assert isinstance(C10Wc2MW(1.0), float)
    assert isinstance(C10Zt2mt(1.0), float)
    assert isinstance(C10Z2tri(1.0), float)
    assert isinstance(F0SP(1.0), float)
