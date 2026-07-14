import pytest
from pyhyperiso.core.Math.Scalar import Scalar, sqrt, exp, sin, pow_scalar


def test_scalar_creation():
    s = Scalar(1.5, -0.5)
    print(s, s.real(), s.imag())
    assert s.real() == pytest.approx(1.5)
    assert s.imag() == pytest.approx(-0.5)

def test_scalar_from_complex():
    z = complex(2, 3)
    s = Scalar.from_complex(z)
    assert s.real() == pytest.approx(2.0)
    assert s.imag() == pytest.approx(3.0)


def test_scalar_add_sub():
    a = Scalar(1, 2)
    b = Scalar(3, -1)
    c = a + b
    d = a - b
    assert c.real() == pytest.approx(4)
    assert c.imag() == pytest.approx(1)
    assert d.real() == pytest.approx(-2)
    assert d.imag() == pytest.approx(3)


def test_scalar_mul_div():
    a = Scalar(1, 1)
    b = Scalar(1, -1)
    c = a * b
    d = a / b
    assert c.real() == pytest.approx(2)
    assert c.imag() == pytest.approx(0)
    assert d.real() == pytest.approx(0)
    assert d.imag() == pytest.approx(1)


def test_scalar_casts():
    s = Scalar(4.5)
    assert float(s) == pytest.approx(4.5)
    assert complex(s) == pytest.approx(complex(4.5, 0.0))


def test_scalar_math_functions():
    s = Scalar(0.0, 1.0)
    r = sqrt(s)
    assert isinstance(r, Scalar)
    assert r.real() > 0 or r.imag() > 0

    s2 = exp(Scalar(0.0, 3.141592653589793))  # e^(iπ) ≈ -1
    assert s2.real() == pytest.approx(-1.0, abs=1e-5)

    s3 = sin(Scalar(0))
    assert s3.real() == pytest.approx(0.0)
    assert s3.imag() == pytest.approx(0.0)


def test_scalar_pow():
    base = Scalar(2)
    exp_s = Scalar(3)
    result = pow_scalar(base, exp_s)
    assert result.real() == pytest.approx(8.0)

    result2 = pow_scalar(base, 4)
    assert result2.real() == pytest.approx(16.0)
