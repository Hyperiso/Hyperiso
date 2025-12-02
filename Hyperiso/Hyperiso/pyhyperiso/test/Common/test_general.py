from pyhyperiso.core.Common.General import ParameterType, PyParamId, PyLhaID
import pytest


def test_lhaid_from_string():
    lhaid = PyLhaID("1_2_3")
    assert lhaid.get_parts() == [1, 2, 3]
    assert str(lhaid) == "PyLhaID(1_2_3)"


def test_lhaid_from_list():
    lhaid = PyLhaID([4, 5])
    assert lhaid.get_parts() == [4, 5]


def test_paramid_default():
    pid = PyParamId()
    assert pid.type is None
    assert pid.block == "NULL"
    assert int(pid.code) == 0


def test_paramid_from_block_and_code():
    pid = PyParamId(block="BLOCK", code=PyLhaID(42))
    assert pid.block == "BLOCK"
    assert int(pid.code) == 42


def test_paramid_with_type():
    lhaid = PyLhaID("10_20")
    pid = PyParamId(ParameterType.BSM, "MSSM", lhaid)
    assert pid.type == ParameterType.BSM
    assert pid.block == "MSSM"
    assert pid.code.get_parts() == [10, 20]


def test_paramid_set_type():
    pid = PyParamId(block="FLAV", code=1)
    assert pid.type is None
    pid.set_parameter_type(ParameterType.FLAVOR)
    assert pid.type == ParameterType.FLAVOR


def test_paramid_repr_and_dict():
    pid = PyParamId(ParameterType.SM, "X", [1, 2])
    repr_str = repr(pid)
    assert "PyParamId" in repr_str
    d = pid.to_dict()
    assert d == {
        "type": "SM",
        "block": "X",
        "code": "1_2"
    }
