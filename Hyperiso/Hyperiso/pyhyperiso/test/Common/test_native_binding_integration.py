"""Integration smoke tests for the compiled pybind11 extension."""

from pyhyperiso.phyperiso import pyhyperiso as native


def test_native_module_imports_and_reports_release_version():
    assert native.__version__ == "1.0.1"
    assert native.common is not None
    assert native.core is not None
    assert native.statistic is not None


def test_block_name_canonical_alias_is_stable():
    block = native.common.BlockName("mass")
    block.add_alias("MASS_ALIAS")

    assert block.canonical() == "mass"
    assert str(block) == "mass"

    block.to_upper()
    assert block.canonical() == "MASS"
    assert str(block) == "MASS"
