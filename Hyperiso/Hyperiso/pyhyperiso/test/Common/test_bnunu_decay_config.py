"""Regression tests for B -> K(*) nu nu form-factor configuration."""

from pyhyperiso.Observable import (
    BKnunuConfig,
    BKstarnunuConfig,
    BPFFSource,
    BVFFSource,
)


def test_bnunu_form_factor_config_defaults_roundtrip():
    bk = BKnunuConfig()
    bk_cpp = bk.to_cpp()
    assert bk.ff_src is BPFFSource.AS
    assert bk_cpp.ff_src == BPFFSource.AS.value
    assert BKnunuConfig.from_cpp(bk_cpp).ff_src is BPFFSource.AS

    bkstar = BKstarnunuConfig()
    bkstar_cpp = bkstar.to_cpp()
    assert bkstar.ff_src is BVFFSource.BSZ_SR_LAT
    assert bkstar_cpp.ff_src == BVFFSource.BSZ_SR_LAT.value
    assert BKstarnunuConfig.from_cpp(bkstar_cpp).ff_src is BVFFSource.BSZ_SR_LAT


def test_bnunu_form_factor_config_custom_sources_roundtrip():
    bk = BKnunuConfig(ff_src=BPFFSource.FLAG24)
    assert BKnunuConfig.from_cpp(bk.to_cpp()).ff_src is BPFFSource.FLAG24

    bkstar = BKstarnunuConfig(ff_src=BVFFSource.GRvDV)
    assert BKstarnunuConfig.from_cpp(bkstar.to_cpp()).ff_src is BVFFSource.GRvDV
