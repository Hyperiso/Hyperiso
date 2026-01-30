import pytest
from pyhyperiso.core.Common.Configs import (
    WilsonBuildConfig,
    WilsonRequest,
    AlphasConfig,
    MassConfig
)
from pyhyperiso.core.Common.GeneralEnum import (
    WGroup, WCoeff, QCDOrder, ContributionType, ScaleType, MassType
)
from pyhyperiso.core.Common.Mapper import GroupMapper

def test_wilson_build_config_to_cpp():
    config = WilsonBuildConfig(
        groups={WGroup.B, WGroup.BScalar},
        matching_scale=1000.0,
        hadronic_scale=2.0,
        order=QCDOrder.NLO
    )

    cpp_config = config.to_cpp()
    assert cpp_config.matching_scale == 1000.0
    assert cpp_config.hadronic_scale == 2.0
    assert cpp_config.order == QCDOrder.NLO.value
    # assert GroupMapper().id_of(WGroup.B) in cpp_config.groups #TODO : check
    # assert GroupMapper().id_of(WGroup.BScalar) in cpp_config.groups


def test_wilson_request_to_cpp():
    request = WilsonRequest(
        group=WGroup.B,
        coefficient=WCoeff.C9,
        order=QCDOrder.LO,
        contribution=ContributionType.BSM,
        scale_type=ScaleType.MATCHING,
        sum_qcd_orders=True
    )

    cpp_request = request.to_cpp()
    assert cpp_request.group == WGroup.B.value
    assert cpp_request.coefficient == WCoeff.C9.value
    assert cpp_request.order == QCDOrder.LO.value
    assert cpp_request.contribution == ContributionType.BSM.value
    assert cpp_request.scale_type == ScaleType.MATCHING.value
    assert cpp_request.sum_qcd_orders is True


def test_alphas_config_to_cpp():
    config = AlphasConfig(
        scale=91.1876,
        m_b_type=MassType.MSBAR,
        m_t_type=MassType.POLE
    )
    cpp = config.to_cpp()
    assert cpp.scale == pytest.approx(91.1876)
    assert cpp.m_b_type == MassType.MSBAR.value
    assert cpp.m_t_type == MassType.POLE.value


def test_mass_config_to_cpp():
    config = MassConfig(
        pdg_id=6,
        scale=173.0,
        m_b_type=MassType.POLE,
        m_t_type=MassType.MSBAR
    )
    cpp = config.to_cpp()
    assert cpp.pdg_id == 6
    assert cpp.scale == pytest.approx(173.0)
    assert cpp.m_b_type == MassType.POLE.value
    assert cpp.m_t_type == MassType.MSBAR.value
