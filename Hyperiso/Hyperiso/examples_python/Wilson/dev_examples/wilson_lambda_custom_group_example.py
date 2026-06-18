"""Example: define a custom Wilson group from Python lambdas."""

from pyhyperiso.Common import Model, QCDOrder, ContributionType, WilsonBasis, ParameterType, ParamId
from pyhyperiso.Common import GroupMapper, WCoefMapper
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Wilson import WilsonInterface, CustomWilsonCoefficientConfig, CustomWilsonGroupConfig


def scalar_real(x):
    """Convert a Hyperiso scalar_t-like object to a Python float for toy formulas."""
    return x.real() if hasattr(x, "real") else float(x)


if __name__ == "__main__":

    config = HyperisoConfig(model=Model.SM)
    hyp = HyperisoMaster()
    hyp.init(lha_file="lha/si_input.flha", config=config)

    # Register runtime Wilson names. They do not need to exist in the static WGroup/WCoeff enums.
    GroupMapper.register_custom("PY_DEV_WILSON", aliases=["py-dev-wilson"])
    WCoefMapper.register_custom("C_PY_DEV_1", aliases=["c_py_dev_1"], flha=(995001, 1))
    WCoefMapper.register_custom("C_PY_DEV_2", aliases=["c_py_dev_2"], flha=(995001, 2))

    group = GroupMapper.id_of("py-dev-wilson")
    c1 = WCoefMapper.id_of("c_py_dev_1")
    c2 = WCoefMapper.id_of("c_py_dev_2")

    # Standard input dependency used by the first matching lambda.
    # Declaring it here means Statistic/Observable dependency tooling can see it later.
    m_top = ParamId(ParameterType.SM, "SMINPUTS", 6)

    # First coefficient: depends on mt, but only weakly in this toy example.
    coef1 = CustomWilsonCoefficientConfig(c1)
    coef1.set_matching(
        QCDOrder.LO,
        {m_top},
        lambda src: 1.25 + 0.0 * scalar_real(src.get_val(m_top)),
        ContributionType.SM,
    )

    # Second coefficient: no input dependency.
    coef2 = CustomWilsonCoefficientConfig(c2)
    coef2.set_matching(
        QCDOrder.LO,
        set(),
        lambda src: -0.40,
        ContributionType.SM,
    )

    wc_group = CustomWilsonGroupConfig(
        group,
        matching_scale=160.0,
        hadronic_scale=4.8,
        order=QCDOrder.LO,
        contribution=ContributionType.SM,
    )
    wc_group.add_coefficient(coef1).add_coefficient(coef2)

    # If no running lambda is supplied, the C++ code installs identity running.
    # Here we keep the example minimal and rely on that default.
    wilson = WilsonInterface()
    wilson.add_custom_group(wc_group)

    print("C_PY_DEV_1 matching LO =", wilson.get_M(group, c1, QCDOrder.LO, ContributionType.SM))
    print("C_PY_DEV_2 running LO  =", wilson.get_R(group, c2, QCDOrder.LO, ContributionType.SM, WilsonBasis.STANDARD))
