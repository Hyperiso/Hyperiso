"""Example: make custom lambda observables visible to Statistic dependencies."""

from pyhyperiso.Common import Model, QCDOrder, ContributionType, ParameterType, ParamId
from pyhyperiso.Common import GroupMapper, WCoefMapper, ObservableMapper, LhaID
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Wilson import CustomWilsonCoefficientConfig, CustomWilsonGroupConfig
from pyhyperiso.Observable import ObservableInterface, LambdaDecayConfig, LambdaObservableConfig
from pyhyperiso.Statistic import StatisticInterface, StatisticConfig


def scalar_real(x):
    """Convert a Hyperiso scalar_t-like object to a Python float."""
    return x.real() if hasattr(x, "real") else float(x)


if __name__ == "__main__":

    config = HyperisoConfig(model=Model.SM)
    hyp = HyperisoMaster()
    hyp.init(lha_file="lha/si_input.flha", config=config)

    # Three ordinary input parameters. The custom Wilson lambdas and observables
    # declare them explicitly, so the statistic layer can discover them.
    f_bd = ParamId(ParameterType.FLAVOR, "FCONST", [511, 1])
    f_bs = ParamId(ParameterType.FLAVOR, "FCONST", [531, 1])
    m_top = ParamId(ParameterType.SM, "SMINPUTS", 6)

    GroupMapper.register_custom("PY_STAT_WILSON", aliases=["py-stat-wilson"])
    WCoefMapper.register_custom("C_PY_STAT_A", aliases=["c_py_stat_a"], flha=(997001, 1))
    WCoefMapper.register_custom("C_PY_STAT_B", aliases=["c_py_stat_b"], flha=(997001, 2))

    group = GroupMapper.id_of("py-stat-wilson")
    c_a = WCoefMapper.id_of("c_py_stat_a")
    c_b = WCoefMapper.id_of("c_py_stat_b")

    wc_a = CustomWilsonCoefficientConfig(c_a)
    wc_a.set_matching(
        QCDOrder.LO,
        {f_bd, m_top},
        lambda src: 0.2 + 0.5 * scalar_real(src.get_val(f_bd)),
        ContributionType.SM,
    )

    wc_b = CustomWilsonCoefficientConfig(c_b)
    wc_b.set_matching(
        QCDOrder.LO,
        {f_bs, m_top},
        lambda src: -0.1 + 0.3 * scalar_real(src.get_val(f_bs)),
        ContributionType.SM,
    )

    wc_group = CustomWilsonGroupConfig(group)
    wc_group.add_coefficient(wc_a).add_coefficient(wc_b)

    def obs_a(ctx, obs_id):
        c = ctx.get_FR(group, c_a, QCDOrder.LO, ContributionType.SM)
        f = ctx.get_flavor_param(f_bd)
        return scalar_real(c) + scalar_real(f)

    def obs_b(ctx, obs_id):
        c = ctx.get_FR(group, c_b, QCDOrder.LO, ContributionType.SM)
        f = ctx.get_flavor_param(f_bs)
        return scalar_real(c) - 0.5 * scalar_real(f)

    lambda_obs_a = LambdaObservableConfig.scalar(
        "PY_STAT_OBS_A",
        obs_a,
        aliases=["py-stat-obs-a"],
        flha=LhaID(998001, 1),
        dependencies={f_bd},
    )

    lambda_obs_b = LambdaObservableConfig.scalar(
        "PY_STAT_OBS_B",
        obs_b,
        aliases=["py-stat-obs-b"],
        flha=LhaID(998001, 2),
        dependencies={f_bs},
    )

    decay = LambdaDecayConfig(
        canonical="PY_STAT_DECAY",
        aliases=["py-stat-decay"],
        custom_wilson_groups=[wc_group],
        observables=[lambda_obs_a, lambda_obs_b],
        propagate_custom_wilson_dependencies=True,
    )

    oi = ObservableInterface()
    oi.add_lambda_decay(decay, add_observables=True)
    oi.enable_obs()

    id_a = ObservableMapper.id_of("py-stat-obs-a")
    id_b = ObservableMapper.id_of("py-stat-obs-b")

    print("Observable deps A:", oi.get_all_ops_deps_id(id_a))
    print("Observable deps B:", oi.get_all_ops_deps_id(id_b))
    print("Predictions:", oi.compute_all())

    stat = StatisticInterface(StatisticConfig(), oi)
    print("Statistic-visible dependencies:", stat.get_active_observable_dependencies())
