"""Example: define a custom decay and observable from Python lambdas."""

from pyhyperiso.Common import Model, QCDOrder, ContributionType, WilsonBasis
from pyhyperiso.Common import GroupMapper, WCoefMapper, ObservableMapper
from pyhyperiso.Core import HyperisoConfig, HyperisoMaster
from pyhyperiso.Wilson import CustomWilsonCoefficientConfig, CustomWilsonGroupConfig
from pyhyperiso.Observable import ObservableInterface, LambdaDecayConfig, LambdaObservableConfig


def scalar_real(x):
    """Convert a Hyperiso scalar_t-like object to a Python float."""
    return x.real() if hasattr(x, "real") else float(x)


if __name__ == "__main__":

    config = HyperisoConfig(model=Model.SM)
    hyp = HyperisoMaster()
    hyp.init(lha_file="lha/si_input.flha", config=config)

    # 1) Define a custom Wilson coefficient used by the custom observable.
    GroupMapper.register_custom("PY_OBS_WILSON", aliases=["py-obs-wilson"])
    WCoefMapper.register_custom("C_PY_OBS", aliases=["c_py_obs"], flha=(996001, 1))

    group = GroupMapper.id_of("py-obs-wilson")
    coeff = WCoefMapper.id_of("c_py_obs")

    cobs = CustomWilsonCoefficientConfig(coeff)
    cobs.set_matching(QCDOrder.LO, set(), lambda src: 0.75, ContributionType.SM)

    wc_group = CustomWilsonGroupConfig(group, matching_scale=160.0, hadronic_scale=4.8)
    wc_group.add_coefficient(cobs)

    # 2) Define a custom decay with one custom observable.
    def compute_dev_obs(ctx, obs_id):
        # ctx gives access to dynamic Wilson getters, SM/FLAVOR parameters and bins.
        c = ctx.get_FM(group, coeff, QCDOrder.LO, ContributionType.SM)
        return 2.0 * scalar_real(c)

    obs = LambdaObservableConfig.scalar(
        "PY_LAMBDA_OBS",
        compute_dev_obs,
        aliases=["py-lambda-obs"],
    )

    decay = LambdaDecayConfig(
        canonical="PY_LAMBDA_DECAY",
        aliases=["py-lambda-decay"],
        custom_wilson_groups=[wc_group],
        observables=[obs],
        order=QCDOrder.LO,
    )

    interface = ObservableInterface()
    interface.add_lambda_decay(decay, add_observables=True)

    obs_id = ObservableMapper.id_of("py-lambda-obs")
    print(interface.compute_observable_id(obs_id))
