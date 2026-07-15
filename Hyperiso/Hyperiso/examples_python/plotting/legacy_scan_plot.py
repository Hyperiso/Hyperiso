"""Plot normalized experimental and theoretical pulls for selected observables."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List
import numpy as np
import matplotlib.pyplot as plt

from pyhyperiso.core.Common.GeneralEnum import Observables, QCDOrder, ParameterType, Model
from pyhyperiso.core.Common.ParamId import ParamId, LhaID
from pyhyperiso.core.Common.SymbolId import ObservableId
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId

from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
from pyhyperiso.core.Core.ParameterProvider import PyParameterProvider

from pyhyperiso.core.Common.Mapper import ObservableMapper
from pyhyperiso.core.Statistic.StatisticInterface import StatisticInterface
from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


def _scalar_to_float(x: Any) -> float:
    if hasattr(x, "value"):
        return float(x.value)
    return float(x)


@dataclass
class Point:
    q2low: float = 0.0
    q2high: float = 0.0
    obs: float = 0.0
    StatUp: float = 0.0
    StatDown: float = 0.0
    SystUp: float = 0.0
    SystDown: float = 0.0
    significance: Optional[float] = None

    def totError(self) -> Tuple[Tuple[float], Tuple[float]]:
        up = np.sqrt(self.StatUp**2 + self.SystUp**2)
        dn = np.sqrt(self.StatDown**2 + self.SystDown**2)
        return (up,), (dn,)


def exp_point_from_fobs(provider: PyParameterProvider, lhaid: LhaID) -> Point:
    pid = ParamId(type=ParameterType.OBSERVABLE, block="FOBS", code=lhaid)
    p = provider.get_parameter(pid)
    stat, syst = p.std

    return Point(
        obs=_scalar_to_float(p.value),
        StatUp=_scalar_to_float(stat),
        StatDown=_scalar_to_float(stat),
        SystUp=_scalar_to_float(syst),
        SystDown=_scalar_to_float(syst),
    )


def theory_point_from_stats(theo_value: float, g: GaussianSummary) -> Point:
    if g.symmetric:
        up = dn = float(g.sigma)
    else:
        up = float(g.sigma_p) if float(g.sigma_p) != 0.0 else float(g.sigma)
        dn = float(g.sigma_m) if float(g.sigma_m) != 0.0 else float(g.sigma)

    return Point(obs=float(theo_value), StatUp=up, StatDown=dn)


def mu_points(exp: Point, theo: Point, flip: bool = True):
    theoOffset = theo.obs

    if exp.obs > theo.obs or not flip:
        scale = np.sqrt(exp.totError()[0][0] ** 2 + theo.totError()[1][0] ** 2)
        newT = Point(obs=0.0, StatUp=theo.StatUp / scale, StatDown=theo.StatDown / scale)
        newE = Point(
            obs=(exp.obs - theoOffset) / scale,
            StatUp=exp.StatUp / scale,
            StatDown=exp.StatDown / scale,
            SystUp=exp.SystUp / scale,
            SystDown=exp.SystDown / scale,
            significance=exp.significance,
        )
    else:
        scale = np.sqrt(exp.totError()[1][0] ** 2 + theo.totError()[0][0] ** 2)
        newT = Point(obs=0.0, StatUp=theo.StatDown / scale, StatDown=theo.StatUp / scale)
        newE = Point(
            obs=-(exp.obs - theoOffset) / scale,
            StatUp=exp.StatDown / scale,
            StatDown=exp.StatUp / scale,
            SystUp=exp.SystDown / scale,
            SystDown=exp.SystUp / scale,
            significance=exp.significance,
        )

    if exp.significance is not None:
        newE.obs = float(exp.significance)
        if not flip and exp.obs < theo.obs:
            newE.obs = -newE.obs

    return newE, newT


def draw_two_points(ax, y: float, exp_and_theo):
    yo = 0.04
    exp_pt, theo_pt = exp_and_theo
    ax.errorbar(
        y=y - yo, x=theo_pt.obs, xerr=theo_pt.totError(), color="orange", marker="d", linestyle=""
    )
    ax.errorbar(y=y + yo, x=exp_pt.obs, xerr=exp_pt.totError(), color="b", marker="o", linestyle="")


def normalize_summaries_keys(
    summaries: Dict[Any, GaussianSummary],
) -> Dict[BinnedObservableId, GaussianSummary]:
    out: Dict[BinnedObservableId, GaussianSummary] = {}
    for k, v in summaries.items():
        if isinstance(k, BinnedObservableId):
            out[k] = v
        else:
            out[BinnedObservableId(str(k))] = v
    return out


def anomalies_plot(
    flip: bool = True,
):
    selected = [
        Observables.BR_BS_MUMU,
        Observables.BR_BD_MUMU,
        Observables.BR_BU_TAU_NU,
        Observables.R_D,
        Observables.R_DSTAR,
    ]
    selected = [BinnedObservableId(ObservableMapper.to_id(x)) for x in selected]

    selected.append(
        BinnedObservableId(
            ObservableId(ObservableMapper.str(Observables.DBR_DQ2_B__K_MU_MU)), (1.1, 6)
        )
    )
    stat_cfg = StatisticConfig(
        obss={o: QCDOrder.NNLO for o in selected},
        p_specs=[],
        MC_draws=200,
        skew_abs_threshold=0.2,
    )
    stat_interface = StatisticInterface(stat_cfg)

    exp_provider = PyParameterProvider(ParameterType.OBSERVABLE)

    summaries = normalize_summaries_keys(stat_interface.compute_uncertainties())

    labels: List[str] = [""]
    points = []

    def add_nonbinned(obs: Observables, latex_label: str):
        obs_2 = BinnedObservableId(ObservableMapper().to_id(obs))
        lhaid: LhaID = obs_2.flha()
        exp = exp_point_from_fobs(exp_provider, lhaid)

        bid = BinnedObservableId(ObservableMapper.to_id(obs))
        theo_val = summaries[bid].mu
        gid: ObservableId = ObservableMapper.id_of(ObservableMapper.str(obs))
        gobs_2 = BinnedObservableId(gid)
        if gobs_2 not in summaries:
            keys_preview = list(summaries.keys())[:10]
            raise KeyError(
                f"ObservableId {gid} is absent from the uncertainty summaries.\n"
                f"Observable={obs}, mapped name={ObservableMapper.str(obs)}\n"
                f"Available-key preview: {keys_preview}\n"
                "Ensure StatisticConfig.obss contains the observable and that the "
                "constructed identifier matches the C++ identifier."
            )

        g = summaries[gobs_2]
        th = theory_point_from_stats(theo_val, g)

        labels.append(latex_label)
        points.append((latex_label, exp, th))

    add_nonbinned(Observables.DBR_DQ2_B__K_MU_MU, r"${\cal B}(B^+\toK^+\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BS_MUMU, r"${\cal B}(B_s^0\to\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BD_MUMU, r"${\cal B}(B^0\to\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BU_TAU_NU, r"${\cal B}(B^+\to\tau^+\nu)$")
    add_nonbinned(Observables.R_D, r"$R_D$")
    add_nonbinned(Observables.R_DSTAR, r"$R_D^*$")

    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.98, right=0.98, bottom=0.1, left=0.42)

    for i, (_, exp, th) in enumerate(points, start=1):
        exp_pull, th_pull = mu_points(exp, th, flip=flip)
        y = len(labels) - 1 - i
        draw_two_points(ax, y=y, exp_and_theo=(exp_pull, th_pull))

    nsigma = 6
    labels.append("")

    n = len(points)
    ypos = np.arange(n)

    ax.set_yticks(ypos)
    ax.set_yticklabels([lbl for (lbl, _, _) in points])
    ax.invert_yaxis()
    ax.grid(axis="x")
    ax.set_xlabel(r"Pull in $\sigma$")
    ax.set_xlim((-0.99, nsigma) if flip else (-nsigma, nsigma))

    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default_lha = (
        Path(__file__).resolve().parents[2] / "pyhyperiso" / "assets" / "lha" / "si_input.flha"
    )
    parser.add_argument(
        "--lha",
        type=Path,
        default=default_lha,
        help="FLHA input file (default: packaged si_input.flha asset)",
    )
    args = parser.parse_args()

    hyp = PyHyperisoMaster()

    config = PyHyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True,
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
    )

    hyp.init(lha_file=str(args.lha), config=config)

    anomalies_plot(
        flip=True,
    )

    plt.show()


if __name__ == "__main__":
    main()
