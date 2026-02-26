# anomalies_plot_hyperiso.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any, List
import numpy as np
import matplotlib.pyplot as plt

# ---- Hyperiso imports (adapte si besoin selon ton arborescence) ----
from pyhyperiso.core.Common.GeneralEnum import Observables, QCDOrder, ParameterType, DataType, Model
from pyhyperiso.core.Common.ParamId import ParamId, LhaID
from pyhyperiso.core.Common.SymbolId import ObservableId
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId

from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
from pyhyperiso.core.Core.HyperisoConfig import PyHyperisoConfig, ExternalFlag
from pyhyperiso.core.Core.ParamaterProvider import PyParameterProvider

from pyhyperiso.core.BusinessLogic.ObservableInterface import PyObservableInterface
from pyhyperiso.core.Common.Mapper import ObservableMapper  # <-- ton mapper static
from pyhyperiso.core.Statistic.StatisticInterface import StatisticInterface
from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary


# -------------------------------------------------------------------
# Helpers float conversions
# -------------------------------------------------------------------
def _scalar_to_float(x: Any) -> float:
    # ton Scalar wrapper a souvent .value
    if hasattr(x, "value"):
        return float(x.value)
    return float(x)


# -------------------------------------------------------------------
# Plot "Point" structure (analogue au Q2point de Patrick, version minimaliste)
# -------------------------------------------------------------------
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


# -------------------------------------------------------------------
# Build exp/theory points
# -------------------------------------------------------------------
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


# -------------------------------------------------------------------
# muPoints + drawTwoPoints (pull plot)
# -------------------------------------------------------------------
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
        if (not flip and exp.obs < theo.obs):
            newE.obs = -newE.obs

    return newE, newT


def draw_two_points(ax, y: float, exp_and_theo):
    yo = 0.04
    exp_pt, theo_pt = exp_and_theo
    ax.errorbar(y=y - yo, x=theo_pt.obs, xerr=theo_pt.totError(), color="orange", marker="d", linestyle="")
    ax.errorbar(y=y + yo, x=exp_pt.obs, xerr=exp_pt.totError(), color="b", marker="o", linestyle="")


# -------------------------------------------------------------------
# Normalize summaries keys (important si les clés viennent en type C++ pybind)
# -------------------------------------------------------------------
def normalize_summaries_keys(summaries: Dict[Any, GaussianSummary]) -> Dict[BinnedObservableId, GaussianSummary]:
    out: Dict[BinnedObservableId, GaussianSummary] = {}
    for k, v in summaries.items():
        if isinstance(k, BinnedObservableId):
            out[k] = v
        else:
            # typiquement: k est un _CppObservableId ou un truc str()-able
            out[BinnedObservableId(str(k))] = v
    return out


# -------------------------------------------------------------------
# Main plotting function
# -------------------------------------------------------------------
def anomalies_plot(
    # hyp_interface: PyObservableInterface,
    # stat_interface: StatisticInterface,
    # exp_provider: PyParameterProvider,
    flip: bool = True,
):
    # Observables.PHI
    selected = [
        Observables.BR_BS_MUMU,
        Observables.BR_BD_MUMU,
        Observables.BR_BU_TAU_NU,
        Observables.R_D,
        Observables.R_DSTAR
    ]
    selected = [BinnedObservableId(ObservableMapper.to_id(x)) for x in selected]
    
    selected.append(BinnedObservableId(ObservableId(ObservableMapper.str(Observables.DBR_DQ2_B__K_MU_MU)), (1.1, 6)))
    # 2) Statistic interface (incertitude théorie)
    stat_cfg = StatisticConfig(
        obss={o: QCDOrder.NNLO for o in selected},
        p_specs=[],
        MC_draws=200,  # tu peux monter
        skew_abs_threshold=0.2,
    )
    stat_interface = StatisticInterface(stat_cfg)

    # 3) Exp provider (FOBS)
    exp_provider = PyParameterProvider(ParameterType.OBSERVABLE)
    
    # # 1) incertitudes théoriques
    summaries_raw = stat_interface.compute_uncertainties()
    print("raw : ", summaries_raw)
    summaries = normalize_summaries_keys(summaries_raw)
    # summaries = {BinnedObservableId(gs._cpp_type): GaussianSummary.from_cpp(gs) for gs in summaries}
    # obs_interface = PyObservableInterface()
    

    # for o in selected:
    #     obs_interface.add_observable(o, QCDOrder.NNLO, add_dependencies=True)

    

    labels: List[str] = [""]
    points = []  # (label, expPoint, theoPoint)

    def add_nonbinned(obs: Observables, latex_label: str):
        # EXP: via FOBS / LhaID
        obs_2 = BinnedObservableId(ObservableMapper().to_id(obs))
        lhaid: LhaID = obs_2.flha()
        exp = exp_point_from_fobs(exp_provider, lhaid)

        # obs_interface.add_observable(obs, QCDOrder.NNLO)
        print(obs)
        # print("aaah : ", obs_interface._cpp_obj.get_current_observables())
        # THEO central
        # theo_val = hyp_interface.compute_observable_central(obs)
        # theo_val = obs_interface.compute_observable(obs)
        print("the sum : ", summaries)
        
        bid = BinnedObservableId(ObservableMapper.to_id(obs))   # (0,0) par défaut
        theo_val = summaries[bid].mu
        print(theo_val)
        # THEO uncertainty id
        gid: ObservableId = ObservableMapper.id_of(ObservableMapper.str(obs))
        gobs_2 = BinnedObservableId(gid)
        glhaid: LhaID = gobs_2.flha()
        if gobs_2 not in summaries:
            keys_preview = list(summaries.keys())[:10]
            raise KeyError(
                f"[anomalies_plot] ObservableId {gid} absent de summaries.\n"
                f"Obs={obs} str={ObservableMapper.str(obs)}\n"
                f"Exemples de clés présentes: {keys_preview}\n"
                f"Ca arrive si StatisticConfig.obss n'inclut pas cet observable, "
                f"ou si l'id construit n'est pas exactement celui utilisé côté C++."
            )

        g = summaries[gobs_2]
        th = theory_point_from_stats(theo_val, g)
        # th = Point()
        # th.obs = theo_val[0].value
        labels.append(latex_label)
        # points.append((latex_label, exp, th))
        points.append((latex_label, exp, th))

    # ---- Ajoute ici ta sélection d'observables (PK-style) ----
    add_nonbinned(Observables.DBR_DQ2_B__K_MU_MU, r"${\cal B}(B^+\toK^+\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BS_MUMU, r"${\cal B}(B_s^0\to\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BD_MUMU, r"${\cal B}(B^0\to\mu^+\mu^-)$")
    add_nonbinned(Observables.BR_BU_TAU_NU, r"${\cal B}(B^+\to\tau^+\nu)$")
    add_nonbinned(Observables.R_D, r"$R_D$")
    add_nonbinned(Observables.R_DSTAR, r"$R_D^*$")
    # Ajoute les autres de ton choix...

    # ---- Plot ----
    fig, ax = plt.subplots()
    fig.subplots_adjust(top=0.98, right=0.98, bottom=0.1, left=0.42)

    # y: du haut vers le bas, façon PK
    for i, (_, exp, th) in enumerate(points, start=1):
        exp_pull, th_pull = mu_points(exp, th, flip=flip)
        y = len(labels) - 1 - i
        draw_two_points(ax, y=y, exp_and_theo=(exp_pull, th_pull))

    # style
    nsigma = 6
    labels.append("")
    # ax.set_yticks(ax.get_yticks().tolist())
    # ax.set_yticklabels(labels)
    n = len(points)
    ypos = np.arange(n)

    ax.set_yticks(ypos)
    ax.set_yticklabels([lbl for (lbl, _, _) in points])
    ax.invert_yaxis()
    ax.grid(axis="x")
    ax.set_xlabel(r"Pull in $\sigma$")
    ax.set_xlim((-0.99, nsigma) if flip else (-nsigma, nsigma))

    return fig


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
def main():
    # 0) init Hyperiso
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
        # mty_model_path=Path("/my/custom/marty/path"),
    )

    lha_file_path = "lha/si_input.flha"  # <-- adapte
    hyp.init(lha_file=lha_file_path, config=config)

    # 1) Observable interface (théorie)
    # obs_interface = PyObservableInterface()
    # selected = [
    #     Observables.BR_BS_MUMU,
    #     Observables.BR_BD_MUMU,
    #     Observables.BR_BU_TAU_NU,
    # ]

    # for o in selected:
    #     obs_interface.add_observable(o, QCDOrder.NNLO, add_dependencies=True)

    # # 2) Statistic interface (incertitude théorie)
    # stat_cfg = StatisticConfig(
    #     obss={o: QCDOrder.NNLO for o in selected},
    #     p_specs=[],
    #     MC_draws=200,  # tu peux monter
    #     skew_abs_threshold=0.2,
    # )
    # stat_interface = StatisticInterface(stat_cfg)

    # # 3) Exp provider (FOBS)
    # exp_provider = PyParameterProvider(ParameterType.OBSERVABLE)

    # 4) Plot
    fig = anomalies_plot(
        # hyp_interface=obs_interface,
        # stat_interface=stat_interface,
        # exp_provider=exp_provider,
        flip=True,
    )

    plt.show()
    # fig.savefig("anomalies_plot.pdf")


if __name__ == "__main__":
    main()
