from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Common.GeneralEnum import QCDOrder
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.MarginalConfig import MarginalKind
from pyhyperiso.core.Statistic.Copula import CopulaKind


_CppStatisticConfig = st.StatisticConfig


def _to_cpp_obj(x: Any) -> Any:
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    if hasattr(x, "_cpp_obj"):
        return x._cpp_obj
    if hasattr(x, "_cpp"):
        return x._cpp
    return x


def _to_cpp_enum(x: Any) -> Any:
    if hasattr(x, "to_cpp") and callable(getattr(x, "to_cpp")):
        return x.to_cpp()
    return x


def _marginal_override_map_to_cpp(values: Mapping[Any, Any]) -> Dict[Any, Any]:
    return {
        _to_cpp_obj(k): _to_cpp_enum(v)
        for k, v in values.items()
    }


@dataclass
class StatisticConfig:
    # Conservé côté Python pour compatibilité avec ton ancienne API.
    # N'est plus envoyé à StatisticConfig C++.
    obss: Dict[BinnedObservableId, QCDOrder] = field(default_factory=dict)

    # Conservé côté Python comme liste par défaut pour compute_MLE().
    # N'est plus envoyé à StatisticConfig C++.
    p_specs: List[ParamId] = field(default_factory=list)

    override_nuisance_marginals: Dict[ParamId, MarginalKind] = field(default_factory=dict)
    override_exp_data_marginals: Dict[Any, MarginalKind] = field(default_factory=dict)

    nuisance_copula_type: CopulaKind = CopulaKind.GAUSSIAN
    exp_data_copula_type: CopulaKind = CopulaKind.GAUSSIAN

    MC_draws: int = 100
    skew_abs_threshold: float = 0.2

    MLE_max_iter: int = 500
    MLE_tol: float = 1e-8
    MLE_strategy: int = 2
    MLE_run_hesse: bool = True
    MLE_request_minos: bool = False
    MLE_verbose: bool = False

    nuisance_relevance_cutoff: float = 1e-8

    nuisance_sensitivity_pruning: bool = True
    nuisance_sensitivity_probe_sigmas: float = 1.0
    nuisance_sensitivity_rel_cutoff: float = 1e-6
    nuisance_sensitivity_abs_cutoff: float = 1e-12
    nuisance_sensitivity_scale_floor: float = 1e-3

    MLE_trace_first_evals: bool = False
    MLE_trace_max_evals: int = 25

    MLE_allow_profile_hessian_fallback: bool = True
    MLE_profile_hessian_step_scale: float = 1.0
    MLE_profile_hessian_eig_floor_rel: float = 1e-8

    # Option Python-side seulement, si tu exposes select_experiments côté C++.
    selected_experiments: Optional[Sequence[str]] = None

    def to_cpp(self) -> _CppStatisticConfig:
        cpp = _CppStatisticConfig()

        cpp.override_nuisance_marginals = _marginal_override_map_to_cpp(
            self.override_nuisance_marginals
        )
        cpp.override_exp_data_marginals = _marginal_override_map_to_cpp(
            self.override_exp_data_marginals
        )

        cpp.nuisance_copula_type = _to_cpp_enum(self.nuisance_copula_type)
        cpp.exp_data_copula_type = _to_cpp_enum(self.exp_data_copula_type)

        cpp.MC_draws = int(self.MC_draws)
        cpp.skew_abs_threshold = float(self.skew_abs_threshold)

        cpp.MLE_max_iter = int(self.MLE_max_iter)
        cpp.MLE_tol = float(self.MLE_tol)
        cpp.MLE_strategy = int(self.MLE_strategy)
        cpp.MLE_run_hesse = bool(self.MLE_run_hesse)
        cpp.MLE_request_minos = bool(self.MLE_request_minos)
        cpp.MLE_verbose = bool(self.MLE_verbose)

        cpp.nuisance_relevance_cutoff = float(self.nuisance_relevance_cutoff)

        cpp.nuisance_sensitivity_pruning = bool(self.nuisance_sensitivity_pruning)
        cpp.nuisance_sensitivity_probe_sigmas = float(self.nuisance_sensitivity_probe_sigmas)
        cpp.nuisance_sensitivity_rel_cutoff = float(self.nuisance_sensitivity_rel_cutoff)
        cpp.nuisance_sensitivity_abs_cutoff = float(self.nuisance_sensitivity_abs_cutoff)
        cpp.nuisance_sensitivity_scale_floor = float(self.nuisance_sensitivity_scale_floor)

        cpp.MLE_trace_first_evals = bool(self.MLE_trace_first_evals)
        cpp.MLE_trace_max_evals = int(self.MLE_trace_max_evals)

        cpp.MLE_allow_profile_hessian_fallback = bool(self.MLE_allow_profile_hessian_fallback)
        cpp.MLE_profile_hessian_step_scale = float(self.MLE_profile_hessian_step_scale)
        cpp.MLE_profile_hessian_eig_floor_rel = float(self.MLE_profile_hessian_eig_floor_rel)

        return cpp

    def __repr__(self) -> str:
        return (
            "StatisticConfig("
            f"MC_draws={self.MC_draws}, "
            f"skew_abs_threshold={self.skew_abs_threshold}, "
            f"MLE_max_iter={self.MLE_max_iter}, "
            f"MLE_tol={self.MLE_tol}, "
            f"MLE_strategy={self.MLE_strategy}, "
            f"MLE_run_hesse={self.MLE_run_hesse}, "
            f"MLE_request_minos={self.MLE_request_minos}, "
            f"MLE_verbose={self.MLE_verbose}, "
            f"nuisance_relevance_cutoff={self.nuisance_relevance_cutoff}, "
            f"nuisance_sensitivity_pruning={self.nuisance_sensitivity_pruning}, "
            f"p_specs={self.p_specs}"
            ")"
        )
if __name__ == "__main__":
    a = StatisticConfig()
    print(a)
    
    a.to_cpp()
    
    import pyhyperiso.phyperiso.pyhyperiso.common as cA
    import pyhyperiso.phyperiso.pyhyperiso.common as cB

    print("A:", cA.__file__)
    print("B:", cB.__file__)
    print("same QCDOrder object?", cA.QCDOrder is cB.QCDOrder)
    print("same LO?", cA.QCDOrder.LO == cB.QCDOrder.LO)