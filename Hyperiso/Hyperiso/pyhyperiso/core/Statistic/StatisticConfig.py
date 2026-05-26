from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.GeneralEnum import QCDOrder
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.MarginalConfig import MarginalKind


class StatisticLikelihoodMode(Enum):
    PROFILED_NUISANCE = st.StatisticLikelihoodMode.PROFILED_NUISANCE
    CHI2_MC_COVARIANCE = st.StatisticLikelihoodMode.CHI2_MC_COVARIANCE

    def to_cpp(self):
        return self.value


def _require(value, typ, name: str):
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_param_id(pid: ParamId):
    return _require(pid, ParamId, "ParamId").to_cpp()


def _cpp_binned_observable_id(obs: BinnedObservableId):
    return _require(obs, BinnedObservableId, "BinnedObservableId").to_cpp()


def _cpp_marginal_kind(kind: MarginalKind):
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_copula_kind(kind: CopulaKind):
    return _require(kind, CopulaKind, "CopulaKind").to_cpp()


def _cpp_likelihood_mode(mode: StatisticLikelihoodMode):
    return _require(mode, StatisticLikelihoodMode, "StatisticLikelihoodMode").to_cpp()


def _cpp_experiment_obs(obs: ExperimentObs):
    return _require(obs, ExperimentObs, "ExperimentObs").to_cpp()


@dataclass
class StatisticConfig:
    # Python-side only : sélection pratique pour tes workflows.
    obss: Dict[BinnedObservableId, QCDOrder] = field(default_factory=dict)
    p_specs: List[ParamId] = field(default_factory=list)
    selected_experiments: Optional[Sequence[str]] = None

    override_nuisance_marginals: Dict[ParamId, MarginalKind] = field(default_factory=dict)
    override_exp_data_marginals: Dict[ExperimentObs, MarginalKind] = field(default_factory=dict)

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

    likelihood_mode: StatisticLikelihoodMode = StatisticLikelihoodMode.PROFILED_NUISANCE
    chi2_covariance_ridge_rel: float = 1e-8
    chi2_covariance_ridge_abs: float = 1e-12

    nuisance_sensitivity_contexts: int = 2
    nuisance_sensitivity_context_sigma: float = 0.35
    nuisance_sensitivity_seed: int = 12345
    nuisance_sensitivity_keep_on_failure: bool = True

    def to_cpp(self):
        cpp = st.StatisticConfig()

        cpp.override_nuisance_marginals = {
            _cpp_param_id(pid): _cpp_marginal_kind(kind)
            for pid, kind in self.override_nuisance_marginals.items()
        }
        cpp.override_exp_data_marginals = {
            _cpp_experiment_obs(obs): _cpp_marginal_kind(kind)
            for obs, kind in self.override_exp_data_marginals.items()
        }

        cpp.nuisance_copula_type = _cpp_copula_kind(self.nuisance_copula_type)
        cpp.exp_data_copula_type = _cpp_copula_kind(self.exp_data_copula_type)

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

        cpp.likelihood_mode = _cpp_likelihood_mode(self.likelihood_mode)
        cpp.chi2_covariance_ridge_rel = float(self.chi2_covariance_ridge_rel)
        cpp.chi2_covariance_ridge_abs = float(self.chi2_covariance_ridge_abs)
        cpp.nuisance_sensitivity_contexts = int(self.nuisance_sensitivity_contexts)
        cpp.nuisance_sensitivity_context_sigma = float(self.nuisance_sensitivity_context_sigma)
        cpp.nuisance_sensitivity_seed = int(self.nuisance_sensitivity_seed)
        cpp.nuisance_sensitivity_keep_on_failure = bool(self.nuisance_sensitivity_keep_on_failure)

        return cpp

    def selected_observables_to_cpp(self):
        return {_cpp_binned_observable_id(obs): _require(order, QCDOrder, "QCDOrder").value for obs, order in self.obss.items()}

    def p_specs_to_cpp(self):
        return [_cpp_param_id(pid) for pid in self.p_specs]


__all__ = ["StatisticConfig", "StatisticLikelihoodMode"]

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