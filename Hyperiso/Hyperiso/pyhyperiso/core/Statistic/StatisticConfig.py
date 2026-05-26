"""High-level configuration objects for the statistic workflow.

``StatisticConfig`` is the Python entry point for controlling the C++ statistic
manager. It groups observable selections, fit-parameter selections, nuisance
marginal overrides, copula choices, Monte-Carlo settings and maximum-likelihood
options.
"""

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
    """Likelihood backend used by ``StatisticManager.compute_MLE``.

    Attributes:
        PROFILED_NUISANCE: Full likelihood with explicit nuisance parameters,
            nuisance copula and experimental-observable copula.
        CHI2_MC_COVARIANCE: Faster chi-square backend with no explicit nuisance
            coordinates. The covariance is estimated from Monte-Carlo theory
            propagation and combined with the experimental covariance.
    """

    PROFILED_NUISANCE = st.StatisticLikelihoodMode.PROFILED_NUISANCE
    CHI2_MC_COVARIANCE = st.StatisticLikelihoodMode.CHI2_MC_COVARIANCE

    def to_cpp(self):
        """Convert to the bound C++ ``StatisticLikelihoodMode`` enum."""
        return self.value


def _require(value, typ, name: str):
    """Validate a typed wrapper argument."""
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_param_id(pid: ParamId):
    """Convert a Python parameter id to C++."""
    return _require(pid, ParamId, "ParamId").to_cpp()


def _cpp_binned_observable_id(obs: BinnedObservableId):
    """Convert a Python binned observable id to C++."""
    return _require(obs, BinnedObservableId, "BinnedObservableId").to_cpp()


def _cpp_marginal_kind(kind: MarginalKind):
    """Convert a Python marginal kind to C++."""
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_copula_kind(kind: CopulaKind):
    """Convert a Python copula kind to C++."""
    return _require(kind, CopulaKind, "CopulaKind").to_cpp()


def _cpp_likelihood_mode(mode: StatisticLikelihoodMode):
    """Convert a Python likelihood mode to C++."""
    return _require(mode, StatisticLikelihoodMode, "StatisticLikelihoodMode").to_cpp()


def _cpp_experiment_obs(obs: ExperimentObs):
    """Convert a Python experimental observable wrapper to C++."""
    return _require(obs, ExperimentObs, "ExperimentObs").to_cpp()


@dataclass
class StatisticConfig:
    """Configuration of the Python/C++ statistic pipeline.

    Attributes:
        obss: Python-side observable selection. Keys are binned observables and
            values are QCD orders. This is convenient for higher-level workflows;
            the C++ ``StatisticConfig`` itself does not store this field.
        p_specs: Fit parameters to use by default when calling
            ``StatisticInterface.compute_MLE`` or likelihood-scan methods without
            an explicit parameter list.
        selected_experiments: Optional initial experiment selection. When set,
            ``StatisticInterface`` applies it after constructing the C++ manager.
        override_nuisance_marginals: Per-parameter overrides for nuisance
            marginal families.
        override_exp_data_marginals: Per-observable overrides for experimental
            data marginal families.
        nuisance_copula_type: Copula used for nuisance-parameter correlations.
        exp_data_copula_type: Copula used for experimental-observable
            correlations.
        MC_draws: Number of Monte-Carlo accepted predictions used in uncertainty
            propagation.
        skew_abs_threshold: Skewness threshold used to choose symmetric versus
            split-Gaussian summaries.
        MLE_max_iter: Maximum number of function calls/iterations passed to the
            minimization backend.
        MLE_tol: Minimizer tolerance.
        MLE_strategy: Backend minimization strategy. For the Minuit backend,
            ``2`` is a more robust strategy.
        MLE_run_hesse: Whether to ask the backend for a covariance/HESSE step.
        MLE_request_minos: Whether to request MINOS when supported by the build.
        MLE_verbose: Whether to run the backend in verbose mode.
        nuisance_relevance_cutoff: Relative uncertainty prefilter threshold for
            nuisance candidates.
        nuisance_sensitivity_pruning: Whether to prune nuisance candidates by
            finite-difference model sensitivity.
        nuisance_sensitivity_probe_sigmas: Step size, in nuisance standard
            deviations, used by the sensitivity probe.
        nuisance_sensitivity_rel_cutoff: Relative observable-shift threshold for
            keeping a nuisance parameter.
        nuisance_sensitivity_abs_cutoff: Absolute observable-shift threshold for
            keeping a nuisance parameter.
        nuisance_sensitivity_scale_floor: Minimum scale used when normalizing
            relative observable shifts.
        MLE_trace_first_evals: Whether the likelihood should print the first NLL
            evaluations for debugging.
        MLE_trace_max_evals: Maximum number of traced evaluations.
        MLE_allow_profile_hessian_fallback: Whether the fitter may reconstruct
            a profile covariance numerically when the backend covariance is
            unavailable.
        MLE_profile_hessian_step_scale: Finite-difference step scale used by the
            profile-Hessian fallback.
        MLE_profile_hessian_eig_floor_rel: Relative eigenvalue floor used when
            regularizing the profile Hessian.
        likelihood_mode: Likelihood backend used by maximum-likelihood fits.
        chi2_covariance_ridge_rel: Relative ridge used by the chi-square
            covariance backend.
        chi2_covariance_ridge_abs: Absolute ridge used by the chi-square
            covariance backend.
        nuisance_sensitivity_contexts: Number of nuisance contexts used by
            sensitivity pruning. A negative value disables the check in C++.
        nuisance_sensitivity_context_sigma: Width, in nuisance sigma units, used
            to randomize non-central sensitivity contexts.
        nuisance_sensitivity_seed: Random seed for sensitivity contexts.
        nuisance_sensitivity_keep_on_failure: Whether a nuisance is kept when its
            sensitivity evaluation fails.

    Examples:
        >>> cfg = StatisticConfig(
        ...     MC_draws=1000,
        ...     likelihood_mode=StatisticLikelihoodMode.CHI2_MC_COVARIANCE,
        ... )
        >>> cpp_cfg = cfg.to_cpp()
    """

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
        """Convert this Python config to the bound C++ ``StatisticConfig``.

        Returns:
            A bound C++ configuration object. Python-only fields such as
            ``obss``, ``p_specs`` and ``selected_experiments`` are intentionally
            not stored in the C++ object.
        """
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
        """Convert the Python-side observable selection to a C++ mapping.

        Returns:
            Mapping from C++ binned observable identifiers to C++ QCD order enum
            values.
        """
        return {_cpp_binned_observable_id(obs): _require(order, QCDOrder, "QCDOrder").value for obs, order in self.obss.items()}

    def p_specs_to_cpp(self):
        """Convert the default fit-parameter list to C++ identifiers.

        Returns:
            List of C++ ``ParamId`` objects.
        """
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