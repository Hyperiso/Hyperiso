"""High-level configuration object for the statistic workflow.

``StatisticConfig`` is the Python entry point for configuring the bound C++
statistic manager. It mirrors the fields exposed by the pybind11
``StatisticConfig`` binding.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.MarginalConfig import MarginalKind


class StatisticLikelihoodMode(Enum):
    """Likelihood backend used by the statistic manager.

    Attributes:
        PROFILED_NUISANCE: Full likelihood backend with explicit nuisance
            parameters, nuisance copula and experimental-observable copula.
        CHI2_MC_COVARIANCE: Chi-square backend without explicit nuisance
            coordinates. The covariance is estimated from Monte-Carlo theory
            propagation and combined with the experimental covariance.
    """

    PROFILED_NUISANCE = st.StatisticLikelihoodMode.PROFILED_NUISANCE
    CHI2_MC_COVARIANCE = st.StatisticLikelihoodMode.CHI2_MC_COVARIANCE

    def to_cpp(self):
        """Convert this Python enum wrapper to the bound C++ enum.

        Returns:
            Bound C++ ``StatisticLikelihoodMode`` value.
        """
        return self.value


def _require(value, typ, name: str):
    """Validate that a value has the expected wrapper type.

    Args:
        value: Value to validate.
        typ: Expected Python type.
        name: Human-readable argument name.

    Returns:
        The validated value.

    Raises:
        TypeError: If ``value`` is not an instance of ``typ``.
    """
    if not isinstance(value, typ):
        raise TypeError(
            f"{name} must be {typ.__name__}, got {type(value).__name__}."
        )
    return value


def _cpp_param_id(pid: ParamId):
    """Convert a Python ``ParamId`` wrapper to C++.

    Args:
        pid: Python parameter identifier.

    Returns:
        Bound C++ ``ParamId`` value.
    """
    return _require(pid, ParamId, "ParamId").to_cpp()


def _cpp_marginal_kind(kind: MarginalKind):
    """Convert a Python ``MarginalKind`` wrapper to C++.

    Args:
        kind: Python marginal kind enum.

    Returns:
        Bound C++ ``MarginalKind`` value.
    """
    return _require(kind, MarginalKind, "MarginalKind").to_cpp()


def _cpp_copula_kind(kind: CopulaKind):
    """Convert a Python ``CopulaKind`` wrapper to C++.

    Args:
        kind: Python copula kind enum.

    Returns:
        Bound C++ ``CopulaKind`` value.
    """
    return _require(kind, CopulaKind, "CopulaKind").to_cpp()


def _cpp_likelihood_mode(mode: StatisticLikelihoodMode):
    """Convert a Python likelihood mode wrapper to C++.

    Args:
        mode: Python likelihood backend enum.

    Returns:
        Bound C++ ``StatisticLikelihoodMode`` value.
    """
    return _require(
        mode,
        StatisticLikelihoodMode,
        "StatisticLikelihoodMode",
    ).to_cpp()


def _cpp_experiment_obs(obs: ExperimentObs):
    """Convert a Python ``ExperimentObs`` wrapper to C++.

    Args:
        obs: Python experimental-observable wrapper.

    Returns:
        Bound C++ ``ExperimentObs`` value.
    """
    return _require(obs, ExperimentObs, "ExperimentObs").to_cpp()


@dataclass
class StatisticConfig:
    """Configuration of the bound C++ statistic pipeline.

    This dataclass mirrors the fields exposed by the pybind11
    ``StatisticConfig`` binding.

    Attributes:
        override_nuisance_marginals: Per-parameter overrides for nuisance
            marginal distributions.
        override_exp_data_marginals: Per-observable overrides for experimental
            data marginal distributions.
        nuisance_copula_type: Copula used for nuisance-parameter correlations.
        exp_data_copula_type: Copula used for experimental-observable
            correlations.
        MC_draws: Number of Monte-Carlo accepted predictions used in uncertainty
            propagation.
        skew_abs_threshold: Skewness threshold used to choose symmetric versus
            split-Gaussian summaries.
        MLE_max_iter: Maximum number of function calls or iterations passed to
            the minimization backend.
        MLE_tol: Minimizer tolerance.
        MLE_strategy: Backend minimization strategy.
        MLE_run_hesse: Whether to request a covariance/HESSE step.
        MLE_request_minos: Whether to request MINOS when supported by the build.
        MLE_verbose: Whether to run the minimization backend in verbose mode.
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
        MLE_trace_first_evals: Whether to print the first likelihood evaluations
            for debugging.
        MLE_trace_max_evals: Maximum number of traced evaluations.
        MLE_allow_profile_hessian_fallback: Whether the fitter may reconstruct a
            profile covariance numerically when the backend covariance is
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
            sensitivity pruning.
        nuisance_sensitivity_context_sigma: Width, in nuisance sigma units, used
            to randomize non-central sensitivity contexts.
        nuisance_sensitivity_seed: Random seed for sensitivity contexts.
        nuisance_sensitivity_keep_on_failure: Whether a nuisance is kept when
            its sensitivity evaluation fails.

    Examples:
        >>> cfg = StatisticConfig(
        ...     MC_draws=1000,
        ...     likelihood_mode=StatisticLikelihoodMode.CHI2_MC_COVARIANCE,
        ... )
        >>> cpp_cfg = cfg.to_cpp()
    """

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
            Bound C++ ``StatisticConfig`` object.
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
        cpp.nuisance_sensitivity_keep_on_failure = bool(
            self.nuisance_sensitivity_keep_on_failure
        )

        return cpp


__all__ = ["StatisticConfig", "StatisticLikelihoodMode"]


if __name__ == "__main__":
    cfg = StatisticConfig()
    print(cfg)

    cpp_cfg = cfg.to_cpp()
    print(cpp_cfg)