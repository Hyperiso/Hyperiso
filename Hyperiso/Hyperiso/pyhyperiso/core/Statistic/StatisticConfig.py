"""High-level configuration object for the statistic workflow.

``StatisticConfig`` is the Python entry point for configuring the bound C++
statistic manager. It mirrors the fields exposed by the pybind11
``StatisticConfig`` binding.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Optional, Tuple

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.Copula import CopulaKind
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.MarginalConfig import MarginalKind


@dataclass(frozen=True)
class StatisticProgressSnapshot:
    """Immutable progress snapshot produced by the C++ statistic pipeline."""

    phase: str = "idle"
    message: str = ""
    completed: int = 0
    total: int = 0
    attempts: int = 0
    failures: int = 0
    fraction: float = 0.0
    elapsed_seconds: float = 0.0
    eta_seconds: float = -1.0
    finished: bool = False
    sequence: int = 0


class StatisticProgressMonitor:
    """Thread-safe bridge exposing C++ progress to Python frontends."""

    def __init__(self):
        self._cpp_obj = st.StatisticProgressMonitor()

    def reset(
        self, phase: str = "preparing", message: str = "Preparing statistic workflow"
    ) -> None:
        """Reset progress reporting and set the initial phase and message."""
        self._cpp_obj.reset(str(phase), str(message))

    def set_progress(
        self,
        phase: str,
        message: str,
        fraction: float,
        *,
        completed: int = 0,
        total: int = 0,
        eta_seconds: float = -1.0,
        finished: bool = False,
    ) -> None:
        """Publish a progress update to Python and graphical frontends."""
        self._cpp_obj.set_progress(
            str(phase),
            str(message),
            float(fraction),
            int(completed),
            int(total),
            float(eta_seconds),
            bool(finished),
        )

    def snapshot(self) -> StatisticProgressSnapshot:
        """Return an immutable snapshot of the latest progress event."""
        event = self._cpp_obj.snapshot()
        return StatisticProgressSnapshot(
            phase=str(event.phase),
            message=str(event.message),
            completed=int(event.completed),
            total=int(event.total),
            attempts=int(event.attempts),
            failures=int(event.failures),
            fraction=float(event.fraction),
            elapsed_seconds=float(event.elapsed_seconds),
            eta_seconds=float(event.eta_seconds),
            finished=bool(event.finished),
            sequence=int(event.sequence),
        )


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
        raise TypeError(f"{name} must be {typ.__name__}, got {type(value).__name__}.")
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
class AdvancedStatisticConfig:
    """Expert configuration for fitting, nuisance pruning and covariance logic.

    ``StatisticConfig`` should remain the simple user-facing object.  Fields that
    affect minimization internals, likelihood backend selection, covariance
    regularization or nuisance-preselection live here.
    """

    override_nuisance_marginals: Dict[ParamId, MarginalKind] = field(default_factory=dict)
    override_exp_data_marginals: Dict[ExperimentObs, MarginalKind] = field(default_factory=dict)

    nuisance_copula_type: CopulaKind = CopulaKind.GAUSSIAN
    exp_data_copula_type: CopulaKind = CopulaKind.GAUSSIAN

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
    nuisance_sensitivity_contexts: int = 2
    nuisance_sensitivity_context_sigma: float = 0.35
    nuisance_sensitivity_seed: int = 12345
    nuisance_sensitivity_keep_on_failure: bool = True

    MLE_trace_first_evals: bool = False
    MLE_trace_max_evals: int = 25
    MLE_allow_profile_hessian_fallback: bool = True
    MLE_profile_hessian_step_scale: float = 1.0
    MLE_profile_hessian_eig_floor_rel: float = 1e-8

    likelihood_mode: StatisticLikelihoodMode = StatisticLikelihoodMode.CHI2_MC_COVARIANCE
    chi2_covariance_ridge_rel: float = 1e-8
    chi2_covariance_ridge_abs: float = 1e-12

    # Advanced MC/decay thread arbitration.  MC_threads itself stays on
    # StatisticConfig because it is the main user-facing parallelism knob.
    MC_force_decay_threads_to_one: bool = True
    MC_forced_decay_threads: int = 1

    def to_cpp(self):
        """Convert this Python config to the bound C++ ``AdvancedStatisticConfig``."""
        cpp = st.AdvancedStatisticConfig()
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
        cpp.nuisance_sensitivity_contexts = int(self.nuisance_sensitivity_contexts)
        cpp.nuisance_sensitivity_context_sigma = float(self.nuisance_sensitivity_context_sigma)
        cpp.nuisance_sensitivity_seed = int(self.nuisance_sensitivity_seed)
        cpp.nuisance_sensitivity_keep_on_failure = bool(self.nuisance_sensitivity_keep_on_failure)
        cpp.MLE_trace_first_evals = bool(self.MLE_trace_first_evals)
        cpp.MLE_trace_max_evals = int(self.MLE_trace_max_evals)
        cpp.MLE_allow_profile_hessian_fallback = bool(self.MLE_allow_profile_hessian_fallback)
        cpp.MLE_profile_hessian_step_scale = float(self.MLE_profile_hessian_step_scale)
        cpp.MLE_profile_hessian_eig_floor_rel = float(self.MLE_profile_hessian_eig_floor_rel)
        cpp.likelihood_mode = _cpp_likelihood_mode(self.likelihood_mode)
        cpp.chi2_covariance_ridge_rel = float(self.chi2_covariance_ridge_rel)
        cpp.chi2_covariance_ridge_abs = float(self.chi2_covariance_ridge_abs)
        cpp.MC_force_decay_threads_to_one = bool(self.MC_force_decay_threads_to_one)
        cpp.MC_forced_decay_threads = int(self.MC_forced_decay_threads)
        return cpp


@dataclass
class StatisticConfig:
    """Basic configuration of the statistic pipeline.

    This is the small, common config used by examples and scripts.  Advanced
    minimization, nuisance-pruning and covariance options are grouped in
    :class:`AdvancedStatisticConfig` under ``advanced``.

    Args:
        MC_draws: Number of accepted MC predictions used in uncertainty
            propagation or chi-square covariance estimation.
        MC_threads: Number of worker threads used by MC propagation.
        MC_seed: RNG seed used for reproducible MC nuisance and experimental-data sampling.
        skew_abs_threshold: Skewness threshold for symmetric vs split-Gaussian
            summaries.
        print_mc_progress: Print a compact progress bar with live ETA during MC
            sampling.
        print_chi2_pipeline_progress: In chi-square MC-covariance mode, print
            the post-MC workflow stages so the MC ETA is not mistaken for the
            full runtime.
        print_mc_config: Print nuisance candidates and retained MC marginal
            configuration.
        print_fit_summary: Print high-level fit diagnostics.
        print_scan_summary: Print likelihood-scan diagnostics.
        print_cache_summary: Enable ``StatisticManager::print_cache`` output.
        print_debug: Master debug flag enabling extra diagnostic output.
        write_mc_samples_csv: Persist accepted MC observable samples.
        mc_samples_csv_path: Output CSV path when sample writing is enabled.
        advanced: Expert configuration object.
    """

    MC_draws: int = 100
    MC_threads: int = 1
    MC_seed: int = 123456
    skew_abs_threshold: float = 0.2

    print_mc_progress: bool = False
    print_chi2_pipeline_progress: bool = False
    print_mc_config: bool = False
    print_fit_summary: bool = False
    print_scan_summary: bool = False
    print_cache_summary: bool = False
    print_debug: bool = False

    write_mc_samples_csv: bool = False
    mc_samples_csv_path: str = "obs_samples.csv"
    mc_progress_probe_draws: int = 5
    mc_progress_update_every: int = 1

    progress_monitor: Optional[StatisticProgressMonitor] = None
    fit_parameter_bounds: Dict[ParamId, Tuple[float, float]] = field(default_factory=dict)
    fit_parameter_offsets: Dict[ParamId, float] = field(default_factory=dict)
    advanced: AdvancedStatisticConfig = field(default_factory=AdvancedStatisticConfig)

    def to_cpp(self):
        """Convert this Python config to the bound C++ ``StatisticConfig``."""
        cpp = st.StatisticConfig()
        cpp.MC_draws = int(self.MC_draws)
        cpp.MC_threads = int(self.MC_threads)
        cpp.MC_seed = int(self.MC_seed)
        cpp.skew_abs_threshold = float(self.skew_abs_threshold)
        cpp.print_mc_progress = bool(self.print_mc_progress)
        cpp.print_chi2_pipeline_progress = bool(self.print_chi2_pipeline_progress)
        cpp.print_mc_config = bool(self.print_mc_config)
        cpp.print_fit_summary = bool(self.print_fit_summary)
        cpp.print_scan_summary = bool(self.print_scan_summary)
        cpp.print_cache_summary = bool(self.print_cache_summary)
        cpp.print_debug = bool(self.print_debug)
        cpp.write_mc_samples_csv = bool(self.write_mc_samples_csv)
        cpp.mc_samples_csv_path = str(self.mc_samples_csv_path)
        cpp.mc_progress_probe_draws = int(self.mc_progress_probe_draws)
        cpp.mc_progress_update_every = int(self.mc_progress_update_every)
        if self.progress_monitor is not None:
            cpp.progress_monitor = self.progress_monitor._cpp_obj
        cpp.fit_parameter_bounds = {
            _cpp_param_id(pid): (float(bounds[0]), float(bounds[1]))
            for pid, bounds in self.fit_parameter_bounds.items()
        }
        cpp.fit_parameter_offsets = {
            _cpp_param_id(pid): float(offset) for pid, offset in self.fit_parameter_offsets.items()
        }
        cpp.advanced = self.advanced.to_cpp()
        return cpp


__all__ = [
    "StatisticConfig",
    "AdvancedStatisticConfig",
    "StatisticLikelihoodMode",
    "StatisticProgressMonitor",
    "StatisticProgressSnapshot",
]
