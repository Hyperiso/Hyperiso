"""Public Python interface to the statistical analysis backend.

``StatisticInterface`` is a high-level wrapper around the C++ statistic manager.
It connects a configured ``ObservableInterface`` to uncertainty propagation,
maximum-likelihood fitting, confidence-contour computation and likelihood scans.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from pyhyperiso.phyperiso.pyhyperiso import statistic as st
from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
from pyhyperiso.core.Common.BinnedObservableId import BinnedObservableId
from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Statistic.ExperimentObs import ExperimentObs
from pyhyperiso.core.Statistic.GaussianSummary import GaussianSummary
from pyhyperiso.core.Statistic.MCResult import MCResult
from pyhyperiso.core.Statistic.StatisticConfig import StatisticConfig


class ProfilingMethod(Enum):
    """Strategy used to reduce a fit with more than two parameters to 2D.

    Attributes:
        SLICE: Vary the selected axes while keeping the other fit parameters at
            their best-fit values.
        FREE_PROJECTION: Profile all non-displayed fit parameters freely.
        PRIOR_CONSTRAINED_PROJECTION: Profile non-displayed fit parameters with
            Gaussian constraints derived from the global fit.
    """

    SLICE = st.ProfilingMethod.SLICE
    FREE_PROJECTION = st.ProfilingMethod.FREE_PROJECTION
    PRIOR_CONSTRAINED_PROJECTION = st.ProfilingMethod.PRIOR_CONSTRAINED_PROJECTION

    def to_cpp(self):
        """Convert to the bound C++ ``ProfilingMethod`` enum."""
        return self.value


class ContourAlgorithm(Enum):
    """Algorithm used to extract the requested confidence contour."""

    AMS = st.ContourAlgorithm.AMS
    MINUIT = st.ContourAlgorithm.MINUIT

    def to_cpp(self):
        """Convert to the bound C++ ``ContourAlgorithm`` enum."""
        return self.value


class ProfilerMode(Enum):
    """Backend used while profiling hidden fit or nuisance parameters."""

    MINUIT = st.ProfileBackend.MINUIT
    LAPLACE_NUISANCE = st.ProfileBackend.LAPLACE_NUISANCE

    def to_cpp(self):
        """Convert to the C++ contour ``ProfileBackend`` enum."""
        return self.value


# Historical public name kept for source compatibility.  ``ProfilerMode`` is
# the terminology used by the contour engine and by the GUI.
ProfileBackend = ProfilerMode


def _require(value, typ, name: str):
    """Validate a typed argument and return it unchanged."""
    if not isinstance(value, typ):
        raise TypeError(f"{name} doit être {typ.__name__}, reçu {type(value)!r}.")
    return value


def _cpp_param_id(pid: ParamId):
    """Convert a Python ``ParamId`` to C++."""
    return _require(pid, ParamId, "ParamId").to_cpp()


def _param_from_cpp(cpp_obj) -> ParamId:
    """Convert a C++ ``ParamId`` to the Python wrapper."""
    return ParamId.from_cpp(cpp_obj)


def _param_float_map_from_cpp(cpp_map) -> Dict[ParamId, float]:
    """Convert a C++ ``map<ParamId, double>`` to a Python dictionary."""
    return {_param_from_cpp(k): float(v) for k, v in dict(cpp_map).items()}


def _param_float_map_to_cpp(values: Mapping[ParamId, float]):
    """Convert a Python parameter-value mapping to C++."""
    return {_cpp_param_id(k): float(v) for k, v in values.items()}


def _param_nested_float_map_from_cpp(cpp_map) -> Dict[ParamId, Dict[ParamId, float]]:
    """Convert a nested C++ parameter matrix map to Python dictionaries."""
    return {_param_from_cpp(k): _param_float_map_from_cpp(v) for k, v in dict(cpp_map).items()}


def _experiment_float_map_from_cpp(cpp_map) -> Dict[ExperimentObs, float]:
    """Convert a C++ ``map<ExperimentObs, double>`` to Python."""
    return {ExperimentObs.from_cpp(k): float(v) for k, v in dict(cpp_map).items()}


def _experiment_nested_float_map_from_cpp(cpp_map) -> Dict[ExperimentObs, Dict[ExperimentObs, float]]:
    """Convert a nested experimental-observable correlation map to Python."""
    return {ExperimentObs.from_cpp(k): _experiment_float_map_from_cpp(v) for k, v in dict(cpp_map).items()}


def _cpp_profiling_method(value: ProfilingMethod):
    """Convert a Python profiling method to C++."""
    return _require(value, ProfilingMethod, "profiling_method").to_cpp()


def _cpp_profile_backend(value: ProfilerMode):
    """Convert a Python profile backend to C++."""
    return _require(value, ProfilerMode, "profile_backend").to_cpp()


def _cpp_contour_algorithm(value: ContourAlgorithm):
    """Convert a Python contour algorithm to C++."""
    return _require(value, ContourAlgorithm, "contour_algorithm").to_cpp()


@dataclass(frozen=True)
class FitResultWithMaps:
    """Maximum-likelihood fit result keyed by Python parameter identifiers.

    Attributes:
        fit_ok: Whether the C++ fit result is considered usable.
        p_hat: Best-fit values for fit parameters.
        eta_hat: Profiled nuisance values at the best-fit point.
        p_hat_std: Standard deviations for fit parameters, usually from the
            profiled covariance.
        p_correlations: Fit-parameter correlation matrix as nested maps.
        ell_hat: Minimum negative log-likelihood value.
    """

    fit_ok: bool = False
    p_hat: Dict[ParamId, float] = field(default_factory=dict)
    eta_hat: Dict[ParamId, float] = field(default_factory=dict)
    p_hat_std: Dict[ParamId, float] = field(default_factory=dict)
    p_correlations: Dict[ParamId, Dict[ParamId, float]] = field(default_factory=dict)
    ell_hat: float = 0.0

    @classmethod
    def from_cpp(cls, cpp_obj) -> "FitResultWithMaps":
        """Create a Python fit result from the bound C++ result."""
        return cls(
            fit_ok=bool(cpp_obj.fit_ok),
            p_hat=_param_float_map_from_cpp(cpp_obj.p_hat),
            eta_hat=_param_float_map_from_cpp(cpp_obj.eta_hat),
            p_hat_std=_param_float_map_from_cpp(cpp_obj.p_hat_std),
            p_correlations=_param_nested_float_map_from_cpp(cpp_obj.p_correlations),
            ell_hat=float(cpp_obj.ell_hat),
        )

    def to_cpp(self):
        """Convert this result to the bound C++ representation."""
        cpp = st.FitResultWithMaps()
        cpp.fit_ok = bool(self.fit_ok)
        cpp.p_hat = _param_float_map_to_cpp(self.p_hat)
        cpp.eta_hat = _param_float_map_to_cpp(self.eta_hat)
        cpp.p_hat_std = _param_float_map_to_cpp(self.p_hat_std)
        cpp.p_correlations = {_cpp_param_id(k): _param_float_map_to_cpp(v) for k, v in self.p_correlations.items()}
        cpp.ell_hat = float(self.ell_hat)
        return cpp


@dataclass(frozen=True)
class MLFitOptions:
    """Low-level maximum-likelihood fitter options.

    This mirrors the C++ ``MLFitOptions`` struct. In most user workflows these
    settings are configured through :class:`StatisticConfig`; this class is kept
    available for direct conversions and advanced integrations.

    Attributes:
        run_hesse: Ask the minimizer to compute a covariance/HESSE estimate.
        request_minos: Request MINOS if the compiled backend supports it.
        verbose: Enable verbose minimizer output.
        strategy: Backend strategy; ``0`` means backend default.
        max_fcn: Maximum number of function calls; ``0`` means backend default.
        tolerance: Minimizer tolerance; ``0.0`` means backend default.
        allow_profile_hessian_fallback: Allow numerical profile-Hessian fallback.
        profile_hessian_step_scale: Step scale for the fallback Hessian.
        profile_hessian_eig_floor_rel: Relative eigenvalue floor for Hessian
            regularization.
        trace_first_evals: Print the first likelihood evaluations.
        trace_max_evals: Maximum number of traced evaluations.
    """

    run_hesse: bool = True
    request_minos: bool = False
    verbose: bool = False
    strategy: int = 0
    max_fcn: int = 0
    tolerance: float = 0.0
    allow_profile_hessian_fallback: bool = True
    profile_hessian_step_scale: float = 1.0
    profile_hessian_eig_floor_rel: float = 1e-8
    trace_first_evals: bool = False
    trace_max_evals: int = 25

    @classmethod
    def from_cpp(cls, cpp_obj) -> "MLFitOptions":
        """Create Python fitter options from C++ options."""
        return cls(
            run_hesse=bool(cpp_obj.run_hesse),
            request_minos=bool(cpp_obj.request_minos),
            verbose=bool(cpp_obj.verbose),
            strategy=int(cpp_obj.strategy),
            max_fcn=int(cpp_obj.max_fcn),
            tolerance=float(cpp_obj.tolerance),
            allow_profile_hessian_fallback=bool(cpp_obj.allow_profile_hessian_fallback),
            profile_hessian_step_scale=float(cpp_obj.profile_hessian_step_scale),
            profile_hessian_eig_floor_rel=float(cpp_obj.profile_hessian_eig_floor_rel),
            trace_first_evals=bool(cpp_obj.trace_first_evals),
            trace_max_evals=int(cpp_obj.trace_max_evals),
        )

    def to_cpp(self):
        """Convert this options object to C++."""
        cpp = st.MLFitOptions()
        cpp.run_hesse = bool(self.run_hesse)
        cpp.request_minos = bool(self.request_minos)
        cpp.verbose = bool(self.verbose)
        cpp.strategy = int(self.strategy)
        cpp.max_fcn = int(self.max_fcn)
        cpp.tolerance = float(self.tolerance)
        cpp.allow_profile_hessian_fallback = bool(self.allow_profile_hessian_fallback)
        cpp.profile_hessian_step_scale = float(self.profile_hessian_step_scale)
        cpp.profile_hessian_eig_floor_rel = float(self.profile_hessian_eig_floor_rel)
        cpp.trace_first_evals = bool(self.trace_first_evals)
        cpp.trace_max_evals = int(self.trace_max_evals)
        return cpp


@dataclass(frozen=True)
class ContourOptions:
    """Options controlling confidence-contour computation.

    Attributes:
        profiling_method: How the two-dimensional contour fixes or profiles
            parameters.
        profile_backend: Backend used to evaluate profiled NLL points.
        primary_contour_method: Contour extraction algorithm.
        fallback_contour_method: Optional fallback algorithm used by the C++
            engine if the primary extraction fails.
        resolution: Resolution parameter passed to the contour extractor.

    Examples:
        >>> opts = ContourOptions(profile_backend=ProfilerMode.LAPLACE_NUISANCE, resolution=60)
    """

    profiling_method: ProfilingMethod = ProfilingMethod.SLICE
    profile_backend: ProfilerMode = ProfilerMode.LAPLACE_NUISANCE
    primary_contour_method: ContourAlgorithm = ContourAlgorithm.MINUIT
    fallback_contour_method: Optional[ContourAlgorithm] = None
    resolution: int = 40

    @classmethod
    def from_cpp(cls, cpp_obj) -> "ContourOptions":
        """Create Python contour options from C++ options."""
        return cls(
            profiling_method=ProfilingMethod(cpp_obj.profiling_method),
            profile_backend=ProfilerMode(cpp_obj.profile_backend),
            primary_contour_method=ContourAlgorithm(cpp_obj.primary_contour_method),
            fallback_contour_method=None if cpp_obj.fallback_contour_method is None else ContourAlgorithm(cpp_obj.fallback_contour_method),
            resolution=int(cpp_obj.resolution),
        )

    def to_cpp(self):
        """Convert this contour options object to C++."""
        cpp = st.ContourOptions()
        cpp.profiling_method = _cpp_profiling_method(self.profiling_method)
        cpp.profile_backend = _cpp_profile_backend(self.profile_backend)
        cpp.primary_contour_method = _cpp_contour_algorithm(self.primary_contour_method)
        cpp.fallback_contour_method = None if self.fallback_contour_method is None else _cpp_contour_algorithm(self.fallback_contour_method)
        cpp.resolution = int(self.resolution)
        return cpp


class Contour:
    """Python view of a C++ confidence contour.

    The bound object exposes one or more disconnected paths.  Each path is a
    sequence of ``(x, y)`` points in the plane of the two displayed fit
    parameters.
    """

    __slots__ = ("_cpp_obj",)

    def __init__(self, cpp_obj) -> None:
        self._cpp_obj = cpp_obj

    @classmethod
    def from_cpp(cls, cpp_obj) -> "Contour":
        return cls(cpp_obj)

    @property
    def paths(self) -> List[List[tuple[float, float]]]:
        return [
            [(float(point[0]), float(point[1])) for point in path]
            for path in self._cpp_obj.paths
        ]

    @property
    def level(self) -> float:
        return float(self._cpp_obj.level)

    @property
    def success(self) -> bool:
        return bool(self._cpp_obj.success)

    def _to_cpp(self):
        """Return the underlying C++ contour object."""
        return self._cpp_obj

    def __repr__(self) -> str:
        return f"Contour(success={self.success}, level={self.level:.6g}, paths={len(self.paths)})"


@dataclass(frozen=True)
class LikelihoodScanPoint:
    """One grid point of a two-dimensional likelihood scan.

    Attributes:
        x: Value of the first scanned parameter.
        y: Value of the second scanned parameter.
        nll: Absolute negative log-likelihood value.
        delta_nll: Difference with respect to the scan reference point.
    """

    x: float = 0.0
    y: float = 0.0
    nll: float = 0.0
    delta_nll: float = 0.0

    @classmethod
    def from_cpp(cls, cpp_obj) -> "LikelihoodScanPoint":
        """Create a Python scan point from a bound C++ scan point."""
        return cls(x=float(cpp_obj.x), y=float(cpp_obj.y), nll=float(cpp_obj.nll), delta_nll=float(cpp_obj.delta_nll))

    def to_cpp(self):
        """Convert this scan point to C++."""
        cpp = st.LikelihoodScanPoint()
        cpp.x = float(self.x)
        cpp.y = float(self.y)
        cpp.nll = float(self.nll)
        cpp.delta_nll = float(self.delta_nll)
        return cpp


@dataclass(frozen=True)
class LikelihoodScanGrid:
    """Rectangular two-dimensional likelihood scan result.

    Attributes:
        x_param: First scanned fit parameter.
        y_param: Second scanned fit parameter.
        x_center: Reference value used to center the scan in x.
        y_center: Reference value used to center the scan in y.
        nx: Number of x grid points.
        ny: Number of y grid points.
        points: Flattened list of scan points.
    """

    x_param: ParamId
    y_param: ParamId
    x_center: float = 0.0
    y_center: float = 0.0
    nx: int = 0
    ny: int = 0
    points: List[LikelihoodScanPoint] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate scan parameter identifiers."""
        _require(self.x_param, ParamId, "x_param")
        _require(self.y_param, ParamId, "y_param")

    @classmethod
    def from_cpp(cls, cpp_obj) -> "LikelihoodScanGrid":
        """Create a Python scan grid from the bound C++ result."""
        return cls(
            x_param=_param_from_cpp(cpp_obj.x_param),
            y_param=_param_from_cpp(cpp_obj.y_param),
            x_center=float(cpp_obj.x_center),
            y_center=float(cpp_obj.y_center),
            nx=int(cpp_obj.nx),
            ny=int(cpp_obj.ny),
            points=[LikelihoodScanPoint.from_cpp(p) for p in cpp_obj.points],
        )

    def to_cpp(self):
        """Convert this scan grid to the bound C++ representation."""
        cpp = st.LikelihoodScanGrid()
        cpp.x_param = _cpp_param_id(self.x_param)
        cpp.y_param = _cpp_param_id(self.y_param)
        cpp.x_center = float(self.x_center)
        cpp.y_center = float(self.y_center)
        cpp.nx = int(self.nx)
        cpp.ny = int(self.ny)
        cpp.points = [p.to_cpp() for p in self.points]
        return cpp


class StatisticInterface:
    """High-level Python facade for statistical fits and scans.

    Args:
        config: Statistic configuration. Python-only fields such as ``p_specs``
            and ``selected_experiments`` are used by this wrapper.
        observable_interface: Observable interface already configured with the
            observables to compute.

    Examples:
        A typical workflow is::

            stat_cfg = StatisticConfig(MC_draws=1000)
            stat_cfg.p_specs = [my_param_id]
            stat = StatisticInterface(stat_cfg, observable_interface=obs_int)

            fit = stat.compute_MLE()
            print(fit.fit_ok, fit.p_hat)

            contour = stat.compute_confidence_contour(
                p1=my_param_id,
                p2=other_param_id,
                z=1.0,
                bounds=[-1.0, 1.0, -1.0, 1.0],
            )
    """

    def __init__(self, config: StatisticConfig, observable_interface: ObservableInterface) -> None:
        """Create the C++ statistic interface and apply initial selections."""
        _require(config, StatisticConfig, "config")
        _require(observable_interface, ObservableInterface, "observable_interface")
        self.config = config
        self._cpp = st.StatisticInterface(config.to_cpp(), observable_interface._to_cpp())
        # if config.selected_experiments is not None:
        #     self.select_experiments(config.selected_experiments)

    def _to_cpp(self):
        """Return the underlying C++ statistic interface."""
        return self._cpp

    def select_experiment(self, experiment: str) -> None:
        """Restrict the statistic manager to one experiment.

        Args:
            experiment: Experiment name as stored in the experimental database.
        """
        self._cpp.select_experiment(str(experiment))

    def select_experiments(self, experiments: Sequence[str]) -> None:
        """Restrict the statistic manager to a set of experiments.

        Args:
            experiments: Experiment names to keep.
        """
        self._cpp.select_experiments([str(x) for x in experiments])

    def select_experiments_all(self) -> None:
        """Clear any experiment selection and use all available experiments."""
        self._cpp.select_experiments_all()

    def has_experiment_selection(self) -> bool:
        """Return whether an experiment filter is currently active."""
        return bool(self._cpp.has_experiment_selection())

    def selected_experiments(self) -> List[str]:
        """Return the currently selected experiment names."""
        return [str(x) for x in self._cpp.selected_experiments()]

    def reload_nuisance_specs(self) -> None:
        """Reload default and user nuisance specifications in the C++ manager."""
        self._cpp.reload_nuisance_specs()

    def set_nuisance_user_file(self, user_yaml_path: str | Path) -> None:
        """Use a custom user nuisance-specification file.

        Args:
            user_yaml_path: Path to the YAML file overriding default nuisance
                specifications.
        """
        self._cpp.set_nuisance_user_file(str(user_yaml_path))

    def clear_nuisance_user_file(self) -> None:
        """Return to the default user nuisance-specification path."""
        self._cpp.clear_nuisance_user_file()

    def update_cache(self, p_specs: Optional[Sequence[ParamId]] = None) -> None:
        """Refresh the C++ statistic cache.

        Args:
            p_specs: Fit parameters to detach from the nuisance set. When
                omitted, ``self.config.p_specs`` is used.
        """
        p_specs = self.config.p_specs if p_specs is None else p_specs
        self._cpp.update_cache([_cpp_param_id(p) for p in p_specs])

    def compute_uncertainties(self) -> Dict[BinnedObservableId, GaussianSummary]:
        """Propagate nuisance uncertainties and return Gaussian summaries.

        Returns:
            Mapping from observable id to Gaussian or split-Gaussian summary.
        """
        return {BinnedObservableId.from_cpp(k): GaussianSummary.from_cpp(v) for k, v in self._cpp.compute_uncertainties().items()}

    def compute_uncertainties_and_sampling(self) -> MCResult:
        """Run Monte-Carlo uncertainty propagation and keep raw samples.

        Returns:
            ``MCResult`` containing raw samples, summaries and covariance.
        """
        return MCResult.from_cpp(self._cpp.compute_uncertainties_and_sampling())

    def compute_MLE(self, p_specs: Optional[Sequence[ParamId]] = None) -> FitResultWithMaps:
        """Compute the maximum-likelihood estimate for selected fit parameters.

        Args:
            p_specs: Fit parameters. When omitted, ``self.config.p_specs`` is
                used.

        Returns:
            Fit result containing best-fit values, profiled nuisances,
            uncertainties, correlations and minimum NLL.

        Raises:
            Exception: Propagates C++ errors raised by cache construction,
                likelihood evaluation or minimization.
        """
        p_specs = self.config.p_specs if p_specs is None else p_specs
        return FitResultWithMaps.from_cpp(self._cpp.compute_MLE([_cpp_param_id(p) for p in p_specs]))

    def compute_confidence_contour(
        self,
        p1: ParamId,
        p2: ParamId,
        z: float,
        bounds: Sequence[float],
        options: Optional[ContourOptions] = None,
    ) -> Contour:
        """Compute a two-dimensional confidence contour.

        Args:
            p1: First fit parameter.
            p2: Second fit parameter.
            z: Gaussian-equivalent confidence level, for example ``1.0`` for an
                approximate one-sigma contour.
            bounds: Four values ``[xmin, xmax, ymin, ymax]`` defining the search
                rectangle.
            options: Optional contour/profiling options. Defaults to
                ``ContourOptions()``.

        Returns:
            ``Contour`` wrapper exposing ``paths``, ``level`` and ``success``.

        Raises:
            ValueError: If ``bounds`` does not contain exactly four values.
            Exception: Propagates backend failures, for example if MLE has not
                been computed before requesting the contour.
        """
        if len(bounds) != 4:
            raise ValueError("bounds doit contenir exactement 4 valeurs : [xmin, xmax, ymin, ymax].")
        cpp_options = ContourOptions().to_cpp() if options is None else _require(options, ContourOptions, "options").to_cpp()
        return Contour.from_cpp(
            self._cpp.compute_confidence_contour(
                _cpp_param_id(p1),
                _cpp_param_id(p2),
                float(z),
                [float(x) for x in bounds],
                cpp_options,
            )
        )

    def prepare_likelihood_for_scan(self, p_specs: Optional[Sequence[ParamId]] = None) -> None:
        """Prepare and cache a likelihood object for manual two-dimensional scans.

        Args:
            p_specs: Fit parameters used by the likelihood. Defaults to
                ``self.config.p_specs``.
        """
        p_specs = self.config.p_specs if p_specs is None else p_specs
        self._cpp.prepare_likelihood_for_scan([_cpp_param_id(p) for p in p_specs])

    def set_manual_scan_point(self, p_hat: Mapping[ParamId, float], eta_hat: Mapping[ParamId, float]) -> None:
        """Override the scan reference point manually.

        Args:
            p_hat: Fit-parameter values used as scan center/reference.
            eta_hat: Nuisance values used as scan reference.
        """
        self._cpp.set_manual_scan_point(_param_float_map_to_cpp(p_hat), _param_float_map_to_cpp(eta_hat))

    def scan_likelihood_around_current_point(
        self,
        p1: ParamId,
        p2: ParamId,
        x_half_width: float,
        y_half_width: float,
        nx: int,
        ny: int,
    ) -> LikelihoodScanGrid:
        """Evaluate the likelihood on a rectangular grid around the current point.

        Args:
            p1: First scanned parameter.
            p2: Second scanned parameter.
            x_half_width: Half-width of the scan range in the first parameter.
            y_half_width: Half-width of the scan range in the second parameter.
            nx: Number of grid points in the x direction.
            ny: Number of grid points in the y direction.

        Returns:
            A ``LikelihoodScanGrid`` containing all evaluated points.
        """
        return LikelihoodScanGrid.from_cpp(
            self._cpp.scan_likelihood_around_current_point(
                _cpp_param_id(p1),
                _cpp_param_id(p2),
                float(x_half_width),
                float(y_half_width),
                int(nx),
                int(ny),
            )
        )

    def save_likelihood_scan_csv(self, path: str | Path, grid: LikelihoodScanGrid) -> None:
        """Save a likelihood scan grid to a CSV file.

        Args:
            path: Output CSV path.
            grid: Scan grid returned by :meth:`scan_likelihood_around_current_point`.
        """
        self._cpp.save_likelihood_scan_csv(str(path), _require(grid, LikelihoodScanGrid, "grid").to_cpp())

    def get_all_obss_deps(self) -> Dict[ParamId, float]:
        """Return all selected observable dependencies after nuisance pruning."""
        return _param_float_map_from_cpp(self._cpp.get_all_obss_deps())

    def get_active_observable_dependencies(self) -> Dict[ParamId, float]:
        """Return dependencies currently visible to the Statistic manager.

        This method is useful for runtime/lambda observables: it verifies that
        dependencies declared on ``LambdaObservableConfig`` and propagated from
        custom Wilson lambdas reached the statistic layer.
        """
        return _param_float_map_from_cpp(self._cpp.get_active_observable_dependencies())

    def get_p_specs(self, p_specs: Optional[Sequence[ParamId]] = None) -> Dict[ParamId, float]:
        """Return central values for selected fit parameters.

        Args:
            p_specs: Fit parameters to query. Defaults to ``self.config.p_specs``.
        """
        p_specs = self.config.p_specs if p_specs is None else p_specs
        return _param_float_map_from_cpp(self._cpp.get_p_specs([_cpp_param_id(p) for p in p_specs]))

    def get_all_correlations(self) -> Dict[ParamId, Dict[ParamId, float]]:
        """Return the nuisance correlation matrix as nested parameter maps."""
        return _param_nested_float_map_from_cpp(self._cpp.get_all_correlations())

    def get_all_obs_correlations(self) -> Dict[ExperimentObs, Dict[ExperimentObs, float]]:
        """Return experimental-observable correlations as nested maps."""
        return _experiment_nested_float_map_from_cpp(self._cpp.get_all_obs_correlations())

    def get_obs_exp(self) -> Dict[ExperimentObs, float]:
        """Return experimental central values used by the statistic manager."""
        return _experiment_float_map_from_cpp(self._cpp.get_obs_exp())

    def print_cache(self) -> None:
        """Print the current C++ statistic cache for debugging."""
        self._cpp.print_cache()


__all__ = [
    "StatisticInterface",
    "FitResultWithMaps",
    "Contour",
    "ContourOptions",
    "ProfilingMethod",
    "ContourAlgorithm",
    "ProfilerMode",
    "ProfileBackend",
    "MLFitOptions",
    "LikelihoodScanPoint",
    "LikelihoodScanGrid",
]


if __name__ == "__main__":
    
    from pyhyperiso.core.Core.HyperisoMaster import HyperisoMaster
    from pyhyperiso.core.Core.HyperisoConfig import HyperisoConfig, ExternalFlag
    from pyhyperiso.core.Common.GeneralEnum import Model, Observables, QCDOrder, ParameterType
    from pyhyperiso.core.Common.ParamId import ParamId
    from pathlib import Path
    config = HyperisoConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: False,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
            ExternalFlag.HYP_AS_SM_MARTY: True
        },
        model=Model.SM,
        mty_model_name="MSSM_UFO",
        mty_model_path=Path("/my/custom/marty/path")
    )
    
    hyp = HyperisoMaster()
    lha_file_path = "/home/theo/hyperiso/Assets/lha/si_input.flha" 
    
    hyp.init(lha_file=lha_file_path, config=config)
    
    obs = {Observables.BR_BS_MUMU : QCDOrder.LO, Observables.BR_BS_MUMU_UNTAG : QCDOrder.LO, Observables.BR_BD_MUMU : QCDOrder.LO}

    from pyhyperiso.core.BusinessLogic.ObservableInterface import ObservableInterface
    
    obs_int = ObservableInterface()
    
    obs_int.add_observables(obs, True)
    
    config_stat = StatisticConfig()
    config_stat.p_specs = [ParamId(ParameterType.DECAY, "B_ll", 1)]

    si = StatisticInterface(config_stat, observable_interface=obs_int)
    fit = si.compute_MLE()

    print(fit.fit_ok)
    print(fit.p_hat)
    print(fit.p_hat_std)