# Statistical workflows {#statistics}

The statistical layer consumes the observables registered in
@ref ObservableInterface. It can propagate nuisance uncertainties, build a
likelihood, fit one or more physics parameters and produce likelihood scans or
confidence contours.

## Default likelihood

The default `StatisticConfig` uses
`StatisticLikelihoodMode::CHI2_MC_COVARIANCE`. Nuisance parameters are sampled
through the observable calculator to estimate a theory covariance matrix,
which is added to the experimental covariance:

\f[
  \Sigma = \Sigma_{\mathrm{exp}} + \Sigma_{\mathrm{th}}^{\mathrm{MC}}.
\f]

For residual vector \f$r(p)=f(p)-O_{\mathrm{exp}}\f$, the minimized objective is
proportional to

\f[
  \chi^2(p) = r(p)^T\,\Sigma^{-1}\,r(p).
\f]

The alternative `PROFILED_NUISANCE` mode retains explicit nuisance coordinates
and profiles them during the fit. It is more flexible but can be substantially
more expensive.

## Main configuration

Frequently adjusted options are intentionally kept at the top level:

| Option | Default | Purpose |
|---|---:|---|
| `MC_draws` | 100 | Accepted Monte-Carlo draws. |
| `MC_threads` | 1 | Parallel Monte-Carlo workers. |
| `MC_seed` | 123456 | Reproducible random seed. |
| `print_mc_progress` | true | Progress and ETA reporting. |
| `fit_parameter_bounds` | empty | Optional explicit minimizer bounds. |

Expert fit, covariance, pruning and likelihood settings are grouped under
`StatisticConfig::advanced` / `StatisticConfig.advanced`.

## Python example

```python
from pyhyperiso.Common import Decays, ParamId, ParameterType, QCDOrder
from pyhyperiso.Observable import ObservableInterface
from pyhyperiso.Statistic import (
    ContourOptions,
    StatisticConfig,
    StatisticInterface,
    StatisticLikelihoodMode,
)

obs = ObservableInterface()
obs.add_observables_from_decay(Decays.B__l_l, QCDOrder.NNLO)

config = StatisticConfig()
config.MC_draws = 500
config.MC_threads = 4
config.advanced.likelihood_mode = (
    StatisticLikelihoodMode.CHI2_MC_COVARIANCE
)

stats = StatisticInterface(config, obs)

f_bd = ParamId(ParameterType.FLAVOR, "FCONST", [511, 1])
f_bs = ParamId(ParameterType.FLAVOR, "FCONST", [531, 1])

uncertainties = stats.compute_uncertainties()
fit = stats.compute_MLE([f_bd, f_bs])
contour = stats.compute_confidence_contour(
    f_bd,
    f_bs,
    1.0,
    [-0.5, 0.5, -0.5, 0.5],
    ContourOptions(),
)
```

## Fit-parameter sensitivity check

Before fitting, HyperIso probes each requested physics parameter (`p_spec`) and
recomputes the selected observables around the central point. If no selected
observable changes beyond the configured absolute and relative thresholds, the
fit is rejected.

This protects users from:

- unconstrained flat directions;
- arbitrary best-fit values returned by a minimizer;
- singular or misleading covariance matrices;
- meaningless confidence contours.

The check is enabled by default through
`advanced.fit_parameter_sensitivity_check`. Its probe fraction and numerical
cutoffs are expert settings. A failure to execute the probe can be configured
to retain the parameter with a warning, but a confirmed absence of sensitivity
is treated as an invalid fit request.

## Nuisance selection

Nuisance dependencies are collected from the selected observable graph and can
be pruned by numerical relevance. Correlations and marginal distributions are
then assembled consistently for the retained nuisance set.

The physics parameters selected for fitting are detached from nuisance
coordinates so the same quantity is not varied twice.

## Performance guidance

- Increase `MC_threads` only when the machine has enough physical cores and
  memory bandwidth.
- Start contours with a modest grid resolution, then increase it for final
  plots.
- The uncertainty calculation and chi-square covariance construction may run
  separate Monte-Carlo passes; omit calculations that are not needed for a
  quick fit test.
- Use a fixed seed and one thread for release-level numerical regression.
- Enable progress and fit summaries when diagnosing a calculation that appears
  idle.

## Interpretation

A successful numerical minimization does not by itself establish a statistically
valid result. Users remain responsible for checking observable selection,
parameter bounds, experimental correlations, covariance conditioning and the
applicability of asymptotic confidence-level approximations.
