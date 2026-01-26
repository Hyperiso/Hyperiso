from __future__ import annotations

import os
import math
from dataclasses import dataclass
from typing import Sequence, Optional, List, Tuple

import numpy as np
import matplotlib.pyplot as plt


from pyhyperiso.core.Statistic.MarginalDistribution import (
    DistributionFactoryWrapper as DF,
    MarginalKind,
    GaussianMarginalConfig,
    SplitGaussianMarginalConfig,
    FlatMarginalConfig,
    LikelihoodMarginalConfig,
)


def make_likelihood_cfg() -> LikelihoodMarginalConfig:
    """Likelihood synthétique (mélange de 2 gaussiennes) sur une grille."""
    xs = np.linspace(-5.0, 5.0, 401)
    w = np.exp(-0.5 * ((xs - (-1.2)) / 0.7) ** 2) + 0.65 * np.exp(-0.5 * ((xs - 1.6)) / 1.1) ** 2
    w = w / np.sum(w)
    return LikelihoodMarginalConfig(values=xs.tolist(), weights=w.tolist())


@dataclass(frozen=True)
class DistSpec:
    name: str
    kind: MarginalKind
    cfg: object


def _safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def _grid_from_ppf(dist, lo=0.001, hi=0.999, n=600) -> np.ndarray:
    """Grille x basée sur des quantiles (robuste et adapté à la distrib)."""
    try:
        x_lo = _safe_float(dist.ppf(lo))
        x_hi = _safe_float(dist.ppf(hi))
        if not (np.isfinite(x_lo) and np.isfinite(x_hi) and x_hi > x_lo):
            raise ValueError("bad ppf range")
        pad = 0.08 * (x_hi - x_lo)
        return np.linspace(x_lo - pad, x_hi + pad, n)
    except Exception:
        return np.linspace(-6.0, 6.0, n)


def _pdf(dist, xs: np.ndarray) -> np.ndarray:
    logp = np.array([_safe_float(dist.logpdf(float(x))) for x in xs], dtype=float)
    return np.exp(logp)


def _cdf(dist, xs: np.ndarray) -> np.ndarray:
    return np.array([_safe_float(dist.cdf(float(x))) for x in xs], dtype=float)


def _ecdf(samples: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    s = np.sort(samples)
    y = np.arange(1, len(s) + 1) / len(s)
    return s, y


def plot_pdf_vs_samples(name: str, dist, xs: np.ndarray, samples: np.ndarray, outpath: Optional[str] = None) -> None:
    pdf = _pdf(dist, xs)

    plt.figure(figsize=(10.2, 5.7))
    plt.hist(samples, bins=60, density=True, alpha=0.35, label="samples (density)")
    plt.plot(xs, pdf, linewidth=2.3, label="pdf = exp(logpdf)")

    # Bande 68% + médiane
    try:
        q16, q50, q84 = dist.ppf(0.16), dist.ppf(0.50), dist.ppf(0.84)
        q16, q50, q84 = float(q16), float(q50), float(q84)
        if all(map(np.isfinite, [q16, q50, q84])) and (q84 > q16):
            plt.axvspan(q16, q84, alpha=0.12, label="68% interval [16%,84%]")
            plt.axvline(q50, linestyle="--", linewidth=1.8, label="median (50%)")
    except Exception:
        pass

    # Moyenne
    try:
        mu = float(dist.mean())
        if np.isfinite(mu):
            plt.axvline(mu, linestyle=":", linewidth=1.8, label="mean")
    except Exception:
        pass

    plt.title(f"{name} — PDF vs samples")
    plt.xlabel("x")
    plt.ylabel("density")
    plt.grid(True, alpha=0.28)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=200)
        plt.close()


def plot_cdf_theory_vs_empirical(name: str, dist, xs: np.ndarray, samples: np.ndarray, outpath: Optional[str] = None) -> None:
    cdf_th = _cdf(dist, xs)
    s, y = _ecdf(samples)

    plt.figure(figsize=(10.2, 5.7))
    plt.plot(xs, cdf_th, linewidth=2.4, label="theoretical CDF")
    plt.step(s, y, where="post", linewidth=1.8, alpha=0.9, label="empirical CDF (ECDF)")

    # Repères (quartiles)
    for p, ls in [(0.25, "--"), (0.50, "--"), (0.75, "--")]:
        try:
            q = float(dist.ppf(p))
            if np.isfinite(q):
                plt.axvline(q, linestyle=ls, linewidth=1.2)
        except Exception:
            pass

    plt.title(f"{name} — CDF: theoretical vs empirical")
    plt.xlabel("x")
    plt.ylabel("CDF(x)")
    plt.ylim(-0.02, 1.02)
    plt.grid(True, alpha=0.28)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=200)
        plt.close()


def plot_global_cdf_ppf_consistency(specs: List[DistSpec], seed: int, outpath: Optional[str] = None) -> None:
    p = np.linspace(0.001, 0.999, 500)

    plt.figure(figsize=(10.4, 5.8))
    for sp in specs:
        dist = DF.create(sp.kind, sp.cfg, seed=seed)
        err = []
        for pi in p:
            try:
                x = float(dist.ppf(float(pi)))
                err.append(_safe_float(dist.cdf(x)) - float(pi))
            except Exception:
                err.append(float("nan"))
        err = np.array(err, dtype=float)
        plt.plot(p, err, linewidth=2.1, label=sp.name)

    plt.plot([0, 1], [0, 0], linestyle="--", linewidth=1.8, label="0 (ideal)")
    plt.title("Sanity check — CDF(PPF(p)) − p (should be near 0)")
    plt.xlabel("p")
    plt.ylabel("error")
    plt.grid(True, alpha=0.28)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=200)
        plt.close()


def main(show: bool = True, outdir: Optional[str] = "marginal_plots", seed: int = 123, n_samples: int = 3500) -> None:
    specs: List[DistSpec] = [
        DistSpec("Gaussian", MarginalKind.GAUSSIAN, GaussianMarginalConfig(mu=0.6, sigma=1.15)),
        DistSpec("SplitGaussian", MarginalKind.HALF_GAUSSIAN, SplitGaussianMarginalConfig(mu=0.0, sigma_p=0.75, sigma_m=1.55)),
        DistSpec("Flat", MarginalKind.FLAT, FlatMarginalConfig(a=-2.5, b=3.2)),
        # DistSpec("Likelihood", MarginalKind.LIKELIHOOD, make_likelihood_cfg()),
    ]

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)

    for sp in specs:
        dist = DF.create(sp.kind, sp.cfg, seed=seed)

        xs = _grid_from_ppf(dist)
        samples = np.array(dist.rvs(n_samples), dtype=float)

        pdf_path = None if outdir is None else os.path.join(outdir, f"{sp.name}_pdf_samples.png")
        cdf_path = None if outdir is None else os.path.join(outdir, f"{sp.name}_cdf_ecdf.png")

        plot_pdf_vs_samples(sp.name, dist, xs, samples, outpath=pdf_path)
        plot_cdf_theory_vs_empirical(sp.name, dist, xs, samples, outpath=cdf_path)

    global_path = None if outdir is None else os.path.join(outdir, "global_cdf_ppf_consistency.png")
    plot_global_cdf_ppf_consistency(specs, seed=seed, outpath=global_path)

    if show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main(show=True, outdir="marginal_plots", seed=123, n_samples=3500)
