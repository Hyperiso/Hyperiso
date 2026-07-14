from __future__ import annotations

import os
import math
from typing import Any, Tuple

import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist

from pyhyperiso.core.Math.RealMatrix import Matrix


from pyhyperiso.core.Statistic.JointDistribution import JointDistributionFactory as JF

from pyhyperiso.core.Statistic.MarginalDistribution import DistributionFactoryWrapper as MF

from pyhyperiso.core.Statistic.MarginalDistribution import MarginalKind

from pyhyperiso.core.Statistic.MarginalConfig import (
    GaussianMarginalConfig as GaussianCfg,
    SplitGaussianMarginalConfig as SplitGaussianCfg,
    LikelihoodMarginalConfig as LikelihoodCfg,
)

from pyhyperiso.core.Statistic.Copula import CopulaKind

from pyhyperiso.core.Statistic.CopulaConfig import (
    GaussianCopulaConfigPy as GaussianCopulaCfg,
    StudentTCopulaConfigPy as StudentTCopulaCfg,
)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def ar1_corr(d: int, rho: float) -> Matrix:
    if not (-0.999 < rho < 0.999):
        raise ValueError("rho doit être dans (-1,1).")
    R = [[(rho ** abs(i - j)) for j in range(d)] for i in range(d)]
    return Matrix(R)


def make_likelihood_cfg() -> Any:
    xs = np.linspace(-5.0, 5.0, 401)
    w = np.exp(-0.5 * ((xs - (-1.2)) / 0.7) ** 2) + 0.65 * np.exp(-0.5 * ((xs - 1.6) / 1.1) ** 2)
    w = w / np.sum(w)
    return LikelihoodCfg(values=xs.tolist(), weights=w.tolist())


def _safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def ecdf(samples: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    s = np.sort(samples)
    y = np.arange(1, len(s) + 1) / len(s)
    return s, y


def gaussianize(U: np.ndarray) -> np.ndarray:
    nd = NormalDist()
    eps = 1e-12
    Uc = np.clip(U, eps, 1 - eps)
    inv = np.vectorize(nd.inv_cdf, otypes=[float])
    return inv(Uc)


def spearman_corr(U: np.ndarray) -> np.ndarray:
    order = np.argsort(U, axis=0)
    ranks = np.empty_like(order, dtype=float)
    for j in range(U.shape[1]):
        ranks[order[:, j], j] = np.arange(1, len(U) + 1, dtype=float)
    return np.corrcoef(ranks, rowvar=False)


def plot_hist_pdf(samples: np.ndarray, dist, name: str, outpath: str) -> None:
    # grille basée sur quantiles pour être “smart”
    try:
        xlo = _safe_float(dist.ppf(0.001))
        xhi = _safe_float(dist.ppf(0.999))
        if not (np.isfinite(xlo) and np.isfinite(xhi) and xhi > xlo):
            raise ValueError
        pad = 0.08 * (xhi - xlo)
        xs = np.linspace(xlo - pad, xhi + pad, 700)
    except Exception:
        xs = np.linspace(np.min(samples), np.max(samples), 700)

    logp = np.array([_safe_float(dist.logpdf(float(x))) for x in xs], dtype=float)
    pdf = np.exp(logp)

    plt.figure(figsize=(10.2, 5.6))
    plt.hist(samples, bins=60, density=True, alpha=0.35, label="samples (density)")
    plt.plot(xs, pdf, linewidth=2.3, label="marginal pdf")
    plt.title(f"{name} — marginal PDF vs samples")
    plt.xlabel("x")
    plt.ylabel("density")
    plt.grid(True, alpha=0.28)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_cdf_ecdf(samples: np.ndarray, dist, name: str, outpath: str) -> None:
    try:
        xlo = _safe_float(dist.ppf(0.001))
        xhi = _safe_float(dist.ppf(0.999))
        if not (np.isfinite(xlo) and np.isfinite(xhi) and xhi > xlo):
            raise ValueError
        pad = 0.08 * (xhi - xlo)
        xs = np.linspace(xlo - pad, xhi + pad, 700)
    except Exception:
        xs = np.linspace(np.min(samples), np.max(samples), 700)

    cdf_th = np.array([_safe_float(dist.cdf(float(x))) for x in xs], dtype=float)
    s, y = ecdf(samples)

    plt.figure(figsize=(10.2, 5.6))
    plt.plot(xs, cdf_th, linewidth=2.3, label="theoretical CDF")
    plt.step(s, y, where="post", linewidth=1.8, alpha=0.9, label="ECDF")
    plt.title(f"{name} — marginal CDF vs ECDF")
    plt.xlabel("x")
    plt.ylabel("CDF(x)")
    plt.ylim(-0.02, 1.02)
    plt.grid(True, alpha=0.28)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_uniformity(u: np.ndarray, name: str, outpath: str) -> None:
    plt.figure(figsize=(7.6, 4.6))
    plt.hist(u, bins=50, density=True, alpha=0.45, label="u = F(x)")
    plt.plot([0, 1], [1, 1], linestyle="--", linewidth=1.6, label="Uniform(0,1) density")
    plt.title(f"{name} — PIT check (u should be uniform)")
    plt.xlabel("u")
    plt.ylabel("density")
    plt.xlim(0, 1)
    plt.ylim(bottom=0)
    plt.grid(True, alpha=0.25)
    plt.minorticks_on()
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_hist2d(u1: np.ndarray, u2: np.ndarray, title: str, outpath: str) -> None:
    plt.figure(figsize=(6.9, 6.5))
    plt.hist2d(u1, u2, bins=80, density=True)
    plt.colorbar(label="density")
    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1.2)
    plt.title(title)
    plt.xlabel("dim 1")
    plt.ylabel("dim 2")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.22)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_heatmap(M: np.ndarray, title: str, outpath: str) -> None:
    plt.figure(figsize=(6.9, 6.0))
    im = plt.imshow(M, vmin=-1, vmax=1, interpolation="nearest")
    plt.colorbar(im, label="correlation")
    plt.title(title)
    plt.xlabel("dim")
    plt.ylabel("dim")

    n = M.shape[0]
    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_joint_density_2d(jd, X: np.ndarray, title: str, outpath: str, grid: int = 120) -> None:
    # heatmap exp(logpdf) sur une grille autour des quantiles empiriques
    x1 = X[:, 0]
    x2 = X[:, 1]
    q1 = np.quantile(x1, [0.002, 0.998])
    q2 = np.quantile(x2, [0.002, 0.998])
    pad1 = 0.08 * (q1[1] - q1[0])
    pad2 = 0.08 * (q2[1] - q2[0])

    xs = np.linspace(q1[0] - pad1, q1[1] + pad1, grid)
    ys = np.linspace(q2[0] - pad2, q2[1] + pad2, grid)

    D = np.empty((grid, grid), dtype=float)
    for i, y in enumerate(ys):
        row = []
        for x in xs:
            lp = float(jd.logpdf([float(x), float(y)]))
            row.append(math.exp(lp))
        D[i, :] = row

    plt.figure(figsize=(7.6, 6.2))
    im = plt.imshow(
        D, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], interpolation="nearest"
    )
    plt.colorbar(im, label="joint density")
    plt.title(title)
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()



def main(
    outdir: str = "joint_plots",
    seed: int = 123,
    n: int = 15000,
    copula_kind: str = "GAUSSIAN",  # "GAUSSIAN" or "STUDENT_T"
    nu: int = 6,
    rho: float = 0.65,
) -> None:
    ensure_dir(outdir)

    marginal_types = [
        MarginalKind.GAUSSIAN,
        MarginalKind.SPLIT_GAUSSIAN
        if hasattr(MarginalKind, "SPLIT_GAUSSIAN")
        else MarginalKind.HALF_GAUSSIAN,
    ]

    marginal_cfgs = [
        GaussianCfg(mu=0.6, sigma=1.15),
        SplitGaussianCfg(mu=0.0, sigma_p=0.75, sigma_m=1.55),
    ]

    d = len(marginal_types)

    marginals = [MF.create(marginal_types[i], marginal_cfgs[i], seed=seed + i) for i in range(d)]

    R = ar1_corr(d, rho)
    if copula_kind.upper() == "GAUSSIAN":
        c_kind = CopulaKind.GAUSSIAN
        c_cfg = GaussianCopulaCfg(R=R)
    else:
        c_kind = CopulaKind.STUDENT_T
        c_cfg = StudentTCopulaCfg(R=R, nu=int(nu))

    jd = JF.create(marginal_types, marginal_cfgs, c_kind, c_cfg, seed=seed)

    X = np.array(jd.sample(int(n)), dtype=float)  # shape (n, d)

    for i in range(d):
        xi = X[:, i]
        dist = marginals[i]

        plot_hist_pdf(
            xi,
            dist,
            name=f"Marginal {i + 1}",
            outpath=os.path.join(outdir, f"marg_{i + 1}_pdf.png"),
        )
        plot_cdf_ecdf(
            xi,
            dist,
            name=f"Marginal {i + 1}",
            outpath=os.path.join(outdir, f"marg_{i + 1}_cdf.png"),
        )

        ui = np.array([_safe_float(dist.cdf(float(x))) for x in xi], dtype=float)
        plot_uniformity(
            ui, name=f"Marginal {i + 1}", outpath=os.path.join(outdir, f"marg_{i + 1}_pit.png")
        )

    U = np.empty_like(X, dtype=float)
    for i in range(d):
        dist = marginals[i]
        U[:, i] = np.array([_safe_float(dist.cdf(float(x))) for x in X[:, i]], dtype=float)
    U = np.clip(U, 1e-12, 1 - 1e-12)

    if d >= 2:
        plot_hist2d(
            U[:, 0],
            U[:, 1],
            title="Dependence check — copula space (u1,u2)",
            outpath=os.path.join(outdir, "copula_u1_u2_hist2d.png"),
        )

    S = spearman_corr(U)
    Z = gaussianize(U)
    C = np.corrcoef(Z, rowvar=False)

    plot_heatmap(S, "Spearman corr of U (ranks)", os.path.join(outdir, "copula_spearman.png"))
    plot_heatmap(
        C,
        "Corr of Z = Phi^-1(U) (compare to target R)",
        os.path.join(outdir, "copula_gaussianized_corr.png"),
    )

    if d == 2:
        plot_joint_density_2d(
            jd,
            X,
            title="Joint density heatmap: exp(logpdf(x))",
            outpath=os.path.join(outdir, "joint_density_heatmap.png"),
            grid=120,
        )

    print(f"[OK] Saved plots in: {outdir}/")
    print(
        f"Joint dim = {jd.dim()} ; mean logpdf(sample) ≈ {float(np.mean([jd.logpdf(row.tolist()) for row in X[:2000]])):.3f}"
    )
    if copula_kind.upper() == "GAUSSIAN":
        print("Tip: for Gaussian copula, Corr(Phi^-1(U)) should be close to your R.")


if __name__ == "__main__":
    main(
        outdir="joint_plots",
        seed=123,
        n=15000,
        copula_kind="GAUSSIAN",  # or "STUDENT_T"
        nu=6,
        rho=0.65,
    )
