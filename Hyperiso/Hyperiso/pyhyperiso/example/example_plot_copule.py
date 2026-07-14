from __future__ import annotations

import os
import math
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist

from pyhyperiso.core.Math.RealMatrix import Matrix


from pyhyperiso.core.Statistic.Copula import CopulaKind 

from pyhyperiso.core.Statistic.Copula import CopulaFactoryWrapper as CF

from pyhyperiso.core.Statistic.CopulaConfig import GaussianCopulaConfigPy as GaussianCfg

from pyhyperiso.core.Statistic.CopulaConfig import StudentTCopulaConfigPy as StudentTCfg


def ar1_corr(d: int, rho: float) -> Matrix:
    """R_ij = rho^{|i-j|}. SPD si |rho|<1."""
    if not (-0.999 < rho < 0.999):
        raise ValueError("rho doit être dans (-1,1) pour une corrélation AR(1) SPD.")
    R = [[(rho ** abs(i - j)) for j in range(d)] for i in range(d)]
    return Matrix(R)


def _rankdata(x: np.ndarray) -> np.ndarray:
    """Ranks 1..n (mid-rank approx via argsort; suffisant pour tests visuels)."""
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    return ranks


def spearman_corr(U: np.ndarray) -> np.ndarray:
    """Corrélation de Spearman via corrcoef des rangs."""
    R = np.vstack([_rankdata(U[:, i]) for i in range(U.shape[1])]).T
    return np.corrcoef(R, rowvar=False)


def gaussianize(U: np.ndarray) -> np.ndarray:
    """Transforme U~(0,1) vers Z~N(0,1) via invCDF (sans scipy)."""
    nd = NormalDist()
    # clamp pour éviter inf
    eps = 1e-12
    Uc = np.clip(U, eps, 1 - eps)
    inv = np.vectorize(nd.inv_cdf, otypes=[float])
    Z = inv(Uc)
    return Z


def corrcoef_Z(U: np.ndarray) -> np.ndarray:
    Z = gaussianize(U)
    return np.corrcoef(Z, rowvar=False)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)




def plot_scatter_or_hist2d(
    U: np.ndarray, title: str, outpath: str, use_hist2d: bool = False
) -> None:
    u1, u2 = U[:, 0], U[:, 1]
    plt.figure(figsize=(6.8, 6.4))
    if use_hist2d:
        plt.hist2d(u1, u2, bins=70, density=True)
        plt.colorbar(label="density")
    else:
        plt.plot(u1, u2, linestyle="none", marker=".", markersize=2, alpha=0.25)

    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1.2)  # diagonal
    plt.title(title)
    plt.xlabel("u1")
    plt.ylabel("u2")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_marginal_uniformity(U: np.ndarray, title_prefix: str, outdir: str) -> None:
    d = U.shape[1]
    for j in range(d):
        plt.figure(figsize=(7.6, 4.6))
        plt.hist(U[:, j], bins=50, density=True, alpha=0.45, label=f"u{j + 1}")

        plt.plot([0, 1], [1, 1], linestyle="--", linewidth=1.6, label="Uniform(0,1) density")
        plt.title(f"{title_prefix} — Marginal uniformity (dim {j + 1})")
        plt.xlabel(f"u{j + 1}")
        plt.ylabel("density")
        plt.xlim(0, 1)
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.25)
        plt.legend(frameon=True)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{title_prefix}_marginal_u{j + 1}.png"), dpi=200)
        plt.close()


def plot_corr_heatmap(M: np.ndarray, title: str, outpath: str) -> None:
    plt.figure(figsize=(6.8, 6.0))
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


def plot_density_surface_2d(copula: Any, title: str, outpath: str, grid: int = 90) -> None:
    """Visualise exp(log_density(u)) sur [0,1]^2. (Seulement 2D)."""
    eps = 1e-4
    xs = np.linspace(eps, 1 - eps, grid)
    ys = np.linspace(eps, 1 - eps, grid)
    D = np.empty((grid, grid), dtype=float)

    for i, y in enumerate(ys):
        row = []
        for x in xs:
            ld = float(copula.log_density([float(x), float(y)]))
            row.append(math.exp(ld))
        D[i, :] = row

    approx_integral = float(D.mean())  # area=1, grid uniforme
    plt.figure(figsize=(7.0, 6.0))
    im = plt.imshow(D, origin="lower", extent=[0, 1, 0, 1], interpolation="nearest")
    plt.colorbar(im, label="density c(u)")
    plt.title(f"{title}\n(mean density ≈ {approx_integral:.3f} ; should be ~1 in 2D)")
    plt.xlabel("u1")
    plt.ylabel("u2")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()




def make_gaussian_copula(R: Matrix, seed: int):
    cfg = GaussianCfg(R=R)
    if hasattr(CF, "gaussian"):
        return CF.gaussian(R, seed=seed)
    if hasattr(CF, "create"):
        return CF.create(CopulaKind.GAUSSIAN, cfg, seed=seed)
    return CF.create_gaussian(cfg, seed)  # type: ignore


def make_student_t_copula(R: Matrix, nu: int, seed: int):
    cfg = StudentTCfg(R=R, nu=int(nu))
    if hasattr(CF, "student_t"):
        return CF.student_t(R, nu=nu, seed=seed)
    if hasattr(CF, "create"):
        return CF.create(CopulaKind.STUDENT_T, cfg, seed=seed)
    return CF.create_student_t(cfg, seed)  # type: ignore


def main(
    outdir: str = "copula_plots",
    seed: int = 123,
    n: int = 12000,
    d: int = 3,
    rho: float = 0.65,
    nu: int = 6,
) -> None:
    ensure_dir(outdir)

    R = ar1_corr(d, rho)

    cop_gauss = make_gaussian_copula(R, seed=seed)
    cop_t = make_student_t_copula(R, nu=nu, seed=seed)

    specs = [
        ("GaussianCopula", cop_gauss),
        (f"StudentTCopula(nu={nu})", cop_t),
    ]

    for name, cop in specs:
        U = np.array(cop.sample_u(int(n)), dtype=float)  # shape (n, d)

        if d >= 2:
            plot_scatter_or_hist2d(
                U[:, :2],
                title=f"{name} — samples in (0,1)^2 (projection)",
                outpath=os.path.join(outdir, f"{name}_hist2d.png"),
                use_hist2d=False,
            )

            if d == 2:
                plot_density_surface_2d(
                    cop,
                    title=f"{name} — density surface in (0,1)^2",
                    outpath=os.path.join(outdir, f"{name}_density_surface.png"),
                    grid=90,
                )

        plot_marginal_uniformity(U, title_prefix=name, outdir=outdir)

        S = spearman_corr(U)
        C = corrcoef_Z(U)

        plot_corr_heatmap(
            S,
            title=f"{name} — Spearman corr (ranks of U)",
            outpath=os.path.join(outdir, f"{name}_spearman.png"),
        )
        plot_corr_heatmap(
            C,
            title=f"{name} — Corr of Z = Phi^-1(U) (compare to target R)",
            outpath=os.path.join(outdir, f"{name}_corr_gaussianized.png"),
        )

    print(f"[OK] Plots saved in: {outdir}/")


if __name__ == "__main__":
    main(
        outdir="copula_plots",
        seed=123,
        n=12000,
        d=3,
        rho=0.65,
        nu=6,
    )
