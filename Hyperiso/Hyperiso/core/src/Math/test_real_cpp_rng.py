#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
from pathlib import Path
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def equicorr(n: int, rho: float) -> np.ndarray:
    R = np.full((n, n), rho, dtype=float)
    np.fill_diagonal(R, 1.0)
    return R


def load_matrix_txt(path: Path) -> np.ndarray:
    txt = Path(path).read_text().strip().split()
    it = iter(txt)
    n = int(next(it))
    vals = [float(next(it)) for _ in range(n * n)]
    R = np.array(vals, dtype=float).reshape(n, n)
    return R


def matrix_to_stdin_payload(R: np.ndarray) -> str:
    n = R.shape[0]
    lines = [str(n)]
    for i in range(n):
        lines.append(" ".join(f"{float(x)}" for x in R[i]))
    return "\n".join(lines) + "\n"



def compile_cpp_if_needed(bin_path: Path, cpp_path: Path) -> None:
    if bin_path.exists():
        return
    if not cpp_path or not cpp_path.exists():
        raise FileNotFoundError(
            f"Binaire introuvable ({bin_path}) et source C++ non trouvée ({cpp_path}). "
            f"Passe --cpp <chemin/vers/correlated_rng.cpp> ou --bin <chemin/vers/binaire>."
        )
    print(f"[INFO] Compilation de {cpp_path} → {bin_path} ...")
    cmd = ["g++", "-std=c++17", "-O2", "-Wall", str(cpp_path), "-o", str(bin_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError("Échec de la compilation. Voir erreurs ci-dessus.")
    print("[INFO] Compilation OK.")


def run_cpp_once(bin_path: Path, R: np.ndarray, dist: str = "gaussian", seed: int = None) -> np.ndarray:
    payload = matrix_to_stdin_payload(R)
    args = [str(bin_path), dist]
    if seed is not None:
        args.append(str(int(seed)))

    proc = subprocess.run(args, input=payload, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Exécution C++ a échoué (code {proc.returncode}): {proc.stderr.strip()}")

    try:
        vals = [float(tok) for tok in proc.stdout.strip().split()]
    except ValueError:
        raise RuntimeError(f"Sortie C++ non numérique:\n{proc.stdout}\nERR:\n{proc.stderr}")
    return np.array(vals, dtype=float)


def sample_many(bin_path: Path, R: np.ndarray, n_samples: int, seed_base: int, dist: str = "gaussian") -> np.ndarray:
    """Appelle le binaire n_samples fois (seed = seed_base + i) et empile les réalisations."""
    ys: List[np.ndarray] = []
    for i in range(n_samples):
        y = run_cpp_once(bin_path, R, dist=dist, seed=seed_base + i)
        ys.append(y)
    return np.vstack(ys)


def compute_stats(Y: np.ndarray, R: np.ndarray) -> Dict[str, float]:
    emp_corr = np.corrcoef(Y, rowvar=False)
    err = emp_corr - R
    means = Y.mean(axis=0)
    stds = Y.std(axis=0, ddof=1)
    return {
        "samples": Y.shape[0],
        "dim": Y.shape[1],
        "max|corr_error|": float(np.max(np.abs(err))),
        "mean|corr_error|": float(np.mean(np.abs(err))),
        "max|mean|": float(np.max(np.abs(means))),
        "max|std-1|": float(np.max(np.abs(stds - 1.0))),
    }


def plot_case(save_dir: Path, name: str, Y: np.ndarray, R: np.ndarray):
    save_dir.mkdir(parents=True, exist_ok=True)
    emp_corr = np.corrcoef(Y, rowvar=False)

    plt.figure()
    plt.imshow(emp_corr, vmin=-1, vmax=1)
    plt.title(f"Empirical correlation – {name}")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(save_dir / f"{name}_emp_corr.png", dpi=160)
    plt.close()

    plt.figure()
    plt.imshow(R, vmin=-1, vmax=1)
    plt.title(f"Target correlation – {name}")
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(save_dir / f"{name}_target_corr.png", dpi=160)
    plt.close()

    subset = Y[: min(4000, len(Y)), :]
    plt.figure()
    plt.scatter(subset[:, 0], subset[:, 1], s=5)
    plt.title(f"Scatter y1 vs y2 – {name}")
    plt.xlabel("y1")
    plt.ylabel("y2")
    plt.tight_layout()
    plt.savefig(save_dir / f"{name}_scatter_y1_y2.png", dpi=160)
    plt.close()

    xs = np.linspace(Y[:, 0].min(), Y[:, 0].max(), 400)
    pdf = 1.0 / np.sqrt(2 * np.pi) * np.exp(-0.5 * xs * xs)
    plt.figure()
    plt.hist(Y[:, 0], bins=60, density=True)
    plt.plot(xs, pdf)
    plt.title(f"Histogram(y1) + N(0,1) – {name}")
    plt.xlabel("y1")
    plt.ylabel("density")
    plt.tight_layout()
    plt.savefig(save_dir / f"{name}_hist_y1.png", dpi=160)
    plt.close()



def build_cases(args) -> List[Tuple[str, np.ndarray]]:
    cases: List[Tuple[str, np.ndarray]] = []
    requested = [c.strip().lower() for c in args.cases.split(",")]

    for c in requested:
        if c == "identity":
            cases.append(("identity", np.eye(args.n_dim, dtype=float)))
        elif c.startswith("equicorr"):
            rho = args.rho
            cases.append((f"equicorr_rho_{rho}".replace(".", "_"), equicorr(args.n_dim, rho)))
        elif c == "block":
            R = np.array([[1.0, args.block_rho, 0.0],
                          [args.block_rho, 1.0, 0.0],
                          [0.0, 0.0, 1.0]], dtype=float)
            cases.append((f"block_rho12_{args.block_rho}".replace(".", "_"), R))
        elif c == "near":
            rho = args.near_rho
            cases.append((f"near_singular_rho_{rho}".replace(".", "_"), equicorr(args.n_dim, rho)))
        elif c == "custom":
            if not args.custom_matrix:
                raise ValueError("--cases inclut 'custom' mais --custom-matrix est vide.")
            R = load_matrix_txt(Path(args.custom_matrix))
            cases.append(("custom", R))
        else:
            raise ValueError(f"Cas inconnu: {c}")
    return cases


def main():
    p = argparse.ArgumentParser(description="Testeur Python pour le générateur C++ 'correlated_rng'.")
    p.add_argument("--bin", type=str, default="./correlated_rng",
                   help="Chemin du binaire C++ (défaut: ./correlated_rng).")
    p.add_argument("--cpp", type=str, default="",
                   help="Chemin du fichier source C++ (optionnel). S'il est fourni et que --bin n'existe pas, on compile.")
    p.add_argument("--dist", type=str, default="gaussian", help="Nom de distribution à passer au binaire.")
    p.add_argument("--samples", type=int, default=15000, help="Nombre d'échantillons par cas.")
    p.add_argument("--seed-base", type=int, default=12345, help="Seed de base (on utilise seed_base+i).")
    p.add_argument("--cases", type=str, default="identity,equicorr,block,near",
                   help="Liste de cas séparés par des virgules parmi: identity,equicorr,block,near,custom")
    p.add_argument("--n-dim", type=int, default=3, help="Dimension pour identity/equicorr/near.")
    p.add_argument("--rho", type=float, default=0.5, help="Rho pour equicorr.")
    p.add_argument("--near-rho", type=float, default=0.999, help="Rho pour near (quasi-singulier).")
    p.add_argument("--block-rho", type=float, default=0.8, help="Corrélation entre y1 et y2 dans le cas 'block'.")
    p.add_argument("--custom-matrix", type=str, default="", help="Chemin vers une matrice custom (format C++).")
    p.add_argument("--out-dir", type=str, default="./cpp_test_out", help="Dossier de sortie (CSV + PNG).")
    p.add_argument("--no-plots", action="store_true", help="Ne pas générer de figures.")
    args = p.parse_args()

    bin_path = Path(args.bin)
    cpp_path = Path(args.cpp) if args.cpp else None
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not bin_path.exists() and cpp_path is not None:
        compile_cpp_if_needed(bin_path, cpp_path)

    if not bin_path.exists():
        raise FileNotFoundError(
            f"Binaire introuvable: {bin_path}\n"
            f"→ Compile-le (g++ -std=c++17 -O2 correlated_rng.cpp -o correlated_rng)\n"
            f"ou passe --cpp pour compiler automatiquement."
        )

    cases = build_cases(args)

    rows = []
    for idx, (name, R) in enumerate(cases, 1):
        print(f"[{idx}/{len(cases)}] Cas '{name}' : génération de {args.samples} échantillons …")
        Y = sample_many(bin_path, R, n_samples=args.samples, seed_base=args.seed_base, dist=args.dist)
        stats = compute_stats(Y, R)
        row = {"case": name, **stats}
        rows.append(row)

        np.save(out_dir / f"{name}_samples.npy", Y)
        np.savetxt(out_dir / f"{name}_target_R.txt", R, fmt="%.6f")
        np.savetxt(out_dir / f"{name}_empirical_corr.txt", np.corrcoef(Y, rowvar=False), fmt="%.6f")

        if not args.no_plots:
            plot_case(out_dir, name, Y, R)

    df = pd.DataFrame(rows).set_index("case")
    csv_path = out_dir / "summary.csv"
    df.to_csv(csv_path, float_format="%.6g")

    print("\n=== Résumé ===")
    print(df.to_string())
    print(f"\nRésultats enregistrés dans: {out_dir.resolve()}")
    print(f"Résumé CSV: {csv_path.resolve()}")


if __name__ == "__main__":
    main()
