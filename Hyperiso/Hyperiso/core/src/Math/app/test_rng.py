import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt, pi, exp


rng = np.random.default_rng(42)

def generate_correlated(R, n_samples, rng):
    R = np.asarray(R, dtype=float)
    L = np.linalg.cholesky(R)
    z = rng.standard_normal(size=(n_samples, R.shape[0]))
    y = z @ L.T
    return y, L

def summarize_case(name, R, n_samples=50000):
    y, L = generate_correlated(R, n_samples, rng)
    sample_corr = np.corrcoef(y, rowvar=False)
    err = sample_corr - R
    means = y.mean(axis=0)
    stds = y.std(axis=0, ddof=1)
    summary = {
        "case": name,
        "n": R.shape[0],
        "samples": n_samples,
        "max|corr_error|": float(np.max(np.abs(err))),
        "mean|corr_error|": float(np.mean(np.abs(err))),
        "max|mean|": float(np.max(np.abs(means))),
        "max|std-1|": float(np.max(np.abs(stds - 1.0))),
        "min_cholesky_diag": float(np.min(np.diag(np.linalg.cholesky(R)))),
    }
    return y, L, sample_corr, err, means, stds, summary

def equicorrelation(n, rho):
    R = np.full((n, n), rho, dtype=float)
    np.fill_diagonal(R, 1.0)
    return R

cases = [
    ("Identity (no correlation)", np.eye(3)),
    ("Equicorrelation rho=0.5", equicorrelation(3, 0.5)),
    ("Block: rho12=0.8, var3 indep", np.array([[1.0, 0.8, 0.0],
                                                [0.8, 1.0, 0.0],
                                                [0.0, 0.0, 1.0]], dtype=float)),
    ("Nearly singular equicorr rho=0.999", equicorrelation(3, 0.999)),
]

all_summaries = []
results = {}

for name, R in cases:
    y, L, sample_corr, err, means, stds, summary = summarize_case(name, R, n_samples=60000)
    all_summaries.append(summary)
    results[name] = {
        "R": R,
        "L": L,
        "y": y,
        "sample_corr": sample_corr,
        "err": err,
        "means": means,
        "stds": stds,
    }

summary_df = pd.DataFrame(all_summaries).set_index("case")

for name in results:
    df_corr = pd.DataFrame(results[name]["sample_corr"])
    df_err = pd.DataFrame(results[name]["err"])

for name, data in results.items():
    sample_corr = data["sample_corr"]
    y = data["y"]
    R = data["R"]
    
    plt.figure()
    plt.imshow(sample_corr, vmin=-1, vmax=1)
    plt.title(f"Empirical correlation matrix – {name}")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.imshow(R, vmin=-1, vmax=1)
    plt.title(f"Target correlation matrix – {name}")
    plt.colorbar()
    plt.tight_layout()
    plt.show()
    
    subset = y[:5000, :]
    plt.figure()
    plt.scatter(subset[:, 0], subset[:, 1], s=5)
    plt.title(f"Scatter: y1 vs y2 – {name}")
    plt.xlabel("y1")
    plt.ylabel("y2")
    plt.tight_layout()
    plt.show()
    
    vals = y[:, 0]
    plt.figure()
    plt.hist(vals, bins=60, density=True)
    xs = np.linspace(vals.min(), vals.max(), 400)
    pdf = (1.0 / sqrt(2 * pi)) * np.exp(-0.5 * xs * xs)
    plt.plot(xs, pdf)
    plt.title(f"Histogram of y1 (with N(0,1) pdf) – {name}")
    plt.xlabel("y1")
    plt.ylabel("density")
    plt.tight_layout()
    plt.show()
