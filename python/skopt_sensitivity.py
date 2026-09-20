"""Per-parameter sensitivity readout for the Phase 1 skopt campaign.

Reads output/skopt_trials.csv (written by optimize_phase1.py) and, for each of the search
dimensions, reports how much RMSE_all moves with it, how tightly the GP has converged on it,
and whether trials are hitting a range edge. Uses several independent readouts because n is
small (20 trials, 10 of them random) and any single one is easy to over-read:

  * Spearman / Pearson correlation and a quadratic-fit R^2 of RMSE_all against each parameter
  * GP with a per-parameter (ARD) Matern length scale, fitted to log(RMSE_all): a length scale
    at its upper bound (50) means the GP sees no dependence on that parameter
  * Random-forest permutation importance
  * Spread of each parameter in the random phase (trials 1-10) vs the GP-guided phase (11-20)

Search bounds come from optimize_phase1.SPACE, so they can't drift from the campaign.

Run from the repo root:
    .venv/bin/python python/skopt_sensitivity.py
"""
import csv
import warnings
from pathlib import Path

import numpy as np
import scipy.stats as st
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from optimize_phase1 import SPACE

REPO_ROOT = Path(__file__).resolve().parent.parent
TRIALS_CSV = REPO_ROOT / "output" / "skopt_trials.csv"
NAMES = [dim.name for dim in SPACE]
LO = np.array([dim.bounds[0] for dim in SPACE], dtype=float)
HI = np.array([dim.bounds[1] for dim in SPACE], dtype=float)
N_RANDOM = 10  # gp_minimize's default n_initial_points


def load_trials(path=TRIALS_CSV):
    """Return raw parameter values X, unit-cube positions U (0..1 within bounds), RMSE_all, RMSE_me3."""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    X = np.array([[float(r[n]) for n in NAMES] for r in rows])
    y = np.array([float(r["RMSE_all"]) for r in rows])
    y_me3 = np.array([float(r["RMSE_me3"]) for r in rows])
    return X, (X - LO) / (HI - LO), y, y_me3


def main():
    warnings.filterwarnings("ignore")
    X, U, y, y_me3 = load_trials()
    n = len(y)
    order = np.argsort(y)
    top = order[:6]
    print(f"n={n}  RMSE_all min/median/max = {y.min():.3f}/{np.median(y):.3f}/{y.max():.2f}")
    print("top6 trials:", [int(i) + 1 for i in top], y[top].round(3))

    print("\nTrial table (u = position within bounds, 0..1):")
    print("  #  RMSE_all me3   " + " ".join(f"{nm[:9]:>9s}" for nm in NAMES))
    for i in range(n):
        print(f"{i + 1:3d}  {y[i]:7.2f} {y_me3[i]:5.2f}  " + " ".join(f"{v:9.2f}" for v in U[i]))

    print(f"\n{'param':19s} {'spear':>6s} {'p':>5s} {'quadR2':>6s} {'u_min':>5s} {'u_max':>5s} "
          f"{'sd rand':>7s} {'sd GP':>6s} {'top6mean':>8s} {'top6sd':>6s}")
    for j, nm in enumerate(NAMES):
        u = U[:, j]
        sp, p = st.spearmanr(u, y)
        A = np.c_[np.ones_like(u), u, u ** 2]
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        r2 = 1 - ((y - A @ coef) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        print(f"{nm:19s} {sp:6.2f} {p:5.2f} {r2:6.2f} {u.min():5.2f} {u.max():5.2f} "
              f"{u[:N_RANDOM].std():7.3f} {u[N_RANDOM:].std():6.3f} {u[top].mean():8.2f} {u[top].std():6.3f}")

    print("\nRaw values, top-6:")
    for j, nm in enumerate(NAMES):
        digits = 3 if X[:, j].max() < 2 else 0
        print(f"{nm:19s}", X[top, j].round(digits).tolist(), "bounds", (LO[j], HI[j]))

    print("\nEdge hits (u<=.03 / u>=.97) over all trials ; over top6:")
    for j, nm in enumerate(NAMES):
        print(f"{nm:19s} all: {(U[:, j] <= .03).sum()}/{(U[:, j] >= .97).sum()}  "
              f"top6: {(U[top, j] <= .03).sum()}/{(U[top, j] >= .97).sum()}  best u={U[order[0], j]:.3f}")

    # GP with a per-parameter (ARD) length scale on log RMSE_all
    ly = np.log(y)
    ly_std = (ly - ly.mean()) / ly.std()
    best = None
    for seed in range(8):
        kernel = (ConstantKernel(1.0, (1e-2, 1e2))
                  * Matern(length_scale=np.ones(len(NAMES)), length_scale_bounds=(0.05, 50), nu=2.5)
                  + WhiteKernel(0.05, (1e-4, 1)))
        gp = GaussianProcessRegressor(kernel, n_restarts_optimizer=3, random_state=seed).fit(U, ly_std)
        if best is None or gp.log_marginal_likelihood_value_ > best.log_marginal_likelihood_value_:
            best = gp
    print(f"\nGP (ARD Matern-2.5, on log RMSE_all); fitted noise = {best.kernel_.k2.noise_level:.3f}")
    print("length scales (unit cube; small = matters, at the 50 cap = GP sees no dependence):")
    for nm, ls in sorted(zip(NAMES, best.kernel_.k1.k2.length_scale), key=lambda t: t[1]):
        print(f"  {nm:19s} {ls:7.2f}")

    # random-forest permutation importance (mean over 30 forests)
    imp = np.zeros(len(NAMES))
    for s in range(30):
        rf = RandomForestRegressor(200, min_samples_leaf=2, random_state=s, oob_score=True).fit(U, ly)
        rng = np.random.RandomState(s)
        for j in range(len(NAMES)):
            Up = U.copy()
            Up[:, j] = rng.permutation(Up[:, j])
            imp[j] += ((rf.predict(Up) - rf.predict(U)) ** 2).mean()
    imp /= 30
    print("\nRF permutation importance (mean sq change in predicted log-RMSE):")
    for nm, v in sorted(zip(NAMES, imp), key=lambda t: -t[1]):
        print(f"  {nm:19s} {v:7.4f}")
    print("RF OOB R2 (log RMSE_all):", round(rf.oob_score_, 2))

    print(f"\nRange of u over the GP-guided trials ({N_RANDOM + 1}-{n}), min..max:")
    for j, nm in enumerate(NAMES):
        print(f"  {nm:19s} {U[N_RANDOM:, j].min():.2f}..{U[N_RANDOM:, j].max():.2f}")


if __name__ == "__main__":
    main()
