"""Can a subset of the Phase 1 search parameters be fixed without losing predictive power?

Companion to skopt_sensitivity.py. Fits a GP on log(RMSE_all) using only a chosen subset of the
parameters and scores it by leave-one-out (LOO) prediction error over the trials in
output/skopt_trials.csv. If the reduced GP predicts held-out trials as well as the full one, the
dropped parameters carry no detectable signal at this sample size and are candidates to fix.

Reports: input correlations in the random phase (confounding check), LOO error for a battery of
subsets, drop-one-from-full deltas (positive = that parameter carries signal), and rank
correlations within the converged basin. Pass --free to score one specific subset.

This tests predictive value under the *placeholder* targets used in Phase 1; it says nothing about
whether a parameter will matter once real targets are used.

Run from the repo root:
    .venv/bin/python python/skopt_subset_loo.py
    .venv/bin/python python/skopt_subset_loo.py --free tfK27,tfK9,K9dn,K4dn,spread
"""
import argparse
import itertools
import warnings

import numpy as np
import scipy.stats as st
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from skopt_sensitivity import N_RANDOM, load_trials

# Short labels, in optimize_phase1.PHASE1_SPACE order.
SHORT = ["K4up", "K4dn", "K9up", "K9dn", "K27up", "K27dn", "spread", "tfK9", "tfK27"]
IDX = {s: i for i, s in enumerate(SHORT)}
BASIN_RMSE = 4.2  # trials below this count as "in the converged basin"

BATTERY = {
    "all 9": SHORT,
    "top5: tfK27,spread,K27up,tfK9,K9dn": ["tfK27", "spread", "K27up", "tfK9", "K9dn"],
    "top4: tfK27,spread,K27up,tfK9": ["tfK27", "spread", "K27up", "tfK9"],
    "6: top5 + K4dn": ["tfK27", "spread", "K27up", "tfK9", "K9dn", "K4dn"],
    "3: tfK27,spread,K27up": ["tfK27", "spread", "K27up"],
    "tfK27 only": ["tfK27"],
    "5 free: tfK27,tfK9,K9dn,K4dn,spread": ["tfK27", "tfK9", "K9dn", "K4dn", "spread"],
    "4 free: tfK27,tfK9,K9dn,K4dn": ["tfK27", "tfK9", "K9dn", "K4dn"],
    "6 free: above 5 + K27up": ["tfK27", "tfK9", "K9dn", "K4dn", "spread", "K27up"],
}


def make_gp(n_dims):
    kernel = (ConstantKernel(1.0, (1e-2, 1e2))
              * Matern(length_scale=np.ones(n_dims), length_scale_bounds=(0.05, 50), nu=2.5)
              + WhiteKernel(0.05, (1e-4, 1)))
    return GaussianProcessRegressor(kernel, n_restarts_optimizer=2, random_state=0)


def loo_rmse(U, ly, cols):
    """LOO prediction RMSE (in log-RMSE_all units) of a GP restricted to the given columns."""
    errs = []
    for i in range(len(ly)):
        keep = np.arange(len(ly)) != i
        mu, sd = ly[keep].mean(), ly[keep].std()
        gp = make_gp(len(cols)).fit(U[keep][:, cols], (ly[keep] - mu) / sd)
        errs.append(gp.predict(U[i:i + 1, cols])[0] * sd + mu - ly[i])
    return float(np.sqrt((np.array(errs) ** 2).mean()))


def cols_of(labels):
    return [IDX[s] for s in labels]


def main():
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--free", help="comma-separated short labels to score as a single subset; labels: "
                    + ",".join(SHORT))
    args = ap.parse_args()

    _, U, y, _ = load_trials()
    ly = np.log(y)

    if args.free:
        labels = args.free.split(",")
        print(f"LOO RMSE for [{args.free}]: {loo_rmse(U, ly, cols_of(labels)):.3f}  "
              f"(baseline sd of log y = {ly.std():.2f})")
        return

    print(f"Input correlations, random-phase trials 1-{N_RANDOM} only, |rho|>=0.5:")
    corr = np.corrcoef(U[:N_RANDOM].T)
    for i, j in itertools.combinations(range(len(SHORT)), 2):
        if abs(corr[i, j]) >= 0.5:
            print(f"  {SHORT[i]:7s} ~ {SHORT[j]:7s} {corr[i, j]:+.2f}")

    print(f"\nLOO-CV RMSE of GP on log(RMSE_all)  [baseline sd of log y = {ly.std():.2f}]")
    for name, labels in BATTERY.items():
        print(f"  {name:38s} {loo_rmse(U, ly, cols_of(labels)):.3f}")

    print("\nDrop-one from all 9 (delta vs full; positive = the parameter carries signal):")
    base = loo_rmse(U, ly, list(range(len(SHORT))))
    for j, s in enumerate(SHORT):
        score = loo_rmse(U, ly, [k for k in range(len(SHORT)) if k != j])
        print(f"  drop {s:7s} {score:.3f}  (delta {score - base:+.3f})")

    basin = y < BASIN_RMSE
    print(f"\nWithin the converged basin (n={basin.sum()} trials with RMSE<{BASIN_RMSE}): "
          "spearman(u, RMSE_all). Weak evidence -- every parameter spans <=20% of its range here.")
    for j, s in enumerate(SHORT):
        print(f"  {s:7s} {st.spearmanr(U[basin, j], y[basin])[0]:+.2f}   "
              f"u-range {U[basin, j].min():.2f}..{U[basin, j].max():.2f}")


if __name__ == "__main__":
    main()
