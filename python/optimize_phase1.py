"""
Phase 1 Bayesian optimization: tune the sacCer3 coupled 4-mark GenomicLayers
rate knobs (R/model5_tf_active_dna_seq.R) to minimize RMSE against the 9
empirical H3K4/H3K9/H3K27 me1/me2/me3 coverage targets.

Each trial shells out to `Rscript R/model5_tf_active_dna_seq.R --key=value ...`,
which runs a full sacCer3 simulation and writes a small JSON objective file to
output/sacCer3_coupled_model/<sim_tag>.objective.json. This script reads that
JSON back as the objective value for scikit-optimize.

Only the SPACE dimensions are searched; FIXED_PARAMS are passed to every trial unchanged.
Outputs are named by --run-tag (output/skopt_<tag>_trials.csv, output/skopt_<tag>_best.json,
sim_tag skopt_<tag>_trial_NNNN) so a new campaign never overwrites an earlier one.

Usage:
    .venv/bin/python python/optimize_phase1.py --n-calls 40 --n-iter 50

Each call runs one full R simulation, so wall-clock time is
roughly n_calls * (time for one Rscript run at the given --n-iter).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from skopt import gp_minimize
from skopt.space import Integer, Real
from skopt.utils import use_named_args

REPO_ROOT = Path(__file__).resolve().parent.parent
R_SCRIPT = REPO_ROOT / "R" / "model5_tf_active_dna_seq.R"
OUTPUT_DIR = REPO_ROOT / "output" / "sacCer3_coupled_model"

# The original 9-dim Phase 1 search space that produced output/skopt_trials.csv (commit ff5179b).
# Not searched any more; kept so skopt_sensitivity.py / skopt_subset_loo.py can read those trials.
PHASE1_SPACE = [
    Integer(3000, 30000, name="meUp_sampler_K4"),
    Integer(3000, 40000, name="meDown_K4"),
    Integer(3000, 30000, name="meUp_sampler_K9"),
    Integer(3000, 30000, name="meDown_K9"),
    Integer(3000, 30000, name="meUp_sampler_K27"),
    Integer(3000, 30000, name="meDown_K27"),
    Integer(500, 20000, name="K27_spread_abu"),
    Real(0.3, 0.95, name="tf_K9_target_frac"),
    Real(0.1, 0.7, name="tf_K27_target_frac"),
]

# ---- Fixed parameters (reduced campaign) --------------------------------------------------
# Flat in the Phase 1 sensitivity analysis, so pinned at their basin values. Edit here to change.
FIXED_PARAMS = {
    "meUp_sampler_K9": 18000,
    "meUp_sampler_K4": 20000,
    "meDown_K27": 23000,
    "meUp_sampler_K27": 13500,
}
# --------------------------------------------------------------------------------------------

# Free parameters searched by skopt.
SPACE = [
    Real(0.15, 0.45, name="tf_K27_target_frac"),
    Real(0.45, 0.80, name="tf_K9_target_frac"),
    Integer(12000, 18000, name="meDown_K9"),
    Integer(20000, 60000, name="meDown_K4"),
    Integer(3000, 10000, name="K27_spread_abu"),
]

assert not set(FIXED_PARAMS) & {dim.name for dim in SPACE}, "a parameter is both fixed and free"


def run_trial(params: dict, sim_tag: str, n_iter: int, seed: int | None) -> dict:
    cmd = ["Rscript", str(R_SCRIPT)]
    for key, value in params.items():
        cmd.append(f"--{key}={value}")
    cmd.append(f"--n_iter={n_iter}")
    cmd.append("--verbose=FALSE")
    cmd.append(f"--sim_tag={sim_tag}")
    if seed is not None:
        cmd.append(f"--seed={seed}")

    result = subprocess.run(
        cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=3600
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Rscript failed for {sim_tag} (exit {result.returncode}).\n"
            f"--- stderr tail ---\n{result.stderr[-3000:]}\n"
            f"--- stdout tail ---\n{result.stdout[-1000:]}"
        )

    obj_path = OUTPUT_DIR / f"{sim_tag}.objective.json"
    with open(obj_path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-calls", type=int, default=20, help="Number of skopt evaluations")
    parser.add_argument("--n-iter", type=int, default=50, help="Inner R simulation loop length per trial")
    parser.add_argument("--seed", type=int, default=None, help="Fixed R RNG seed per trial (default: stochastic)")
    parser.add_argument("--random-state", type=int, default=0, help="skopt random_state for reproducible search")
    parser.add_argument("--objective", choices=["RMSE_all", "RMSE_me3"], default="RMSE_all")
    parser.add_argument("--run-tag", default="reduced", help="Names this campaign's output files")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    trials_path = REPO_ROOT / "output" / f"skopt_{args.run_tag}_trials.csv"
    best_path = REPO_ROOT / "output" / f"skopt_{args.run_tag}_best.json"
    if trials_path.exists():
        parser.error(f"{trials_path} already exists; pick a new --run-tag or move it aside")
    trial_counter = {"n": 0}

    param_names = [dim.name for dim in SPACE]
    fixed_names = list(FIXED_PARAMS)
    with open(trials_path, "w") as f:
        f.write(",".join(param_names + fixed_names) + ",RMSE_all,RMSE_me3,sim_tag,elapsed_s\n")

    @use_named_args(SPACE)
    def objective(**params):
        trial_counter["n"] += 1
        sim_tag = f"skopt_{args.run_tag}_trial_{trial_counter['n']:04d}"
        t0 = time.time()
        result = run_trial({**FIXED_PARAMS, **params}, sim_tag, args.n_iter, args.seed)
        elapsed = time.time() - t0

        with open(trials_path, "a") as f:
            row = [str(params[name]) for name in param_names] + [str(FIXED_PARAMS[name]) for name in fixed_names]
            row += [str(result["RMSE_all"]), str(result["RMSE_me3"]), sim_tag, f"{elapsed:.1f}"]
            f.write(",".join(row) + "\n")

        print(
            f"[trial {trial_counter['n']}/{args.n_calls}] "
            f"RMSE_all={result['RMSE_all']:.3f} RMSE_me3={result['RMSE_me3']:.3f} "
            f"({elapsed:.0f}s) params={params}",
            flush=True,
        )
        return result[args.objective]

    res = gp_minimize(
        objective,
        SPACE,
        n_calls=args.n_calls,
        n_initial_points=min(10, args.n_calls),  # skopt default is 10; lets a short smoke test run
        random_state=args.random_state,
    )

    # skopt returns numpy scalars (np.int64 / np.float64), which json cannot serialise.
    best_params = {name: v.item() if hasattr(v, "item") else v for name, v in zip(param_names, res.x)}
    best = {"best_value": float(res.fun), "objective": args.objective, "params": best_params,
            "fixed_params": FIXED_PARAMS}

    with open(best_path, "w") as f:
        json.dump(best, f, indent=2)

    print("\n=== Best trial ===")
    print(json.dumps(best, indent=2))
    print(f"\nAll trials logged to {trials_path}")
    print(f"Best params saved to {best_path}")


if __name__ == "__main__":
    main()
