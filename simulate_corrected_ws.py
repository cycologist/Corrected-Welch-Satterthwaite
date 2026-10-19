#!/usr/bin/env python3
"""Simulation study for the corrected Welch-Satterthwaite estimator.

This is a Python translation of simulate_corrected_ws.R.

The script simulates K independent variance components,

    S_k^2 = chi-square(nu_k) / nu_k,

and compares the corrected estimator with the target degrees of freedom
K * nu. Run it with:

    python3 simulate_corrected_ws.py
"""

from __future__ import annotations

import math
import random
import statistics
import sys
from pathlib import Path

# Import the function file relative to this script, so the script works
# from any working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from corrected_welch_satterthwaite import cws  # noqa: E402

SEED = 20260715
K_VALUES = (2, 4, 6, 8, 10, 20, 30, 40)
NU_VALUES = (1, 2, 3, 4, 5, 10)
NSIM = 10000


def simulate_one(K: int, nu: int, nsim: int) -> dict[str, float]:
    """Run one simulation cell and return its summary statistics."""
    target = K * nu
    nu_values = [nu] * K
    half_nu = nu / 2.0

    df_est = []
    for _ in range(nsim):
        s2 = [random.gammavariate(half_nu, 2.0) / nu for _ in range(K)]
        df_est.append(cws(s2, nu_values))

    mean_df = statistics.fmean(df_est)
    bias = mean_df - target
    relative_bias = bias / target
    mc_se = statistics.stdev(df_est) / math.sqrt(nsim)
    rmse = math.sqrt(
        math.fsum((value - target) ** 2 for value in df_est) / nsim
    )

    return {
        "K": K,
        "nu": nu,
        "target_df": target,
        "mean_df": mean_df,
        "bias": bias,
        "relative_bias": relative_bias,
        "mc_se": mc_se,
        "rmse": rmse,
    }


def main() -> None:
    random.seed(SEED)

    header = (
        f"{'K':>3} {'nu':>2} {'target_df':>9} {'mean_df':>7} "
        f"{'bias':>8} {'relative_bias':>13} {'mc_se':>7} {'rmse':>7}"
    )
    print(header)

    for K in K_VALUES:
        for nu in NU_VALUES:
            result = simulate_one(K, nu, NSIM)
            print(
                f"{result['K']:3d} {result['nu']:2d} "
                f"{result['target_df']:9.0f} {result['mean_df']:7.3f} "
                f"{result['bias']:8.4f} {result['relative_bias']:13.6f} "
                f"{result['mc_se']:7.5f} {result['rmse']:7.3f}"
            )


if __name__ == "__main__":
    main()
