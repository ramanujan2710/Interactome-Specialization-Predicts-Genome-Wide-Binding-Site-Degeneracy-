"""
Stage 14 -- Adversarial artifact checks: is the core result tautological
with motif information content, driven by a handful of extreme
observations, or an artifact of Pearson's sensitivity to skew?

These checks are what caught an invalid score-threshold bug during
development of this pipeline (see config.py, SCORE_THRESHOLD) -- keep them
as a standing regression test whenever the upstream data or methodology
changes.
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .config import out_path
from .utils import rankz, resid, winsorize, post_hoc_power


def run():
    print("\n### STAGE 14: adversarial artifact / outlier / tautology checks ###")
    m = pd.read_csv(out_path("physical_module_full.tsv"), sep="\t")
    ok = (~m["phys_module_fraction"].isna()) & (~m["site_count_sample"].isna())
    sub = m[ok].copy()

    print("\n-- 1. tautology check: is module_fraction just avg_ic in disguise? --")
    r, p = pearsonr(sub["avg_ic"], sub["site_count_sample"])
    print(f"avg_ic vs site_count (expected weak/nonsignificant if not tautological): r={r:.3f} p={p:.3g}")
    r2, p2 = pearsonr(sub["phys_module_fraction"], sub["avg_ic"])
    print(f"phys_module_fraction vs avg_ic (should be near-zero, i.e. not collinear): r={r2:.3f} p={p2:.3g}")

    print("\n-- 2. Spearman vs Pearson (rank-robustness) --")
    for name, df in [("full", sub), ("non-homeo", sub[sub["family"] != "homeodomain"])]:
        rp, pp = pearsonr(df["phys_module_fraction"], df["site_count_sample"])
        rs, ps = spearmanr(df["phys_module_fraction"], df["site_count_sample"])
        print(f"[{name}] Pearson r={rp:.3f} p={pp:.3g} | Spearman rho={rs:.3f} p={ps:.3g} n={len(df)}")

    print("\n-- 3. leverage/outlier sensitivity: drop top-k most extreme site_count values --")
    for k in [5, 10, 20]:
        trimmed = sub.sort_values("site_count_sample", ascending=False).iloc[k:]
        r, p = pearsonr(trimmed["phys_module_fraction"], trimmed["site_count_sample"])
        ctrl = np.column_stack([rankz(trimmed["avg_ic"].values), rankz(trimmed["phys_total_degree"].values), rankz(trimmed["annot_len"].values)])
        yr = resid(rankz(trimmed["site_count_sample"].values), ctrl)
        xr = resid(rankz(trimmed["phys_module_fraction"].values), ctrl)
        rp, pp = pearsonr(xr, yr)
        print(f"drop top-{k} highest-density TFs: raw r={r:.3f} p={p:.3g} | partial r={rp:.3f} p={pp:.3g} n={len(trimmed)}")

    print("\n-- 4. winsorized (5%/95%) version --")
    df = sub.copy()
    df["site_count_w"] = winsorize(df["site_count_sample"])
    df["modfrac_w"] = winsorize(df["phys_module_fraction"])
    r, p = pearsonr(df["modfrac_w"], df["site_count_w"])
    ctrl = np.column_stack([rankz(df["avg_ic"].values), rankz(df["phys_total_degree"].values), rankz(df["annot_len"].values)])
    yr = resid(rankz(df["site_count_w"].values), ctrl)
    xr = resid(rankz(df["modfrac_w"].values), ctrl)
    rp, pp = pearsonr(xr, yr)
    print(f"winsorized: raw r={r:.3f} p={p:.3g} | partial r={rp:.3f} p={pp:.3g} n={len(df)}")

    print("\n-- 5. post-hoc power at observed effect size --")
    r_obs, _ = pearsonr(sub["phys_module_fraction"], sub["site_count_sample"])
    print(f"observed r={r_obs:.3f}, n={len(sub)} -> post-hoc power = {post_hoc_power(r_obs, len(sub)):.3f}")

    print("\nIf ANY of the above look off (near-zero rank correlation, effect driven entirely by")
    print("dropping outliers, power far below 0.8), stop and re-examine before trusting Stage 10-13.")


if __name__ == "__main__":
    run()
