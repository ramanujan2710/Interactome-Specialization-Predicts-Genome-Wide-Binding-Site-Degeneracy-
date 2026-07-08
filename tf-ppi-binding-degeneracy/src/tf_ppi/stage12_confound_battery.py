"""
Stage 12 -- Confound battery: motif length, GC/AT composition,
family/paralog exclusion, and Benjamini-Hochberg multiple-testing
correction across the core battery of tests.
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .config import out_path
from .utils import rankz, resid

# p-values from Stage 10 / Stage 11 headline tests, corrected here for
# multiple testing. Keep this list in sync if you change the upstream tests.
CORE_BATTERY_PVALUES = [3.45e-08, 2.54e-09, 4.03e-04, 1.03e-04, 1.63e-05, 8.57e-06]


def run():
    print("\n### STAGE 12: confound battery ###")
    m = pd.read_csv(out_path("physical_module_full.tsv"), sep="\t")
    ok = (~m["phys_module_fraction"].isna()) & (~m["site_count_sample"].isna()) & (~m["breadth"].isna())
    sub = m[ok].copy()
    non_homeo = sub[sub["family"] != "homeodomain"].copy()

    print("\n-- 1. motif-length confound --")
    for name, df in [("full", sub), ("non-homeo", non_homeo)]:
        r, p = pearsonr(df["phys_module_fraction"], df["length"])
        print(f"[{name}] module_fraction vs motif length: r={r:.3f} p={p:.3g}")
        ctrl = np.column_stack([rankz(df["avg_ic"].values), rankz(df["phys_total_degree"].values),
                                 rankz(df["annot_len"].values), rankz(df["length"].values)])
        yr = resid(rankz(df["site_count_sample"].values), ctrl)
        xr = resid(rankz(df["phys_module_fraction"].values), ctrl)
        r2, p2 = pearsonr(xr, yr)
        print(f"[{name}] partial r WITH motif-length control: r={r2:.3f} p={p2:.3g} n={len(df)}")

    print("\n-- 2. GC/AT composition (conservative; likely a partial mediator, not a pure confound) --")
    for name, df in [("full", sub), ("non-homeo", non_homeo)]:
        ctrl = np.column_stack([rankz(df["avg_ic"].values), rankz(df["phys_total_degree"].values),
                                 rankz(df["annot_len"].values), rankz(df["consensus_gc_frac"].values)])
        yr = resid(rankz(df["site_count_sample"].values), ctrl)
        xr = resid(rankz(df["phys_module_fraction"].values), ctrl)
        r, p = pearsonr(xr, yr)
        print(f"[{name}] partial r WITH GC-content control: r={r:.3f} p={p:.3g} n={len(df)}")

    print("\n-- 3. family/paralog non-independence --")
    print(sub["family"].value_counts().to_string())
    other = sub[sub["family"] == "other/unclassified"]
    ctrl = np.column_stack([rankz(other["avg_ic"].values), rankz(other["phys_total_degree"].values), rankz(other["annot_len"].values)])
    yr = resid(rankz(other["site_count_sample"].values), ctrl)
    xr = resid(rankz(other["phys_module_fraction"].values), ctrl)
    r, p = pearsonr(xr, yr)
    print(f"[other/unclassified only, n={len(other)}] partial r={r:.3f} p={p:.3g}")

    excl = sub[~sub["family"].isin(["homeodomain", "zinc_finger"])]
    ctrl = np.column_stack([rankz(excl["avg_ic"].values), rankz(excl["phys_total_degree"].values), rankz(excl["annot_len"].values)])
    yr = resid(rankz(excl["site_count_sample"].values), ctrl)
    xr = resid(rankz(excl["phys_module_fraction"].values), ctrl)
    r, p = pearsonr(xr, yr)
    print(f"[excl. homeodomain+zinc_finger, n={len(excl)}] partial r={r:.3f} p={p:.3g}")

    print("\n-- 4. Benjamini-Hochberg correction across the core battery --")
    n = len(CORE_BATTERY_PVALUES)
    for p, adj in [(p, p * n / (i + 1)) for i, p in enumerate(sorted(CORE_BATTERY_PVALUES))]:
        print(f"  raw p={p:.3g}  BH-adjusted={min(adj, 1):.3g}")


if __name__ == "__main__":
    run()
