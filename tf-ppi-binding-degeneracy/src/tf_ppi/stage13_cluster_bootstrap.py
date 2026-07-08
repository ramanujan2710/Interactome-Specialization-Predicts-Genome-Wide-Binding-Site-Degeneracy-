"""
Stage 13 -- TF-family cluster bootstrap: the most conservative robustness
test, resampling whole DNA-binding-domain families (rather than individual
TFs) with replacement, to account for phylogenetic/paralog non-independence.

NOTE: with only ~12 distinguishable DBD family clusters in this dataset,
this test is imperfectly powered for a fully decisive verdict. Report the
resulting confidence interval honestly rather than treating it as settled
(see the manuscript Limitations section).
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .config import out_path
from .utils import rankz, resid


def _bootstrap(df, rng, n_boot):
    fams = df["family"].unique()
    rs = []
    for _ in range(n_boot):
        chosen = rng.choice(fams, size=len(fams), replace=True)
        b = pd.concat([df[df["family"] == f] for f in chosen], ignore_index=True)
        if len(b) < 10 or b["site_count_sample"].std() == 0 or b["phys_module_fraction"].std() == 0:
            continue
        ctrl = np.column_stack([rankz(b["avg_ic"].values), rankz(b["phys_total_degree"].values), rankz(b["annot_len"].values)])
        yr = resid(rankz(b["site_count_sample"].values), ctrl)
        xr = resid(rankz(b["phys_module_fraction"].values), ctrl)
        r, _ = pearsonr(xr, yr)
        if not np.isnan(r):
            rs.append(r)
    return np.array(rs)


def run(n_boot=2000, seed=42):
    print("\n### STAGE 13: TF-family cluster bootstrap (final robustness test) ###")
    m = pd.read_csv(out_path("physical_module_full.tsv"), sep="\t")
    ok = (~m["phys_module_fraction"].isna()) & (~m["site_count_sample"].isna())
    sub = m[ok].copy()
    rng = np.random.default_rng(seed)

    rs = _bootstrap(sub, rng, n_boot)
    lo, hi = np.percentile(rs, [2.5, 97.5])
    print(f"[full cohort] n_boot={len(rs)} mean r={rs.mean():.3f} 95% CI=[{lo:.3f},{hi:.3f}] frac>0={np.mean(rs>0):.3f}")

    non_homeo = sub[sub["family"] != "homeodomain"]
    rs2 = _bootstrap(non_homeo, rng, n_boot)
    lo2, hi2 = np.percentile(rs2, [2.5, 97.5])
    print(f"[non-homeodomain] n_boot={len(rs2)} mean r={rs2.mean():.3f} 95% CI=[{lo2:.3f},{hi2:.3f}] frac>0={np.mean(rs2>0):.3f}")
    print("NOTE: with ~12 distinguishable DBD family clusters, this test is underpowered for a fully")
    print("      decisive verdict -- report the wide CI honestly rather than treating it as settled.")


if __name__ == "__main__":
    run()
