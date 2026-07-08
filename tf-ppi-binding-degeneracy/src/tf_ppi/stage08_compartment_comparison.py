"""
Stage 8 -- Genic vs. intergenic compartment comparison; builds the main
"everything merged" table (compartment_compare_full.tsv) used by all
downstream tests.

Output: compartment_compare_full.tsv
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .config import out_path
from .utils import rankz, resid, fisher_z_test


def _load_density(state_file):
    state = json.load(open(state_file))
    total_bp = state["total_bp_done"]
    rows = [{"name_lc": name, "site_count": cnt, "density_per_kb": cnt / (total_bp / 1000.0)}
            for name, cnt in state["site_count"].items()]
    return pd.DataFrame(rows), total_bp


def run():
    print("\n### STAGE 8: genic vs. intergenic compartment comparison ###")
    genic_df, genic_bp = _load_density(out_path("gene_sample_agg_state.json"))
    inter_df, inter_bp = _load_density(out_path("intergenic_sample_agg_state.json"))
    print("genic bp:", genic_bp, "intergenic bp:", inter_bp)

    m = pd.read_csv(out_path("tf_with_family_full.tsv"), sep="\t")
    d = np.load(out_path("coexpr_matrix.npz"), allow_pickle=True)
    bdf = pd.DataFrame({"name_lc": d["names"], "breadth": d["breadth"]})
    m = m.merge(bdf, on="name_lc", how="left")
    m = m.merge(genic_df.rename(columns={"site_count": "site_count_genic", "density_per_kb": "density_genic"}),
                on="name_lc", how="inner")
    m = m.merge(inter_df.rename(columns={"site_count": "site_count_inter", "density_per_kb": "density_inter"}),
                on="name_lc", how="inner")
    # site_count_sample = the gene-body-restricted sample, used as the
    # primary "genome-wide predicted site density" measure downstream.
    m["site_count_sample"] = m["site_count_genic"]
    print("merged n=", len(m))

    print("\n=== module_fraction vs density, GENIC vs INTERGENIC ===")
    for label, col in [("GENIC", "site_count_genic"), ("INTERGENIC", "site_count_inter")]:
        r, p = pearsonr(m["tf_module_fraction"], m[col])
        ctrl = np.column_stack([rankz(m["avg_ic"].values), rankz(m["total_degree"].values), rankz(m["annot_len"].values)])
        yr = resid(rankz(m[col].values), ctrl)
        xr = resid(rankz(m["tf_module_fraction"].values), ctrl)
        rp, pp = pearsonr(xr, yr)
        print(f"{label}: raw r={r:.3f} p={p:.3g} | partial r={rp:.3f} p={pp:.3g} n={len(m)}")

    non_homeo = m[m["family"] != "homeodomain"]
    r_genic, _ = pearsonr(non_homeo["tf_module_fraction"], non_homeo["site_count_genic"])
    r_inter, _ = pearsonr(non_homeo["tf_module_fraction"], non_homeo["site_count_inter"])
    z, p = fisher_z_test(r_inter, len(non_homeo), r_genic, len(non_homeo))
    print(f"\nFisher z-test (non-homeodomain): intergenic r={r_inter:.3f} vs genic r={r_genic:.3f} -> z={z:.2f} p={p:.3g}")

    m.to_csv(out_path("compartment_compare_full.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    run()
