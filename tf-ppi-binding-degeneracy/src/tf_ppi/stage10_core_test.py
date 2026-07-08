"""
Stage 10 -- CORE TEST: physical-evidence-only PPI module fraction vs.
genome-wide predicted binding-site density. This is the headline result of
the study.

Output: physical_module_full.tsv
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .config import out_path
from .utils import rankz, resid


def run():
    print("\n### STAGE 10: physical-evidence module fraction vs. site density (CORE TEST) ###")
    d = np.load(out_path("physical_edges.npz"))
    i1, i2 = d["i1"], d["i2"]

    flag = pd.read_csv(out_path("protein_tf_flag.tsv"), sep="\t")
    n_nodes = len(flag)
    all_ids = flag["protein_id"].tolist()

    net = pd.read_csv(out_path("compartment_compare_full.tsv"), sep="\t")
    tf_pids = set(net["protein_id"])
    tf_mask = np.array([pid in tf_pids for pid in all_ids])
    pid_to_row = {pid: idx for idx, pid in enumerate(net["protein_id"])}

    deg_to_tf = np.zeros(n_nodes)
    total_deg = np.zeros(n_nodes)
    np.add.at(total_deg, i1, 1)
    np.add.at(total_deg, i2, 1)
    mask_i2_tf, mask_i1_tf = tf_mask[i2], tf_mask[i1]
    np.add.at(deg_to_tf, i1, mask_i2_tf.astype(float))
    np.add.at(deg_to_tf, i2, mask_i1_tf.astype(float))

    rows = []
    for idx, pid in enumerate(all_ids):
        if pid in pid_to_row:
            td = total_deg[idx]
            rows.append({"protein_id": pid, "phys_total_degree": td, "phys_deg_to_tf": deg_to_tf[idx],
                         "phys_module_fraction": (deg_to_tf[idx] / td) if td > 0 else np.nan})
    pdf = pd.DataFrame(rows)
    m = net.merge(pdf, on="protein_id", how="inner")
    print("n=", len(m), " (nonzero phys degree:", (m["phys_total_degree"] > 0).sum(), ")")

    r1, p1 = spearmanr(m["annot_len"], m["phys_total_degree"])
    print(f"[literature-bias check] annot_len vs phys_total_degree: r={r1:.3f} p={p1:.3g} "
          "(should be weaker than the combined-score-network version)")
    r1b, p1b = spearmanr(m["annot_len"], m["total_degree"])
    print(f"[for comparison] annot_len vs combined-score total_degree: r={r1b:.3f} p={p1b:.3g}")

    r2, p2 = spearmanr(m["phys_module_fraction"], m["tf_module_fraction"])
    print(f"[validity check] phys_module_fraction vs combined-score module_fraction: r={r2:.3f} p={p2:.3g}")

    ok = (~m["phys_module_fraction"].isna()) & (~m["site_count_sample"].isna())
    sub = m[ok]
    r, p = pearsonr(sub["phys_module_fraction"], sub["site_count_sample"])
    print(f"\n[HEADLINE, full cohort] raw: r={r:.3f} p={p:.3g} n={len(sub)}")
    ctrl = np.column_stack([rankz(sub["avg_ic"].values), rankz(sub["phys_total_degree"].values), rankz(sub["annot_len"].values)])
    yr = resid(rankz(sub["site_count_sample"].values), ctrl)
    xr = resid(rankz(sub["phys_module_fraction"].values), ctrl)
    rp, pp = pearsonr(xr, yr)
    print(f"[HEADLINE, full cohort] partial (avg_ic, phys_degree, annot_len): r={rp:.3f} p={pp:.3g}")

    non_homeo = sub[sub["family"] != "homeodomain"]
    r3, p3 = pearsonr(non_homeo["phys_module_fraction"], non_homeo["site_count_sample"])
    ctrl3 = np.column_stack([rankz(non_homeo["avg_ic"].values), rankz(non_homeo["phys_total_degree"].values), rankz(non_homeo["annot_len"].values)])
    yr3 = resid(rankz(non_homeo["site_count_sample"].values), ctrl3)
    xr3 = resid(rankz(non_homeo["phys_module_fraction"].values), ctrl3)
    rp3, pp3 = pearsonr(xr3, yr3)
    print(f"[non-homeodomain] raw r={r3:.3f} p={p3:.3g} | partial r={rp3:.3f} p={pp3:.3g} n={len(non_homeo)}")

    m.to_csv(out_path("physical_module_full.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    run()
