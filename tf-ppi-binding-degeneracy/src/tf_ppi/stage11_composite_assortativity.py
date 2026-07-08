"""
Stage 11 -- Composite specialization index (module fraction + inverse
expression breadth), and TF-TF assortativity permutation tests (naive and
degree-matched), both computed on the physical-evidence-only network.

No file outputs; prints results to stdout (consumed by figures/make_figures.py
for plotting).
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .config import out_path
from .utils import rankz, resid


def _mask(rng, idx_all, n_nodes, n_tf):
    chosen = rng.choice(idx_all, size=n_tf, replace=False)
    mask = np.zeros(n_nodes, dtype=bool)
    mask[chosen] = True
    return mask


def run():
    print("\n### STAGE 11: composite index + physical-network TF-TF assortativity ###")
    m = pd.read_csv(out_path("physical_module_full.tsv"), sep="\t")

    non_homeo = m[(m["family"] != "homeodomain") & (~m["breadth"].isna())]
    comp = (rankz(non_homeo["phys_module_fraction"].values) + rankz(-non_homeo["breadth"].values)) / 2
    y = rankz(non_homeo["site_count_sample"].values)
    r, p = pearsonr(comp, y)
    print(f"[composite, physical-only] raw: r={r:.3f} p={p:.3g} n={len(non_homeo)}")
    ctrl = np.column_stack([rankz(non_homeo["avg_ic"].values), rankz(non_homeo["phys_total_degree"].values), rankz(non_homeo["annot_len"].values)])
    yr = resid(y, ctrl)
    cr = resid(comp, ctrl)
    r2, p2 = pearsonr(cr, yr)
    print(f"[composite, physical-only] controlled: r={r2:.3f} p={p2:.3g}")

    print("\n-- TF-TF assortativity, physical-evidence network --")
    de = np.load(out_path("physical_edges.npz"))
    i1p, i2p = de["i1"], de["i2"]
    flag = pd.read_csv(out_path("protein_tf_flag.tsv"), sep="\t")
    n_nodes = len(flag)
    tf_pids = set(m["protein_id"])
    all_ids = flag["protein_id"].tolist()
    tf_mask = np.array([pid in tf_pids for pid in all_ids])
    n_tf = tf_mask.sum()
    print("n_tf:", n_tf, "n_nodes:", n_nodes)

    def stat(mask, i_, j_):
        return int((mask[i_] & mask[j_]).sum())

    obs = stat(tf_mask, i1p, i2p)
    print("observed TF-TF physical edges:", obs)

    rng = np.random.default_rng(7)
    idx_all = np.arange(n_nodes)
    n_perm = 500
    null_counts = np.array([stat(_mask(rng, idx_all, n_nodes, n_tf), i1p, i2p) for _ in range(n_perm)])
    z = (obs - null_counts.mean()) / null_counts.std()
    pval = (np.sum(null_counts >= obs) + 1) / (n_perm + 1)
    print(f"NAIVE PERM null: mean={null_counts.mean():.1f} sd={null_counts.std():.1f} obs={obs} z={z:.2f} p={pval:.4f}")

    deg = np.zeros(n_nodes)
    np.add.at(deg, i1p, 1)
    np.add.at(deg, i2p, 1)
    for n_bins in [5, 10, 20, 50]:
        order = np.argsort(deg)
        rank = np.empty(n_nodes, dtype=int)
        rank[order] = np.arange(n_nodes)
        bin_id = rank * n_bins // n_nodes
        tf_bin_counts = {b: int(np.sum(tf_mask & (bin_id == b))) for b in np.unique(bin_id)}
        bin_pools = {b: idx_all[bin_id == b] for b in np.unique(bin_id)}
        null2 = np.zeros(n_perm)
        for pidx in range(n_perm):
            mask = np.zeros(n_nodes, dtype=bool)
            for b, k in tf_bin_counts.items():
                if k == 0:
                    continue
                mask[rng.choice(bin_pools[b], size=k, replace=False)] = True
            null2[pidx] = stat(mask, i1p, i2p)
        z2 = (obs - null2.mean()) / null2.std()
        print(f"DEGREE-MATCHED PERM (n_bins={n_bins}): mean={null2.mean():.1f} sd={null2.std():.1f} obs={obs} z={z2:.2f}")


if __name__ == "__main__":
    run()
