"""
Stage 3 -- Match JASPAR motifs to STRING protein IDs by preferred name; build
the confirmed-TF table with combined-score PPI degree, degree-to-TF-module,
and the (combined-score) "module fraction". This combined-score network is
used only for sanity/validity checks later; the headline results use the
evidence-purified physical network built in Stage 9.

Outputs: arrays.npz, merged_full.tsv, flag_with_jaspar_confirmed.tsv
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .config import PATHS, out_path


def run():
    print("\n### STAGE 3: JASPAR<->STRING merge + combined-score network stats ###")
    jaspar = json.load(open(out_path("jaspar_ic.json")))
    info = pd.read_csv(PATHS["info"], sep="\t")
    info["name_lc"] = info["preferred_name"].str.lower()
    name_to_pid = dict(zip(info["name_lc"], info["#string_protein_id"]))

    rows = []
    for m in jaspar:
        key = m["name"].lower()
        if key in name_to_pid:
            rows.append({"ma_id": m["ma_id"], "name": m["name"],
                         "protein_id": name_to_pid[key], "length": m["length"],
                         "total_ic": m["total_ic"], "avg_ic": m["avg_ic"],
                         "consensus_gc_frac": m["consensus_gc_frac"]})
    jdf = pd.DataFrame(rows)
    print("JASPAR motifs matched to STRING IDs:", len(jdf),
          "unique proteins:", jdf["protein_id"].nunique())

    jic = jdf.groupby("protein_id").agg(
        avg_ic=("avg_ic", "mean"), total_ic=("total_ic", "mean"),
        length=("length", "mean"), n_motifs=("ma_id", "count"),
        consensus_gc_frac=("consensus_gc_frac", "mean")).reset_index()

    flag = pd.read_csv(out_path("protein_tf_flag.tsv"), sep="\t")
    flag["protein_id"] = flag["protein_id"].astype(str)
    all_ids = flag["protein_id"].tolist()
    idx_map = {pid: i for i, pid in enumerate(all_ids)}
    n_nodes = len(all_ids)

    # --- build deduplicated (i<j) combined-score network as integer arrays ---
    i_list, j_list, w_list = [], [], []
    with open(PATHS["links"]) as f:
        next(f)
        for line in f:
            p1, p2, score = line.split()
            if p1 not in idx_map or p2 not in idx_map:
                continue
            i1, i2 = idx_map[p1], idx_map[p2]
            if i1 < i2:
                i_list.append(i1); j_list.append(i2); w_list.append(int(score))
    i_arr = np.array(i_list, dtype=np.int32)
    j_arr = np.array(j_list, dtype=np.int32)
    sc_arr = np.array(w_list, dtype=np.int32)
    np.savez(out_path("arrays.npz"), i1=i_arr, i2=j_arr, sc=sc_arr, n_nodes=n_nodes)
    print("combined-score network: nodes =", n_nodes, "deduped edges =", len(i_arr))

    jaspar_tf_pids = set(jdf["protein_id"].unique())
    jaspar_mask = np.array([pid in jaspar_tf_pids for pid in all_ids])

    deg_to_tf = np.zeros(n_nodes)
    total_deg = np.zeros(n_nodes)
    mask_i2_tf = jaspar_mask[j_arr]
    mask_i1_tf = jaspar_mask[i_arr]
    np.add.at(deg_to_tf, i_arr, mask_i2_tf.astype(float))
    np.add.at(deg_to_tf, j_arr, mask_i1_tf.astype(float))
    np.add.at(total_deg, i_arr, 1); np.add.at(total_deg, j_arr, 1)

    flag["deg_to_jaspar_tf"] = deg_to_tf
    flag["total_degree"] = total_deg
    flag["jaspar_tf"] = jaspar_mask.astype(int)

    merged = flag.merge(jic, on="protein_id", how="inner")
    merged["tf_module_fraction"] = merged["deg_to_jaspar_tf"] / merged["total_degree"].replace(0, np.nan)
    merged["name_lc"] = merged["name"].str.lower()

    r, p = spearmanr(merged["avg_ic"], merged["tf_module_fraction"])
    print(f"[sanity] avg_ic vs combined-score module_fraction: r={r:.3f} p={p:.3g} "
          "(reported only as a sanity check, NOT a headline result)")

    merged.to_csv(out_path("merged_full.tsv"), sep="\t", index=False)
    flag.to_csv(out_path("flag_with_jaspar_confirmed.tsv"), sep="\t", index=False)
    print("confirmed TFs with PWM + network data:", len(merged))


if __name__ == "__main__":
    run()
