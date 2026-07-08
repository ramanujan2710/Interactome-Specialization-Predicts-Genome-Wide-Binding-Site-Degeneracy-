#!/usr/bin/env python3
"""
Generate the five manuscript figures directly from pipeline output tables.

Run this AFTER completing all 14 pipeline stages (run_pipeline.py). By
default it reads from ../results and writes PNGs into this directory.

    python figures/make_figures.py
    python figures/make_figures.py --results-dir /path/to/results --out-dir /path/to/figures
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def rankz(x):
    r = rankdata(x)
    return (r - r.mean()) / r.std()


def resid(y, X):
    A = np.column_stack([X, np.ones(len(y))])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ b


FAMILY_COLORS = {
    "homeodomain": "#d62728", "zinc_finger": "#1f77b4", "bHLH": "#2ca02c",
    "bZIP": "#9467bd", "nuclear_receptor": "#ff7f0e", "forkhead": "#17becf",
    "other/unclassified": "#7f7f7f", "multiple": "#8c564b", "T_box": "#e377c2",
    "Pax_paired": "#bcbd22", "POU": "#aec7e8", "HMG_box": "#c49c94",
}


def colors_for(df):
    return [FAMILY_COLORS.get(f, "#7f7f7f") for f in df["family"]]


def bootstrap(df, rng, n_boot=1500):
    fams = df["family"].unique()
    rs = []
    for _ in range(n_boot):
        chosen = rng.choice(fams, size=len(fams), replace=True)
        b = pd.concat([df[df["family"] == f] for f in chosen], ignore_index=True)
        if len(b) < 10 or b["site_count_sample"].std() == 0 or b["phys_module_fraction"].std() == 0:
            continue
        c = np.column_stack([rankz(b["avg_ic"].values), rankz(b["phys_total_degree"].values), rankz(b["annot_len"].values)])
        yb = resid(rankz(b["site_count_sample"].values), c)
        xb = resid(rankz(b["phys_module_fraction"].values), c)
        rr, _ = pearsonr(xb, yb)
        if not np.isnan(rr):
            rs.append(rr)
    return np.array(rs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="../results")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()
    RES, OUTDIR = args.results_dir, args.out_dir
    os.makedirs(OUTDIR, exist_ok=True)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "figure.dpi": 200, "savefig.dpi": 300,
    })

    m = pd.read_csv(f"{RES}/physical_module_full.tsv", sep="\t")
    ok = (~m["phys_module_fraction"].isna()) & (~m["site_count_sample"].isna())
    sub = m[ok].copy()
    non_homeo = sub[sub["family"] != "homeodomain"].copy()

    # ---------------- Figure 1: core scatter ----------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, df, title in [(axes[0], sub, f"All confirmed TFs (n={len(sub)})"),
                           (axes[1], non_homeo, f"Excluding homeodomain family (n={len(non_homeo)})")]:
        x = df["phys_module_fraction"].values
        y = df["site_count_sample"].values
        r, p = pearsonr(x, y)
        ax.scatter(x, y, c=colors_for(df), s=22, alpha=0.75, edgecolor="white", linewidth=0.3)
        b = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, np.polyval(b, xs), color="black", lw=1.6, ls="--")
        ax.set_xlabel("Physical-evidence PPI module fraction\n(share of physical interactome that is other TFs)")
        ax.set_ylabel("Genome-wide predicted site count\n(18 Mb gene-body sample)")
        ax.set_title(f"{title}\nPearson r={r:.3f}, p={p:.2g}", fontsize=10)
    handles = [mpatches.Patch(color=c, label=f.replace("_", " ")) for f, c in FAMILY_COLORS.items() if f in sub["family"].unique()]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Figure 1. PPI interactome specialization predicts genome-wide binding-site density", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figure1_core_scatter.png", bbox_inches="tight")
    plt.close(fig)

    # ---------------- Figure 2: confound forest plot ----------------
    rows = []

    def add(label, r, p, n, lo=None, hi=None):
        rows.append({"label": label, "r": r, "p": p, "n": n, "lo": lo, "hi": hi})

    ctrl = np.column_stack([rankz(sub["avg_ic"].values), rankz(sub["phys_total_degree"].values), rankz(sub["annot_len"].values)])
    yr = resid(rankz(sub["site_count_sample"].values), ctrl); xr = resid(rankz(sub["phys_module_fraction"].values), ctrl)
    r_p, p_p = pearsonr(xr, yr)
    r_raw, p_raw = pearsonr(sub["phys_module_fraction"], sub["site_count_sample"])
    add("Raw correlation (full cohort)", r_raw, p_raw, len(sub))
    add("Partial: + info content, degree, annotation length", r_p, p_p, len(sub))

    ctrl_l = np.column_stack([rankz(sub["avg_ic"].values), rankz(sub["phys_total_degree"].values),
                               rankz(sub["annot_len"].values), rankz(sub["length"].values)])
    yrl = resid(rankz(sub["site_count_sample"].values), ctrl_l); xrl = resid(rankz(sub["phys_module_fraction"].values), ctrl_l)
    r_l, p_l = pearsonr(xrl, yrl)
    add("+ motif length control", r_l, p_l, len(sub))

    ctrl_g = np.column_stack([rankz(sub["avg_ic"].values), rankz(sub["phys_total_degree"].values),
                               rankz(sub["annot_len"].values), rankz(sub["consensus_gc_frac"].values)])
    yrg = resid(rankz(sub["site_count_sample"].values), ctrl_g); xrg = resid(rankz(sub["phys_module_fraction"].values), ctrl_g)
    r_g, p_g = pearsonr(xrg, yrg)
    add("+ GC/AT composition control (conservative)", r_g, p_g, len(sub))

    ctrl_nh = np.column_stack([rankz(non_homeo["avg_ic"].values), rankz(non_homeo["phys_total_degree"].values), rankz(non_homeo["annot_len"].values)])
    yr_nh = resid(rankz(non_homeo["site_count_sample"].values), ctrl_nh); xr_nh = resid(rankz(non_homeo["phys_module_fraction"].values), ctrl_nh)
    r_nh, p_nh = pearsonr(xr_nh, yr_nh)
    add("Non-homeodomain stratum, partial", r_nh, p_nh, len(non_homeo))

    other = sub[sub["family"] == "other/unclassified"]
    ctrl_o = np.column_stack([rankz(other["avg_ic"].values), rankz(other["phys_total_degree"].values), rankz(other["annot_len"].values)])
    yro = resid(rankz(other["site_count_sample"].values), ctrl_o); xro = resid(rankz(other["phys_module_fraction"].values), ctrl_o)
    r_o, p_o = pearsonr(xro, yro)
    add("Unclassified-family-only subset, partial", r_o, p_o, len(other))

    excl = sub[~sub["family"].isin(["homeodomain", "zinc_finger"])]
    ctrl_e = np.column_stack([rankz(excl["avg_ic"].values), rankz(excl["phys_total_degree"].values), rankz(excl["annot_len"].values)])
    yre = resid(rankz(excl["site_count_sample"].values), ctrl_e); xre = resid(rankz(excl["phys_module_fraction"].values), ctrl_e)
    r_e, p_e = pearsonr(xre, yre)
    add("Excl. homeodomain + zinc-finger, partial", r_e, p_e, len(excl))

    rng = np.random.default_rng(42)
    rs_boot = bootstrap(sub, rng)
    lo, hi = np.percentile(rs_boot, [2.5, 97.5])
    add("Family-cluster bootstrap mean (95% CI)", rs_boot.mean(), np.nan, len(sub), lo, hi)

    fdf = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 5))
    ypos = np.arange(len(fdf))[::-1]
    for i, row in zip(ypos, fdf.itertuples()):
        color = "#1f77b4" if row.r > 0 else "#d62728"
        ax.plot(row.r, i, "o", color=color, ms=7, zorder=3)
        if row.lo is not None and not np.isnan(row.lo):
            ax.plot([row.lo, row.hi], [i, i], color=color, lw=2, alpha=0.5, zorder=2)
    ax.axvline(0, color="grey", lw=1, ls=":")
    ax.set_yticks(ypos); ax.set_yticklabels(fdf["label"])
    ax.set_xlabel("Effect size (Pearson r, rank-transformed where partialled)")
    ax.set_title("Figure 2. Effect size across the full confound and stratification battery", fontsize=11)
    for i, row in zip(ypos, fdf.itertuples()):
        txt = f"r={row.r:.2f}" + (f", p={row.p:.1e}" if not np.isnan(row.p) else " (bootstrap mean)")
        ax.text(0.62, i, txt, fontsize=7.5, va="center")
    ax.set_xlim(-0.15, 0.85)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figure2_confound_forest.png", bbox_inches="tight")
    plt.close(fig)

    # ---------------- Figure 3: composite index + compartment replication ----------------
    comp_df = pd.read_csv(f"{RES}/compartment_compare_full.tsv", sep="\t")
    modfrac_col = "phys_module_fraction" if "phys_module_fraction" in comp_df.columns else "tf_module_fraction"
    non_homeo_c = comp_df[(comp_df["family"] != "homeodomain") & (~comp_df["breadth"].isna())].copy()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    comp = (rankz(non_homeo_c[modfrac_col].values) + rankz(-non_homeo_c["breadth"].values)) / 2
    y = rankz(non_homeo_c["site_count_sample"].values)
    r, p = pearsonr(comp, y)
    ax.scatter(comp, y, s=20, alpha=0.7, color="#2ca02c", edgecolor="white", linewidth=0.3)
    b = np.polyfit(comp, y, 1); xs = np.linspace(comp.min(), comp.max(), 50)
    ax.plot(xs, np.polyval(b, xs), color="black", lw=1.6, ls="--")
    ax.set_xlabel("Composite specialization index\n(rank-z module fraction + rank-z inverse expression breadth)")
    ax.set_ylabel("Rank-z genome-wide site density")
    ax.set_title(f"Composite index, non-homeodomain (n={len(non_homeo_c)})\nr={r:.3f}, p={p:.2g}", fontsize=10)

    ax2 = axes[1]
    genic_r, _ = pearsonr(non_homeo_c[modfrac_col], non_homeo_c["site_count_genic"])
    inter_r, _ = pearsonr(non_homeo_c[modfrac_col], non_homeo_c["site_count_inter"])
    bars = ax2.bar(["Gene-body\ncompartment", "Intergenic\ncompartment"], [genic_r, inter_r],
                    color=["#1f77b4", "#ff7f0e"], width=0.5)
    for b_, val in zip(bars, [genic_r, inter_r]):
        ax2.text(b_.get_x() + b_.get_width() / 2, val + 0.01, f"r={val:.3f}", ha="center", fontsize=9)
    ax2.set_ylabel("Module fraction vs. site density (r)")
    ax2.set_title("Independent replication across\ngenomic compartments (Fisher z n.s.)", fontsize=10)
    ax2.set_ylim(0, max(genic_r, inter_r) * 1.3)

    fig.suptitle("Figure 3. Composite specialization index and cross-compartment replication", fontsize=11, y=1.03)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figure3_composite_compartment.png", bbox_inches="tight")
    plt.close(fig)

    # ---------------- Figure 4: TF-TF assortativity ----------------
    de = np.load(f"{RES}/physical_edges.npz")
    i1p, i2p = de["i1"], de["i2"]
    flag = pd.read_csv(f"{RES}/protein_tf_flag.tsv", sep="\t")
    n_nodes = len(flag)
    tf_pids = set(sub["protein_id"])
    all_ids = flag["protein_id"].tolist()
    tf_mask = np.array([pid in tf_pids for pid in all_ids])
    n_tf = tf_mask.sum()

    def stat(mask, i_, j_):
        return int((mask[i_] & mask[j_]).sum())

    obs = stat(tf_mask, i1p, i2p)
    rng2 = np.random.default_rng(7)
    idx_all = np.arange(n_nodes)
    n_perm = 1000
    naive_null = np.zeros(n_perm)
    for pi in range(n_perm):
        chosen = rng2.choice(idx_all, size=n_tf, replace=False)
        mask = np.zeros(n_nodes, dtype=bool); mask[chosen] = True
        naive_null[pi] = stat(mask, i1p, i2p)

    deg = np.zeros(n_nodes)
    np.add.at(deg, i1p, 1); np.add.at(deg, i2p, 1)
    n_bins = 20
    order = np.argsort(deg); rank = np.empty(n_nodes, dtype=int); rank[order] = np.arange(n_nodes)
    bin_id = rank * n_bins // n_nodes
    tf_bin_counts = {b: int(np.sum(tf_mask & (bin_id == b))) for b in np.unique(bin_id)}
    bin_pools = {b: idx_all[bin_id == b] for b in np.unique(bin_id)}
    matched_null = np.zeros(n_perm)
    for pi in range(n_perm):
        mask = np.zeros(n_nodes, dtype=bool)
        for b, k in tf_bin_counts.items():
            if k == 0:
                continue
            mask[rng2.choice(bin_pools[b], size=k, replace=False)] = True
        matched_null[pi] = stat(mask, i1p, i2p)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, null, ttl, z in [
        (axes[0], naive_null, "Naive random-relabeling null", (obs - naive_null.mean()) / naive_null.std()),
        (axes[1], matched_null, "Degree-matched (20-bin) null", (obs - matched_null.mean()) / matched_null.std()),
    ]:
        ax.hist(null, bins=30, color="#7f7f7f", alpha=0.8)
        ax.axvline(obs, color="#d62728", lw=2)
        ax.text(obs, ax.get_ylim()[1] * 0.9, f"  observed = {obs}\n  z = {z:.1f}", color="#d62728", fontsize=9)
        ax.set_xlabel("TF-TF physical-evidence edge count")
        ax.set_ylabel("Permutations")
        ax.set_title(ttl, fontsize=10)
    fig.suptitle("Figure 4. TF-TF hyper-assortativity in the physical-evidence PPI network", fontsize=11, y=1.03)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figure4_assortativity.png", bbox_inches="tight")
    plt.close(fig)

    # ---------------- Figure 5: family-cluster bootstrap ----------------
    rs_nh = bootstrap(non_homeo, rng, 1500)
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    ax.hist(rs_boot, bins=40, alpha=0.6, label=f"Full cohort (mean={rs_boot.mean():.2f})", color="#1f77b4")
    ax.hist(rs_nh, bins=40, alpha=0.6, label=f"Non-homeodomain (mean={rs_nh.mean():.2f})", color="#ff7f0e")
    ax.axvline(0, color="black", lw=1, ls=":")
    ax.set_xlabel("Partial correlation r (per bootstrap resample)")
    ax.set_ylabel("Resamples")
    ax.set_title("Figure 5. TF-family cluster bootstrap\n(resampling whole DNA-binding-domain families)", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/figure5_cluster_bootstrap.png", bbox_inches="tight")
    plt.close(fig)

    print("All figures written to", OUTDIR)
    for f in sorted(os.listdir(OUTDIR)):
        if f.endswith(".png"):
            print(" -", f)

    fdf.to_csv(f"{RES}/figure2_forest_table.tsv", sep="\t", index=False)
    summary = {
        "n_full": len(sub), "n_non_homeo": len(non_homeo),
        "r_raw_full": r_raw, "p_raw_full": p_raw, "r_partial_full": r_p, "p_partial_full": p_p,
        "r_partial_nonhomeo": r_nh, "p_partial_nonhomeo": p_nh,
        "assort_obs": int(obs), "assort_naive_mean": float(naive_null.mean()), "assort_naive_sd": float(naive_null.std()),
        "assort_matched_mean": float(matched_null.mean()), "assort_matched_sd": float(matched_null.std()),
        "boot_full_mean": float(rs_boot.mean()), "boot_full_lo": float(lo), "boot_full_hi": float(hi),
        "boot_nonhomeo_mean": float(rs_nh.mean()),
    }
    json.dump(summary, open(f"{RES}/manuscript_numbers.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
