"""
Stage 7 -- Checkpointed fetch of JASPAR2026 predicted sites and phyloP
conservation over each genomic sample built in Stage 6.

This stage is resumable and idempotent by design: it processes regions in
TIME_BUDGET_SEC-second slices, writes progress to a JSON checkpoint file
after each slice, and can simply be re-invoked (`run()` or `python -m
tf_ppi.stage07_fetch_genomic_data`) until it reports "DONE" for both
samples. This is useful when running under a sandboxed or timeout-limited
shell that cannot complete a multi-minute fetch in a single call.

IMPORTANT: no PWM match-score threshold is applied here (see config.py,
SCORE_THRESHOLD) -- every JASPAR2026 bigBed entry is counted. Applying a
fixed absolute threshold reintroduces a motif-length-dependent artifact
that suppresses the core result (see the project README / manuscript
Results for the empirical diagnosis).

Outputs: gene_sample_agg_state.json, intergenic_sample_agg_state.json
"""
import json
import os
import time
from collections import defaultdict

import numpy as np

from .config import PATHS, out_path, SCORE_THRESHOLD, TIME_BUDGET_SEC


def _fetch_and_aggregate(regions_file, state_file, tf_names):
    import pyBigWig

    t0 = time.time()
    with open(regions_file) as f:
        regions = json.load(f)

    if os.path.exists(state_file):
        state = json.load(open(state_file))
        next_idx = state["next_idx"]
        site_count = defaultdict(int, state["site_count"])
        phylop_sum = defaultdict(float, state["phylop_sum"])
        phylop_n = defaultdict(int, state["phylop_n"])
        total_bp_done = state["total_bp_done"]
    else:
        next_idx = 0
        site_count, phylop_sum, phylop_n = defaultdict(int), defaultdict(float), defaultdict(int)
        total_bp_done = 0

    bb = pyBigWig.open(PATHS["bigbed"])
    bw = pyBigWig.open(PATHS["phylop"])
    chroms = bb.chroms()

    n = len(regions)
    idx = next_idx
    while idx < n:
        if time.time() - t0 > TIME_BUDGET_SEC:
            break
        c, s, e = regions[idx]
        L = chroms.get(c)
        if L is None:
            idx += 1
            continue
        s2, e2 = max(s, 0), min(e, L)
        if e2 <= s2:
            idx += 1
            continue
        entries = bb.entries(c, s2, e2)
        if entries:
            vals = bw.values(c, s2, e2, numpy=True)
            for (est, eend, rest) in entries:
                fields = rest.split("\t")
                sc = float(fields[1])
                if sc < SCORE_THRESHOLD:
                    continue
                name = fields[3].lower() if len(fields) > 3 else "?"
                if name not in tf_names:
                    continue
                site_count[name] += 1
                cs, ce = est - s2, eend - s2
                seg = vals[max(cs, 0):min(ce, e2 - s2)]
                seg = seg[~np.isnan(seg)]
                if len(seg):
                    phylop_sum[name] += float(seg.sum())
                    phylop_n[name] += len(seg)
        total_bp_done += (e2 - s2)
        idx += 1
    bb.close()
    bw.close()

    state = {"next_idx": idx, "site_count": dict(site_count), "phylop_sum": dict(phylop_sum),
             "phylop_n": dict(phylop_n), "total_bp_done": total_bp_done}
    json.dump(state, open(state_file, "w"))
    done = idx >= n
    print(f"  processed to {idx}/{n} regions, total_bp={total_bp_done}, "
          f"TFs with hits={len(site_count)}, t={time.time()-t0:.1f}s, {'DONE' if done else 'RESUME NEEDED'}")
    return done


def run():
    print("\n### STAGE 7: fetch predicted sites + phyloP over both genomic samples ###")
    import pandas as pd
    tf_names = set(pd.read_csv(out_path("tf_with_family_full.tsv"), sep="\t")["name_lc"])

    print("-- gene-body sample --")
    while not _fetch_and_aggregate(out_path("gene_sample_regions.json"),
                                    out_path("gene_sample_agg_state.json"), tf_names):
        pass
    print("-- intergenic sample --")
    while not _fetch_and_aggregate(out_path("intergenic_sample_regions.json"),
                                    out_path("intergenic_sample_agg_state.json"), tf_names):
        pass


if __name__ == "__main__":
    run()
