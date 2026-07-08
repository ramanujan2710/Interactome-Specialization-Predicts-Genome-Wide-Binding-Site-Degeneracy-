"""
Stage 6 -- Build two independent ~18 Mb genomic sampling regions from the
GTF: (a) gene-body-restricted, (b) intergenic-restricted. These are two
statistically independent replication samples used throughout the rest of
the pipeline.

Outputs: gene_sample_regions.json, intergenic_sample_regions.json
"""
import json
from collections import defaultdict

import numpy as np

from .config import PATHS, out_path, CHROM_MAP, SAMPLE_CAP_BP, SAMPLE_TARGET_BP, INTERGENIC_BUFFER_BP


def _build_gene_body_sample():
    genes = []
    with open(PATHS["gtf"]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fs = line.rstrip("\n").split("\t")
            if fs[2] != "gene" or 'gene_biotype "protein_coding"' not in fs[8]:
                continue
            chrom, start, end = fs[0], int(fs[3]), int(fs[4])
            if chrom not in CHROM_MAP:
                continue
            genes.append((CHROM_MAP[chrom], start, end))
    print("protein-coding genes on main chroms:", len(genes))

    rng = np.random.default_rng(2024)
    rng.shuffle(genes)
    selected, total = [], 0
    for (c, s, e) in genes:
        e2 = min(e, s + SAMPLE_CAP_BP)
        span = e2 - s
        if span <= 0:
            continue
        selected.append((c, s, e2))
        total += span
        if total >= SAMPLE_TARGET_BP:
            break
    print("gene-body sample:", len(selected), "regions,", total, "bp")
    json.dump(selected, open(out_path("gene_sample_regions.json"), "w"))


def _build_intergenic_sample():
    import pyBigWig

    genes_by_chrom = defaultdict(list)
    with open(PATHS["gtf"]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fs = line.rstrip("\n").split("\t")
            if fs[2] != "gene":
                continue
            chrom, start, end = fs[0], int(fs[3]), int(fs[4])
            if chrom not in CHROM_MAP:
                continue
            genes_by_chrom[CHROM_MAP[chrom]].append((start, end))

    bw = pyBigWig.open(PATHS["phylop"])
    chrom_lens = bw.chroms()
    bw.close()

    intergenic = []
    for c, ivs in genes_by_chrom.items():
        ivs.sort()
        merged = []
        cs, ce = ivs[0][0] - INTERGENIC_BUFFER_BP, ivs[0][1] + INTERGENIC_BUFFER_BP
        for (s, e) in ivs[1:]:
            s2, e2 = s - INTERGENIC_BUFFER_BP, e + INTERGENIC_BUFFER_BP
            if s2 <= ce:
                ce = max(ce, e2)
            else:
                merged.append((cs, ce))
                cs, ce = s2, e2
        merged.append((cs, ce))
        L = chrom_lens.get(c, 10**9)
        prev_end = 0
        for (s, e) in merged:
            s = max(s, 0)
            if s > prev_end:
                intergenic.append((c, prev_end, min(s, L)))
            prev_end = max(prev_end, e)
        if prev_end < L:
            intergenic.append((c, prev_end, L))
    intergenic = [(c, s, e) for (c, s, e) in intergenic if e - s > 500]
    print("intergenic intervals:", len(intergenic), "total bp:", sum(e - s for _, s, e in intergenic))

    rng = np.random.default_rng(555)
    rng.shuffle(intergenic)
    selected, total = [], 0
    for (c, s, e) in intergenic:
        span = e - s
        if span > SAMPLE_CAP_BP:
            offset = rng.integers(0, span - SAMPLE_CAP_BP)
            s2, e2 = s + offset, s + offset + SAMPLE_CAP_BP
        else:
            s2, e2 = s, e
        selected.append((c, int(s2), int(e2)))
        total += (e2 - s2)
        if total >= SAMPLE_TARGET_BP:
            break
    print("intergenic sample:", len(selected), "regions,", total, "bp")
    json.dump(selected, open(out_path("intergenic_sample_regions.json"), "w"))


def run():
    print("\n### STAGE 6: build gene-body and intergenic genomic samples ###")
    _build_gene_body_sample()
    _build_intergenic_sample()


if __name__ == "__main__":
    run()
