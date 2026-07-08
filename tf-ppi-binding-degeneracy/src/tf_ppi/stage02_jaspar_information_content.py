"""
Stage 2 -- Parse JASPAR PFMs, compute per-position and average Shannon
information content (bits) with a standard 0.25 pseudocount, plus consensus
GC fraction (used later as a composition control distinct from information
content).

Output: jaspar_ic.json
"""
import json
import math
import re

import numpy as np

from .config import PATHS, out_path


def _parse_jaspar(path):
    motifs = []
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if line.startswith(">"):
            parts = line[1:].split("\t")
            ma_id = parts[0]
            name = parts[1] if len(parts) > 1 else ""
            rows = {}
            i += 1
            for _ in "ACGT":
                row = lines[i].strip()
                m = re.match(r"^([ACGT])\s*\[(.*)\]$", row)
                base = m.group(1)
                vals = [float(x) for x in m.group(2).split()]
                rows[base] = vals
                i += 1
            motifs.append((ma_id, name, rows))
        else:
            i += 1
    return motifs


def _info_content(rows):
    """Per-position Shannon information content (bits), standard PWM->PPM->IC
    conversion with a 0.25-per-base pseudocount."""
    L = len(rows["A"])
    bpp = []
    for p in range(L):
        col = [rows[b][p] for b in "ACGT"]
        col_pc = [c + 0.25 for c in col]
        s = sum(col_pc)
        probs = [c / s for c in col_pc]
        H = -sum(pi * math.log2(pi) for pi in probs if pi > 0)
        bpp.append(2 - H)
    return bpp


def run():
    print("\n### STAGE 2: JASPAR PWM parsing + information content ###")
    motifs = _parse_jaspar(PATHS["jaspar"])
    print("total motifs parsed:", len(motifs))

    results = []
    for ma_id, name, rows in motifs:
        L = len(rows["A"])
        bpp = _info_content(rows)
        total_ic = sum(bpp)
        avg_ic = total_ic / L if L else 0
        colsum = np.array([sum(rows[b]) for b in "ACGT"])
        gc_frac = float((colsum[1] + colsum[2]) / colsum.sum()) if colsum.sum() > 0 else np.nan
        results.append({"ma_id": ma_id, "name": name, "length": L,
                         "total_ic": total_ic, "avg_ic": avg_ic,
                         "consensus_gc_frac": gc_frac})

    with open(out_path("jaspar_ic.json"), "w") as out:
        json.dump(results, out)
    print(f"avg_ic range: {min(r['avg_ic'] for r in results):.2f} - "
          f"{max(r['avg_ic'] for r in results):.2f} bits/position")


if __name__ == "__main__":
    run()
