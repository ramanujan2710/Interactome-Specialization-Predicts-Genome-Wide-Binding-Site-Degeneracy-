"""
Stage 4 -- DNA-binding-domain family classification, independent of Stage 1's
keyword set. Used for the homeodomain-exclusion robustness check and the
TF-family cluster bootstrap (Stage 13).

Output: tf_with_family_full.tsv
"""
import re

import pandas as pd

from .config import PATHS, out_path

FAMILY_PATTERNS = {
    "homeodomain": [r"homeobox", r"homeodomain"],
    "zinc_finger": [r"zinc finger", r"zinc-finger", r"C2H2", r"GATA-type", r"AP-2 domain"],
    "bHLH": [r"basic helix-loop-helix", r"\bbHLH\b", r"bHLH-PAS", r"bHLH-zip"],
    "bZIP": [r"basic leucine zipper", r"\bbZIP\b", r"leucine zipper"],
    "nuclear_receptor": [r"nuclear receptor", r"nuclear hormone receptor"],
    "forkhead": [r"forkhead", r"winged helix", r"winged-helix"],
    "ETS": [r"ETS domain", r"ETS-domain"],
    "HMG_box": [r"HMG box", r"HMG-box", r"Sox domain"],
    "MADS_box": [r"MADS-box", r"MADS box"],
    "T_box": [r"T-box"],
    "Runt": [r"Runt domain"],
    "POU": [r"POU domain"],
    "Pax_paired": [r"paired box", r"Pax domain"],
}


def _classify_family(annot, name):
    text = str(annot) + " " + str(name)
    hits = []
    for fam, pats in FAMILY_PATTERNS.items():
        for pat in pats:
            if re.search(pat, text, re.IGNORECASE):
                hits.append(fam)
                break
    return hits[0] if len(hits) == 1 else ("multiple" if len(hits) > 1 else "other/unclassified")


def run():
    print("\n### STAGE 4: DNA-binding-domain family classification ###")
    info = pd.read_csv(PATHS["info"], sep="\t")
    info["name_lc"] = info["preferred_name"].str.lower()
    info["family"] = info.apply(lambda r: _classify_family(r["annotation"], r["preferred_name"]), axis=1)

    m = pd.read_csv(out_path("merged_full.tsv"), sep="\t")
    m = m.merge(info[["name_lc", "family"]], on="name_lc", how="left")
    print(m["family"].value_counts().to_string())
    m.to_csv(out_path("tf_with_family_full.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    run()
