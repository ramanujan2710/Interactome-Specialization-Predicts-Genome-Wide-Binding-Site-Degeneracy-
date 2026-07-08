"""
Stage 9 -- STRING evidence-channel decomposition: build a PHYSICAL-EVIDENCE-
ONLY network (experiments + curated-database channels only), explicitly
excluding textmining, co-expression, genomic-context, and homology channels.
This is the evidence-purified network underlying the headline result.

Output: physical_edges.npz
"""
import numpy as np
import pandas as pd

from .config import PATHS, out_path
from .utils import combine_probs

EVIDENCE_COLS = ["p1", "p2", "neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
                  "homology", "coexpression", "coexpression_transferred", "experiments",
                  "experiments_transferred", "database", "database_transferred", "textmining",
                  "textmining_transferred", "combined_score"]


def run():
    print("\n### STAGE 9: STRING evidence-channel decomposition (physical-only network) ###")
    df = pd.read_csv(PATHS["links_full"], sep=" ", header=0, names=EVIDENCE_COLS,
                      dtype={c: (str if c in ["p1", "p2"] else np.int32) for c in EVIDENCE_COLS})
    print("rows:", len(df))

    df["physical_score"] = combine_probs(df["experiments"], df["experiments_transferred"],
                                          df["database"], df["database_transferred"]) * 1000
    df["textmining_only_score"] = combine_probs(df["textmining"], df["textmining_transferred"]) * 1000
    print("physical_score>0 rows:", (df["physical_score"] > 0).sum())

    flag = pd.read_csv(out_path("protein_tf_flag.tsv"), sep="\t")
    idx_map = {pid: i for i, pid in enumerate(flag["protein_id"])}
    df["i1"] = df["p1"].map(idx_map)
    df["i2"] = df["p2"].map(idx_map)
    df = df.dropna(subset=["i1", "i2"])
    df["i1"] = df["i1"].astype(np.int32)
    df["i2"] = df["i2"].astype(np.int32)
    df = df[df["i1"] < df["i2"]]
    print("deduped rows:", len(df))

    phys = df[df["physical_score"] > 0][["i1", "i2", "physical_score", "combined_score", "textmining_only_score"]]
    print("physical edges kept:", len(phys))
    np.savez(out_path("physical_edges.npz"), i1=phys["i1"].values, i2=phys["i2"].values,
             physical_score=phys["physical_score"].values, combined_score=phys["combined_score"].values,
             textmining_only_score=phys["textmining_only_score"].values)


if __name__ == "__main__":
    run()
