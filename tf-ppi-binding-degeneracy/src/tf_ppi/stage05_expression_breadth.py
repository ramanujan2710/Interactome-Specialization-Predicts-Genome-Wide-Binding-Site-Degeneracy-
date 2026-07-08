"""
Stage 5 -- Expression breadth and pairwise co-expression from the FlyBase
RPKM matrix.

Output: coexpr_matrix.npz (coexpr_corr, names, breadth)
"""
import numpy as np
import pandas as pd

from .config import PATHS, out_path


def run():
    print("\n### STAGE 5: FlyBase RPKM -> expression breadth / co-expression ###")
    path = PATHS["rpkm"]
    with open(path) as f:
        lines = f.readlines()
    header_line = [l for l in lines if l.startswith("#gene_primary_id")][0]
    cols = header_line.strip("#\n").split("\t")
    df = pd.read_csv(path, sep="\t", skiprows=5, header=None)
    df.columns = cols
    print("RPKM matrix shape:", df.shape)

    net = pd.read_csv(out_path("tf_with_family_full.tsv"), sep="\t")
    tf_symbols = set(net["name_lc"])
    df["sym_lc"] = df["gene_symbol"].str.lower()
    matched = df[df["sym_lc"].isin(tf_symbols)].copy()
    print("matched TFs in RPKM matrix:", len(matched), "/", len(tf_symbols))

    meta_cols = ["gene_primary_id", "gene_symbol", "gene_fullname", "gene_type", "sym_lc"]
    expr_cols = [c for c in matched.columns if c not in meta_cols]
    expr_mat = np.log2(matched[expr_cols].values.astype(float) + 1)
    names = matched["sym_lc"].values

    X = expr_mat - expr_mat.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1e-9
    Xn = X / norms
    coexpr_corr = Xn @ Xn.T
    breadth = (matched[expr_cols].values.astype(float) > 1).mean(axis=1)

    np.savez(out_path("coexpr_matrix.npz"), coexpr_corr=coexpr_corr, names=names, breadth=breadth)
    print("co-expression matrix:", coexpr_corr.shape, " breadth computed for", len(names), "TFs")


if __name__ == "__main__":
    run()
