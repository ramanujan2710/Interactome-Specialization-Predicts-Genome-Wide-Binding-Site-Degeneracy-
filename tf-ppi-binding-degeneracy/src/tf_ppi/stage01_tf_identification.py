"""
Stage 1 -- Identify TFs from STRING annotation via DNA-binding-domain
keyword regex (structural keywords, not the bare phrase "transcription
factor" alone, to reduce annotation-text circularity).

Output: protein_tf_flag.tsv
"""
import re
from .config import PATHS, out_path

TF_KEYWORDS = [
    r"zinc finger", r"zinc-finger", r"homeobox", r"homeodomain", r"bHLH",
    r"basic helix-loop-helix", r"bZIP", r"basic leucine zipper", r"forkhead",
    r"POU domain", r"HMG box", r"HMG-box", r"nuclear receptor", r"T-box",
    r"ETS domain", r"ETS-domain", r"winged helix", r"paired box", r"Pax domain",
    r"bHLH-PAS", r"GATA-type", r"C2H2", r"leucine zipper", r"AP-2 domain",
    r"Runt domain", r"Sox domain", r"MADS-box", r"MADS box", r"basic-leucine zipper",
    r"transcription factor", r"DNA-binding domain", r"sequence-specific DNA binding",
]
TF_RE = re.compile("|".join(TF_KEYWORDS), re.IGNORECASE)


def run():
    print("\n### STAGE 1: TF identification from STRING annotation ###")
    rows = []
    with open(PATHS["info"], encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            pid, name, size, annot = parts[0], parts[1], parts[2], parts[3]
            is_tf = bool(TF_RE.search(annot)) or bool(TF_RE.search(name))
            rows.append((pid, name, size, len(annot), is_tf))

    print("total proteins:", len(rows))
    tf_count = sum(r[4] for r in rows)
    print("TF-like annotated proteins:", tf_count, f"({tf_count/len(rows):.1%})")

    with open(out_path("protein_tf_flag.tsv"), "w") as out:
        out.write("protein_id\tname\tsize\tannot_len\tis_tf\n")
        for r in rows:
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{int(r[4])}\n")
    print("wrote protein_tf_flag.tsv")


if __name__ == "__main__":
    run()
