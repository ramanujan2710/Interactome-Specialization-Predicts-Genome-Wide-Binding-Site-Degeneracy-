"""
Central configuration: input/output paths and pipeline constants.

Edit INPUT_DIR / OUTPUT_DIR (or set the TF_PPI_INPUT_DIR / TF_PPI_OUTPUT_DIR
environment variables) to point at your local copies of the required data
files. See data/README.md for exactly what to download and from where.
"""
import os

INPUT_DIR = os.environ.get("TF_PPI_INPUT_DIR", "./data")
OUTPUT_DIR = os.environ.get("TF_PPI_OUTPUT_DIR", "./results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PATHS = {
    "info":       os.path.join(INPUT_DIR, "7227.protein.info.v12.0.txt"),
    "links":      os.path.join(INPUT_DIR, "7227.protein.links.v12.0.txt"),
    "links_full": os.path.join(INPUT_DIR, "7227.protein.links.full.v12.0.txt"),
    "jaspar":     os.path.join(INPUT_DIR, "Jaspar.txt"),
    "bigbed":     os.path.join(INPUT_DIR, "JASPAR2026_dm6.bb"),
    "phylop":     os.path.join(INPUT_DIR, "dm6.phyloP124way.bw"),
    "gtf":        os.path.join(INPUT_DIR, "Drosophila_melanogaster.BDGP6.54.63.gtf"),
    "rpkm":       os.path.join(INPUT_DIR, "gene_rpkm_matrix_fb_2026_02.tsv"),
}


def out_path(name: str) -> str:
    """Resolve a filename to a path inside OUTPUT_DIR."""
    return os.path.join(OUTPUT_DIR, name)


# ----------------------------------------------------------------------------
# Score threshold applied to JASPAR2026 bigBed predictions when counting a
# genomic position as a "predicted binding site" for a TF.
#
# THIS MUST BE 0 (no filtering -- every predicted entry counts).
#
# An earlier iteration of this pipeline applied a fixed threshold of 400
# here, which reintroduces -- in a new guise -- the exact cross-motif score
# comparability confound identified and retracted elsewhere in this project:
# raw JASPAR PWM match scores are not comparable across motifs of different
# length or information content. Empirically, on an identical 18 Mb
# gene-body genomic sample, threshold=400 collapses the core correlation
# (partial r drops from ~0.35 to ~0.13, and loses significance under
# further controls), because the fraction of a motif's sites discarded by
# ANY fixed absolute score cutoff correlates strongly with motif length
# (r=-0.73). threshold=0 reproduces the original, independently-replicated
# result. Do not reintroduce a nonzero threshold without first deriving a
# length-normalized, background-calibrated significance model (e.g. a
# FIMO-style p-value) for cross-motif comparability.
# ----------------------------------------------------------------------------
SCORE_THRESHOLD = 0

# Genomic sample size (bp) for each of the two independent site-density samples.
SAMPLE_TARGET_BP = 18_000_000

# Per-gene / per-interval cap so no single locus (e.g. Dscam) dominates a sample.
SAMPLE_CAP_BP = 40_000

# Buffer excluded on each side of a gene when defining "intergenic" sequence.
INTERGENIC_BUFFER_BP = 2_000

# Time slice (seconds) used per invocation of the checkpointed bigBed/phyloP
# fetch stage; lower this if your shell/runner has a hard per-call timeout.
# The fetch stage is resumable and idempotent, so it is safe to re-invoke it
# repeatedly until it reports "DONE".
TIME_BUDGET_SEC = int(os.environ.get("TF_PPI_TIME_BUDGET_SEC", "36"))

CHROM_MAP = {"2L": "chr2L", "2R": "chr2R", "3L": "chr3L", "3R": "chr3R",
             "4": "chr4", "X": "chrX", "Y": "chrY"}
