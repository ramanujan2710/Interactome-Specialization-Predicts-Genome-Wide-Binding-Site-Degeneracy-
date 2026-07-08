#!/usr/bin/env python3
"""
CLI entrypoint for the tf_ppi analysis pipeline.

    python run_pipeline.py                # run every stage in order
    python run_pipeline.py --stage 8       # run (or resume) a single stage
    python run_pipeline.py --list          # list all stages and exit

Stages checkpoint their outputs to OUTPUT_DIR (default ./results), so
re-running an already-completed stage is cheap and safe, and Stage 7 (the
slowest step, fetching genome-wide predicted sites) is internally
resumable -- simply re-invoke `--stage 7` until it reports DONE for both
genomic samples.
"""
import argparse
import sys

sys.path.insert(0, "src")

from tf_ppi import (
    stage01_tf_identification,
    stage02_jaspar_information_content,
    stage03_network_construction,
    stage04_family_classification,
    stage05_expression_breadth,
    stage06_genomic_sampling,
    stage07_fetch_genomic_data,
    stage08_compartment_comparison,
    stage09_evidence_decomposition,
    stage10_core_test,
    stage11_composite_assortativity,
    stage12_confound_battery,
    stage13_cluster_bootstrap,
    stage14_adversarial_checks,
)

STAGES = [
    ("1",  "TF identification",                              stage01_tf_identification.run),
    ("2",  "JASPAR information content",                     stage02_jaspar_information_content.run),
    ("3",  "JASPAR<->STRING merge + combined-score network",  stage03_network_construction.run),
    ("4",  "DNA-binding-domain family classification",       stage04_family_classification.run),
    ("5",  "Expression breadth / co-expression",              stage05_expression_breadth.run),
    ("6",  "Build genomic samples (genic / intergenic)",      stage06_genomic_sampling.run),
    ("7",  "Fetch bigBed/phyloP over genomic samples",        stage07_fetch_genomic_data.run),
    ("8",  "Compartment comparison (genic vs. intergenic)",   stage08_compartment_comparison.run),
    ("9",  "STRING evidence-channel decomposition",           stage09_evidence_decomposition.run),
    ("10", "Physical module fraction (CORE TEST)",            stage10_core_test.run),
    ("11", "Composite index + TF-TF assortativity",           stage11_composite_assortativity.run),
    ("12", "Confound battery",                                stage12_confound_battery.run),
    ("13", "Family-cluster bootstrap",                        stage13_cluster_bootstrap.run),
    ("14", "Adversarial artifact / outlier / tautology checks", stage14_adversarial_checks.run),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", default=None, help="Run only this stage number (e.g. 8); omit to run all stages in order.")
    ap.add_argument("--list", action="store_true", help="List all stages and exit.")
    args = ap.parse_args()

    if args.list:
        for num, desc, _ in STAGES:
            print(f"  {num:>2}  {desc}")
        return

    for num, desc, fn in STAGES:
        if args.stage is not None and args.stage != num:
            continue
        print(f"\n{'='*80}\nSTAGE {num}: {desc}\n{'='*80}")
        fn()


if __name__ == "__main__":
    main()
