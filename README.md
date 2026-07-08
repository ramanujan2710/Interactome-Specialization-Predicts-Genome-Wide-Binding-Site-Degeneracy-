# Interactome Specialization Predicts Genome-Wide Binding-Site Degeneracy in *Drosophila melanogaster* Transcription Factors

This repository contains the complete, end-to-end analysis pipeline behind the manuscript *"Interactome Specialization Predicts Genome-Wide Binding-Site Degeneracy in Drosophila melanogaster Transcription Factors"* (working draft). It reproduces every reported number and figure directly from public STRING, JASPAR, phyloP, Ensembl, and FlyBase source files.

## What this shows

We test whether the fraction of a transcription factor's (TF's) physical protein-protein interaction (PPI) partners that are themselves other TFs — the "module fraction," computed from an evidence-channel-purified subset of STRING (experimental + curated-database evidence only) — predicts that TF's genome-wide predicted DNA-binding-site density. Across 279 confirmed *D. melanogaster* TFs, module fraction correlates with binding-site degeneracy (partial r ≈ 0.28–0.38 depending on control set), survives an extensive confound and robustness battery, and is accompanied by an independent finding that TFs are strongly hyper-assortative in the physical PPI network (z up to ~180 against a degree-matched permutation null).

The manuscript is not included in this repository; this repo is the code only.

## Repository contents

- `tf_ppi_full_pipeline.py` — the full analysis, as 14 sequential, checkpointed stages (see below).
- `make_figures.py` — regenerates the manuscript figures (scatter plots, confound forest plot, permutation-null histograms, bootstrap distributions) from pipeline output tables.
- `references.json` — the manuscript's reference list, if you want the citation context alongside the code.

## Pipeline stages

| Stage | What it does |
|---|---|
| 1–2 | Identify TFs from STRING annotation (DNA-binding-domain keyword regex) and parse JASPAR PFMs into per-position Shannon information content. |
| 3–4 | Match JASPAR motifs to STRING protein IDs; build the combined-score PPI network; classify DNA-binding-domain family. |
| 5–6 | Compute expression breadth / co-expression from FlyBase RPKM; build two independent ~18 Mb genomic samples (gene-body-restricted and intergenic-restricted) from the Ensembl GTF. |
| 7–8 | Fetch JASPAR2026 genome-wide predicted sites and phyloP conservation over both samples; compare genic vs. intergenic compartments. |
| 9–10 | Decompose STRING evidence channels into a physical-evidence-only network; run the core module-fraction-vs-site-density test. |
| 11 | Composite specialization index (module fraction + inverse expression breadth); TF-TF assortativity via permutation nulls. |
| 12–13 | Confound battery (motif length, GC/AT composition, family/paralog exclusion); TF-family cluster bootstrap. |
| 14 | Adversarial robustness checks: tautology/collinearity, Spearman vs. Pearson, outlier leverage, winsorizing, post-hoc power. |

Each stage writes its output to `OUTPUT_DIR` and reads only from earlier stages' outputs, so any stage can be re-run in isolation once its dependencies exist.

## A pitfall worth knowing about before you modify this

Stage 7 deliberately applies **no PWM match-score threshold** when counting predicted binding sites (`SCORE_THRESHOLD = 0` in the CONFIG block). An earlier version of this pipeline used a fixed threshold of 400, which silently destroyed the core result: the fraction of a motif's sites discarded by any fixed absolute score cutoff correlates strongly with motif length (r ≈ −0.73), because raw PWM scores are not comparable across motifs of different length or information content. If you reintroduce score filtering, do it with a length-normalized or background-calibrated significance model (e.g. FIMO-style p-values), not a raw score cutoff — see the inline comment above `SCORE_THRESHOLD` in the script.

## Requirements

```bash
pip install numpy pandas scipy pyBigWig matplotlib
```

Python 3.9+ recommended. `pyBigWig` requires a working C build toolchain on some platforms (see [pyBigWig's install notes](https://github.com/deeptools/pyBigWig)).

## Input data

Download these into a local `data/` directory (not included in this repo):

| File | Source |
|---|---|
| `7227.protein.info.v12.0.txt` | [STRING](https://string-db.org/cgi/download) — *D. melanogaster* (taxon 7227) protein info, v12.0 |
| `7227.protein.links.v12.0.txt` | STRING — combined-score network, v12.0 |
| `7227.protein.links.full.v12.0.txt` | STRING — full evidence-channel-split network, v12.0 |
| `Jaspar.txt` | [JASPAR](https://jaspar.elixir.no/downloads/) — CORE *D. melanogaster* PFMs (flat-file / count-matrix format) |
| `JASPAR2026_dm6.bb` | JASPAR — genome-wide predicted TFBS track for dm6 (bigBed) |
| `dm6.phyloP124way.bw` | [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/) — 124-way phyloP conservation (bigWig) |
| `Drosophila_melanogaster.BDGP6.54.63.gtf` | [Ensembl](https://ftp.ensembl.org/pub/) — BDGP6 gene annotation |
| `gene_rpkm_matrix_fb_2026_02.tsv` | [FlyBase](https://flybase.org/downloads/bulkdata) — high-throughput gene expression RPKM matrix |

## Running it

```bash
python tf_ppi_full_pipeline.py                 # run every stage in order
python tf_ppi_full_pipeline.py --stage 8        # run/resume a single stage
```

`INPUT_DIR` and `OUTPUT_DIR` are set at the top of `tf_ppi_full_pipeline.py` (default `./data` and `./results`); edit them or symlink your data directory into place.

Stages 6–7 (genome-wide bigBed/phyloP fetching over ~18 Mb of sequence, twice) are the slowest part of the run. They checkpoint to disk in `TIME_BUDGET`-second slices and print `RESUME NEEDED` vs. `DONE`, so if you're running under a timeout-limited shell (CI, a notebook cell, a sandbox), just re-invoke the same `--stage` call until it reports `DONE`. On an unconstrained machine the whole pipeline runs in one pass without any special handling.

Once the pipeline has produced `results/physical_module_full.tsv`, `results/compartment_compare_full.tsv`, `results/physical_edges.npz`, `results/coexpr_matrix.npz`, and `results/protein_tf_flag.tsv`, run:

```bash
python make_figures.py
```

to regenerate the manuscript figures into `figures/`.

## Reproducing specific numbers

- Core result (Stage 10 printout): physical-evidence module fraction vs. genome-wide site density, raw and partial correlations, full cohort and non-homeodomain stratum.
- Confound battery (Stage 12): motif-length control, GC/AT composition control, family-exclusion subsets, Benjamini-Hochberg-corrected p-values.
- Robustness battery (Stage 14): tautology check against motif information content, Spearman vs. Pearson, outlier-trimming and winsorizing sensitivity, post-hoc power.
- TF-TF assortativity (Stage 11): naive and degree-matched permutation nulls, with bin-count sensitivity.

## License

Add a license of your choosing (e.g. MIT) before publishing this repository publicly. No license is specified by default.

## Citation

If you use this pipeline, please cite the accompanying manuscript (details to be added upon submission/preprint posting).
