# tf-ppi-binding-degeneracy

Analysis code for **"Interactome Specialization Predicts Genome-Wide
Binding-Site Degeneracy in *Drosophila melanogaster* Transcription
Factors."**

This repository reproduces, end to end and from raw source files, every
statistical result and figure reported in the associated manuscript: a
confound-controlled, multiply-replicated analysis testing whether the
fraction of a transcription factor's physical protein-protein interaction
(PPI) partners that are themselves other TFs ("module fraction") predicts
that TF's genome-wide predicted DNA-binding-site density.

## Headline result

Using an evidence-channel-purified STRING network (experimental and
curated-database evidence only) and 279 confirmed *D. melanogaster* TFs,
module fraction correlates with genome-wide binding-site density
(partial r = 0.38, p = 5.6×10⁻¹¹), survives a full battery of confound
controls (motif length, GC/AT composition, TF-family clustering), an
adversarial robustness battery (rank correlation, outlier trimming,
winsorizing, post-hoc power), and independent replication across two
genomic compartments. See the manuscript for full results, discussion, and
an explicitly documented limitations section.

**This repository also documents a methodological pitfall discovered
during development:** applying a fixed PWM match-score threshold when
counting predicted binding sites reintroduces a severe motif-length-
dependent artifact that suppresses the core result. See `config.py` and
the manuscript Results section for the full diagnosis. The pipeline is
configured correctly (`SCORE_THRESHOLD = 0`) by default — this note exists
so nobody re-introduces the bug.

## Repository structure

```
tf-ppi-binding-degeneracy/
├── run_pipeline.py              # CLI entrypoint: runs all 14 stages in order
├── src/tf_ppi/
│   ├── config.py                 # paths, constants (incl. the score-threshold note)
│   ├── utils.py                  # shared statistical helpers
│   ├── stage01_tf_identification.py
│   ├── stage02_jaspar_information_content.py
│   ├── stage03_network_construction.py
│   ├── stage04_family_classification.py
│   ├── stage05_expression_breadth.py
│   ├── stage06_genomic_sampling.py
│   ├── stage07_fetch_genomic_data.py
│   ├── stage08_compartment_comparison.py
│   ├── stage09_evidence_decomposition.py
│   ├── stage10_core_test.py               # <- headline result
│   ├── stage11_composite_assortativity.py
│   ├── stage12_confound_battery.py
│   ├── stage13_cluster_bootstrap.py
│   └── stage14_adversarial_checks.py      # <- regression test for the bug above
├── figures/
│   └── make_figures.py           # regenerates all 5 manuscript figures
├── tests/
│   └── test_utils.py             # fast synthetic-data sanity tests (pytest)
├── data/
│   └── README.md                 # what to download, and from where
├── results/                      # pipeline outputs land here (gitignored)
├── requirements.txt
├── pyproject.toml
└── LICENSE
```

## Installation

```bash
git clone <this-repo>
cd tf-ppi-binding-degeneracy
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Data

This pipeline requires eight raw input files (STRING, JASPAR, phyloP, GTF,
FlyBase — none of them ours to redistribute). See
[`data/README.md`](data/README.md) for exact filenames and download links.
Place them in `./data/`, or point `TF_PPI_INPUT_DIR` at wherever you keep
them.

## Running the pipeline

```bash
# run every stage in order
python run_pipeline.py

# list all stages
python run_pipeline.py --list

# run (or resume) a single stage
python run_pipeline.py --stage 10
```

Every stage checkpoints its output to `./results/`, so re-running an
already-completed stage is cheap. **Stage 7** (fetching genome-wide
predicted binding sites and conservation scores over two independent
18 Mb genomic samples) is the slowest step and is internally resumable in
`TF_PPI_TIME_BUDGET_SEC`-second slices — if you're running under a
timeout-limited shell, just re-invoke `--stage 7` until it prints `DONE`
for both samples.

After all 14 stages complete:

```bash
cd figures
python make_figures.py
```

This regenerates `figure1_core_scatter.png` through
`figure5_cluster_bootstrap.png` directly from the pipeline's output
tables, plus a `manuscript_numbers.json` summary in `results/`.

## Tests

```bash
pip install -r requirements.txt pytest
pytest tests/
```

These are fast synthetic-data sanity checks on the shared statistical
helpers (rank-z transform, OLS residualization, STRING-style probabilistic
evidence combination, Fisher z-test, post-hoc power, winsorizing). They do
not require the raw data files. The real regression test for the analysis
itself is **Stage 14** (`run_pipeline.py --stage 14`), which runs the
adversarial artifact/outlier/tautology battery against the actual pipeline
output and is what originally caught the score-threshold bug described
above.

## Pipeline stages at a glance

| # | Stage | Key output |
|---|---|---|
| 1 | TF identification (STRING annotation regex) | `protein_tf_flag.tsv` |
| 2 | JASPAR PWM parsing + information content | `jaspar_ic.json` |
| 3 | JASPAR↔STRING merge, combined-score network | `merged_full.tsv`, `arrays.npz` |
| 4 | DNA-binding-domain family classification | `tf_with_family_full.tsv` |
| 5 | Expression breadth / co-expression (FlyBase RPKM) | `coexpr_matrix.npz` |
| 6 | Build 18 Mb genic + intergenic genomic samples | `*_sample_regions.json` |
| 7 | Fetch predicted sites + phyloP (checkpointed) | `*_sample_agg_state.json` |
| 8 | Genic vs. intergenic compartment comparison | `compartment_compare_full.tsv` |
| 9 | STRING evidence-channel decomposition | `physical_edges.npz` |
| 10 | **Core test**: module fraction vs. site density | `physical_module_full.tsv` |
| 11 | Composite index + TF-TF assortativity permutation tests | stdout |
| 12 | Confound battery (motif length, GC content, family exclusion) | stdout |
| 13 | TF-family cluster bootstrap | stdout |
| 14 | Adversarial robustness battery | stdout |

## Citation

If you use this code, please cite the associated manuscript (see
`CITATION.cff` if provided, or cite this repository directly by URL and
commit hash pending formal publication).

## License

MIT — see [`LICENSE`](LICENSE).
