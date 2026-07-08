# Input data

Raw third-party data files are not committed to this repository (they are
large, and most are redistributable-but-not-ours-to-host). Download the
following into this `data/` directory before running the pipeline:

| File | Source | Notes |
|---|---|---|
| `7227.protein.info.v12.0.txt` | [STRING](https://string-db.org/cgi/download) | *D. melanogaster* (taxon 7227), protein annotation table |
| `7227.protein.links.v12.0.txt` | [STRING](https://string-db.org/cgi/download) | combined-score network |
| `7227.protein.links.full.v12.0.txt` | [STRING](https://string-db.org/cgi/download) | full evidence-channel-split network |
| `Jaspar.txt` | [JASPAR](https://jaspar.elixir.no/downloads/) | CORE *D. melanogaster* PFM/count-matrix flat file |
| `JASPAR2026_dm6.bb` | [JASPAR UCSC track hub](https://jaspar.elixir.no/genome-tracks/) | genome-wide predicted TFBS, dm6, bigBed |
| `dm6.phyloP124way.bw` | [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/) | 124-way phyloP conservation, bigWig |
| `Drosophila_melanogaster.BDGP6.54.63.gtf` | [Ensembl](https://ftp.ensembl.org/pub/release-63/gtf/drosophila_melanogaster/) (or [Ensembl Metazoa](https://metazoa.ensembl.org/Drosophila_melanogaster/Info/Index)) | gene annotation, BDGP6 assembly |
| `gene_rpkm_matrix_fb_2026_02.tsv` | [FlyBase](https://flybase.org/downloads) | high-throughput gene expression RPKM matrix |

Set `TF_PPI_INPUT_DIR` (or edit `src/tf_ppi/config.py`) if you keep these
files somewhere other than `./data`.
