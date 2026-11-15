# Seed microbiome host breadth

Analysis code for **Host-associated community structure and microbial host breadth in the seed microbiome**.

## Data and version

Code version: v1.0.0 (software archive pending). Processed inputs, precomputed
results and reference outputs are publicly archived as dataset v1.0.0:
[doi:10.5281/zenodo.24241936](https://doi.org/10.5281/zenodo.24241936).
The dataset record URL and archive checksums are recorded in `release-assets.json`.
The dataset DOI identifies the data, not the software archive.

The synthesis contains 5,962 marker-gene profiles (3,652 16S and 2,310 ITS),
169 plant species and 45 countries. The literature list contains 52 eligible
publications (S1–S53 excluding S6); S54 and S55 are additional dataset sources.
The mapped 16S set contains 22,749 bacterial and 96 archaeal features.

## Reproduce the downstream analysis

1. Obtain the code for the tagged release.
2. Download the three companion data ZIP files from the linked Zenodo dataset.
   Extract **all three into this repository root**. Their members create `data/`,
   `server_result/` and `results/`. They do not contain raw reads.
3. Run `python scripts/tools/verify_assets.py` to verify files and SHA-256 hashes.
4. Prepare R 4.5.2 and Python dependencies as described in `environment/README.md`.
5. Run `Rscript -e 'rmarkdown::render("Analysis_pipeline.Rmd")'`.
6. Run `python scripts/tools/check_outputs.py` and `python scripts/tools/validate_results.py`.

Alternatively, download and extract the published data with:

```sh
python scripts/tools/fetch_data.py --record https://zenodo.org/records/24241936
```

The downloader validates archive checksums and refuses unsafe ZIP paths.
`SEED_THREADS` caps threads at 16. Default runs reuse archived expensive computations.
The 84-row output map includes the original 83 panel/table mappings and Table S7.

## What this release reproduces

This is a downstream analysis release from filtered count tables and prediction
outputs. It regenerates statistical summaries and source panels using supplied
precomputed assets. It does not rerun raw-read retrieval, denoising or functional
prediction. The historical modified iCAMP computation engine is not included;
the pair tables used in the manuscript are included. Environment restoration
has not been validated on a clean machine. See `docs/VALIDATION.md`.

Source plots may retain technical identifiers such as ASV in column names.
These can represent reference-mapped features; representative-host breadth counts
membership in host-specific representative sets, not all detected hosts.
Final manuscript illustrations include manual typography and label adjustments.
Source-panel regeneration is therefore not a pixel-identical artwork rebuild.

Only pooled source-bootstrap stability is used in the manuscript. Archived S/P
bootstrap rows are historical supplementary computations; full-data S/P estimates
provide the reported comparison. See the precomputed-results README.

## Repository layout

- `analysis/`, `R/`: analysis chapters and shared functions.
- `scripts/tools/`: portable download, integrity and validation helpers.
- `scripts/local_maintenance/`: original workstation staging/comparison scripts;
  these are retained for provenance and are not required by readers.
- `environment/`: dependency records and setup notes.
- `manifests/`: figure map and data inventory.
- `data/`, `server_result/`, `results/`: populated by the data archives, excluded from Git.

## Licensing and citation

Original code: MIT. Original processed-data contributions: CC BY 4.0.
Third-party material retains its original terms; see `THIRD_PARTY_NOTICES.md`.
Use `CITATION.cff` for software citation and cite the separate dataset DOI.

New raw sequences: NGDC/GSA PRJCA066385 / CRA044477.
Contact: Ning Ling, nling@njau.edu.cn.
