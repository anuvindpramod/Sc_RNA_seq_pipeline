# scRNA-seq YAP1/WWTR1 co-expression analysis

This project looks at how YAP (YAP1) and TAZ (WWTR1) — the two main effectors of the Hippo signalling pathway — are expressed across different cell types in breast cancer tissue. The analysis uses the SCP1039 dataset from Wu et al. 2021 (*A single-cell and spatially resolved atlas of human breast cancers*, Nature Genetics), which covers 26 primary breast tumours with expert cell type annotations.

The exploratory analysis lives in `analysis.ipynb`. The `analysis.py` script runs the same pipeline end-to-end from the command line, which is useful if you want to apply it to multiple patients without stepping through the notebook each time.

## Getting the data

The dataset is publicly available on the Broad Single Cell Portal:
https://singlecell.broadinstitute.org/single_cell/study/SCP1039

You need two things from there:
- The per-patient `.h5ad` expression files (under the Downloads tab)
- `Whole_miniatlas_meta.csv` — the cell-level metadata file with cell type annotations

## Setup

```bash
pip install -r requirements.txt
```

## Running the pipeline

```bash
python analysis.py \
  --dataset path/to/CID4290A.h5ad \
  --metadata path/to/Whole_miniatlas_meta.csv \
  --patient_id CID4290A
```

All outputs are written to a folder called `<dataset_name>_output/` in the current directory. This includes the UMAP plots, cluster marker genes, expression group percentages, and the cell type heatmaps.

The `--patient_id` flag should match the barcode prefix in the `.h5ad` file (e.g. `CID4290A`, `CID44991`). See the metadata CSV for the full list of patient IDs.

Other parameters and their defaults:

| Flag | Default | What it controls |
|------|---------|-----------------|
| `--n_pcs` | 18 | PCs used to build the neighbour graph |
| `--leiden_resolution` | 0.6 | Clustering granularity (higher = more clusters) |
| `--n_top_genes` | 2000 | Highly variable genes passed to PCA |
