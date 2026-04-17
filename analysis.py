"""
Single-Cell RNA-seq Analysis Pipeline
======================================
This script walks through a standard scRNA-seq analysis workflow using the
SCP1039 dataset (Wu et al. 2021 — "A single-cell and spatially resolved atlas
of human breast cancers", Nature Genetics).

Dataset: https://singlecell.broadinstitute.org/single_cell/study/SCP1039

The pipeline covers:
    1. Data loading and subsetting
    2. Quality control
    3. Normalisation and feature selection
    4. Dimensionality reduction (PCA → UMAP)
    5. Clustering (Leiden)
    6. Gene expression analysis (YAP1 / WWTR1 co-expression)
    7. Differential expression per cluster
    8. Cell type annotation from metadata

Usage
-----
    python analysis.py --dataset path/to/CID4290A.h5ad \\
                       --metadata path/to/Whole_miniatlas_meta.csv \\
                       --patient_id CID4290A
"""

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns


def run_pipeline(
    dataset_path: str,
    metadata_path: str,
    patient_id: str,
    n_pcs: int = 18,
    leiden_resolution: float = 0.6,
    n_top_genes: int = 2000,
) -> None:

    dataset_name = os.path.basename(dataset_path).split(".")[0]
    output_dir = os.path.join(".", f"{dataset_name}_output")
    os.makedirs(output_dir, exist_ok=True)

    # Tell Scanpy where to save all figures
    sc.settings.figdir = output_dir
    sc.set_figure_params(dpi=150)

    # ==========================================================================
    # 1. DATA LOADING
    # ==========================================================================
    # The metadata CSV maps each cell barcode to its patient ID, cell type,
    # cancer subtype, and other annotations produced by the original study.
    # The first data row contains column type labels ("group", "numeric", ...)
    # rather than real data — these are dropped before use.
    metadata = pd.read_csv(metadata_path, low_memory=False)
    metadata = metadata[metadata["Patient"] != "group"].copy()
    metadata_patient = metadata[metadata["Patient"] == patient_id]

    # The expression data is stored in the AnnData (.h5ad) format, which holds
    # the count matrix, cell barcodes, and gene names in a single file.
    adata = sc.read_h5ad(dataset_path)

    # Each cell barcode is prefixed with the patient ID (e.g. "CID4290A_ACGT..."),
    # so we can isolate one patient's cells from the full atlas.
    adata_patient = adata[adata.obs_names.str.startswith(patient_id), :].copy()
    print(f"Cells loaded for {patient_id}: {adata_patient.n_obs}")

    # ==========================================================================
    # 2. QUALITY CONTROL
    # ==========================================================================
    # Low-quality cells — damaged, dying, or poorly captured — must be
    # identified before analysis. The most informative QC metrics are:
    #
    #   - n_genes_by_counts : number of detected genes per cell
    #       Too few → empty droplet or dead cell
    #       Too many → potential doublet (two cells captured as one)
    #
    #   - total_counts : total UMI counts per cell
    #       Reflects sequencing depth; very low = low-quality cell
    #
    #   - pct_counts_mt : percentage of counts from mitochondrial genes (MT-)
    #       High mitochondrial content indicates a cell whose cytoplasmic
    #       mRNA has leaked out, leaving only the mitochondria intact — a
    #       hallmark of cell damage or death.
    #
    # Thresholds should be chosen by inspecting violin plots and scatter plots
    # of these metrics (see analysis.ipynb for the exploratory plots).
    adata_patient.var["mt"] = adata_patient.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata_patient, qc_vars=["mt"], inplace=True)

    sc.pl.violin(
        adata_patient,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save=f"_{dataset_name}_qc.png",
    )

    # ==========================================================================
    # 3. NORMALISATION AND FEATURE SELECTION
    # ==========================================================================
    # Raw UMI counts vary between cells simply because some cells were sequenced
    # more deeply than others — this is technical variation, not biology.
    # Normalising to a fixed total (10,000 counts per cell, similar to CPM)
    # puts all cells on the same scale so expression values are comparable.
    sc.pp.normalize_total(adata_patient, target_sum=1e4, inplace=True)

    # Gene expression data is highly right-skewed: a handful of highly expressed
    # genes dominate the signal. Log-transforming (log1p = log(x+1), to handle
    # zeros) compresses this dynamic range and makes the data more amenable to
    # linear methods like PCA.
    sc.pp.log1p(adata_patient)

    # Most of the ~30,000 genes in the genome are not informative for
    # distinguishing cell types — they are either not expressed or expressed at
    # the same level in every cell. We keep only the 2,000 most variable genes
    # (those whose variance across cells is highest relative to their mean),
    # which captures the biologically meaningful signal while reducing noise
    # and computation.
    sc.pp.highly_variable_genes(adata_patient, flavor="seurat", n_top_genes=n_top_genes)

    # ==========================================================================
    # 4. DIMENSIONALITY REDUCTION
    # ==========================================================================
    # Even with 2,000 genes, the data still lives in a very high-dimensional
    # space. PCA finds the directions of maximum variance and projects the data
    # onto a much smaller set of principal components (PCs).
    # We use only the highly variable genes for this step.
    sc.pp.pca(adata_patient, svd_solver="arpack", mask_var="highly_variable")

    # Build a k-nearest-neighbour (kNN) graph in PC space. This graph is the
    # shared foundation for both clustering and UMAP visualisation — cells that
    # are transcriptionally similar are connected by edges.
    # n_pcs controls how many PCs are used; this is chosen by inspecting the
    # variance explained by each PC (elbow plot in analysis.ipynb).
    sc.pp.neighbors(adata_patient, n_pcs=n_pcs)

    # UMAP projects the high-dimensional kNN graph down to 2D for visualisation.
    # Important: UMAP preserves local neighbourhood structure (nearby cells in
    # the plot are transcriptionally similar) but distances between clusters are
    # not meaningful — do not interpret the space between clusters.
    sc.tl.umap(adata_patient, n_components=2)
    sc.pl.umap(adata_patient, show=False, save=f"_{dataset_name}_raw.png")

    # ==========================================================================
    # 5. CLUSTERING
    # ==========================================================================
    # Leiden is a community detection algorithm that partitions the kNN graph
    # into groups of highly connected cells (i.e. transcriptionally similar
    # populations). The resolution parameter controls granularity:
    # higher values → more, smaller clusters.
    sc.tl.leiden(
        adata_patient,
        resolution=leiden_resolution,
        key_added="leiden_clusters",
        n_iterations=-1,
        flavor="igraph",
    )
    sc.pl.umap(
        adata_patient,
        color="leiden_clusters",
        title=f"Leiden Clusters — {patient_id}",
        show=False,
        save=f"_{dataset_name}_leiden.png",
    )

    # ==========================================================================
    # 6. GENE EXPRESSION ANALYSIS — YAP1 / WWTR1 CO-EXPRESSION
    # ==========================================================================
    # YAP1 (YAP) and WWTR1 (TAZ) are the downstream effectors of the Hippo
    # signalling pathway and are frequently dysregulated in cancer.
    #
    # We classify each cell into one of four mutually exclusive groups based on
    # whether YAP1 and/or WWTR1 is detected (normalised log expression > 0):
    #
    #   Double_negative : neither YAP1 nor WWTR1 detected
    #   YAP1_only       : only YAP1 detected
    #   WWTR1_only      : only WWTR1 (TAZ) detected
    #   Double_positive : both detected simultaneously
    yap1_expr = adata_patient[:, "YAP1"].X.toarray().flatten()
    wwtr1_expr = adata_patient[:, "WWTR1"].X.toarray().flatten()

    adata_patient.obs["expression_group"] = "Double_negative"
    adata_patient.obs.loc[
        (yap1_expr > 0) & (wwtr1_expr == 0), "expression_group"
    ] = "YAP1_only"
    adata_patient.obs.loc[
        (yap1_expr == 0) & (wwtr1_expr > 0), "expression_group"
    ] = "WWTR1_only"
    adata_patient.obs.loc[
        (yap1_expr > 0) & (wwtr1_expr > 0), "expression_group"
    ] = "Double_positive"

    group_percentages = (
        adata_patient.obs["expression_group"].value_counts(normalize=True) * 100
    )
    group_percentages.to_csv(
        os.path.join(output_dir, f"{dataset_name}_pct_metrics.csv")
    )
    print(f"\nExpression group distribution for {patient_id}:")
    print(group_percentages)

    expression_colors = {
        "Double_negative": "none",
        "YAP1_only": "#1f77b4",
        "WWTR1_only": "#2ca02c",
        "Double_positive": "#d62728",
    }
    sc.pl.umap(
        adata_patient,
        color="expression_group",
        title=f"YAP1/WWTR1 Expression Groups — {patient_id}",
        palette=expression_colors,
        show=False,
        save=f"_{dataset_name}_expression_groups.png",
    )

    # ==========================================================================
    # 7. DIFFERENTIAL EXPRESSION PER CLUSTER
    # ==========================================================================
    # To understand what each Leiden cluster represents biologically, we find
    # genes that are significantly more expressed in one cluster compared to all
    # others (one-vs-rest). This is done using a t-test on the log-normalised
    # expression values.
    #
    # The top marker genes per cluster can then be looked up in databases such
    # as CellMarker or PanglaoDB to manually assign cell type identities.
    sc.tl.rank_genes_groups(adata_patient, groupby="leiden_clusters", method="t-test")
    sc.pl.rank_genes_groups(
        adata_patient,
        n_genes=10,
        title="Top marker genes per cluster",
        show=False,
        save=f"_{dataset_name}.png",
    )

    result = adata_patient.uns["rank_genes_groups"]
    cluster_ids = result["names"].dtype.names
    top_genes_per_cluster = {
        cid: list(result["names"][cid][:10]) for cid in cluster_ids
    }

    top_genes_path = os.path.join(output_dir, f"{dataset_name}_top_genes.txt")
    with open(top_genes_path, "w") as f:
        for cluster_id, genes in top_genes_per_cluster.items():
            f.write(f"Cluster {cluster_id} top features:\n")
            f.write("\n".join(genes))
            f.write("\n\n")

    # ==========================================================================
    # 8. CELL TYPE ANNOTATION
    # ==========================================================================
    # The Wu et al. study provides expert cell type annotations for every cell
    # in the atlas (stored in the metadata CSV). We map these labels back onto
    # our cells using the barcode as the key. This gives us two levels of
    # resolution:
    #   - celltype_major : broad categories (e.g. Cancer Epithelial, T-cells)
    #   - celltype_minor : finer subtypes (e.g. Cancer Her2 SC, T cells CD8+)
    cell_type_major = metadata_patient.set_index("NAME")["celltype_major"].to_dict()
    cell_type_minor = metadata_patient.set_index("NAME")["celltype_minor"].to_dict()
    adata_patient.obs["celltype_major"] = adata_patient.obs.index.map(cell_type_major)
    adata_patient.obs["celltype_minor"] = adata_patient.obs.index.map(cell_type_minor)

    sc.pl.umap(
        adata_patient,
        color="celltype_major",
        title=f"Major Cell Types — {patient_id}",
        show=False,
        save=f"_{dataset_name}_celltype_major.png",
    )
    sc.pl.umap(
        adata_patient,
        color="celltype_minor",
        title=f"Minor Cell Types — {patient_id}",
        show=False,
        save=f"_{dataset_name}_celltype_minor.png",
    )

    # Heatmap showing what percentage of total cells fall into each
    # (cell type × expression group) combination. This reveals whether
    # YAP1/WWTR1 expression is enriched in specific cell populations.
    n_cells = len(adata_patient)

    for level, col in [("major", "celltype_major"), ("minor", "celltype_minor")]:
        cross_tab = pd.crosstab(
            adata_patient.obs[col], adata_patient.obs["expression_group"]
        )
        cross_tab["Total"] = cross_tab.sum(axis=1)
        cross_tab.loc["Total"] = cross_tab.sum()
        cross_tab_pct = (cross_tab / n_cells) * 100

        plt.figure(figsize=(12, 8))
        sns.heatmap(
            cross_tab_pct,
            annot=True,
            fmt=".2f",
            cmap="YlOrRd",
            cbar_kws={"label": "Percentage of Cells"},
        )
        plt.title(
            f"Distribution of Expression Groups across Cell Types"
            f" normalised by total cell count — {patient_id} ({level})"
        )
        plt.xlabel("Expression Group")
        plt.ylabel("Cell Type")
        plt.tight_layout()
        plt.savefig(
            os.path.join(output_dir, f"{dataset_name}_heatmap_{level}.png")
        )
        plt.close()

    print(f"\nPipeline complete. All outputs saved to: {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "scRNA-seq analysis pipeline for the SCP1039 breast cancer atlas. "
            "Reproduces the workflow from analysis.ipynb."
        )
    )
    parser.add_argument("--dataset", required=True, help="Path to .h5ad file")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV")
    parser.add_argument(
        "--patient_id", required=True, help="Patient ID prefix (e.g. CID4290A)"
    )
    parser.add_argument(
        "--n_pcs", type=int, default=18,
        help="Number of PCs for neighbour graph (default: 18)",
    )
    parser.add_argument(
        "--leiden_resolution", type=float, default=0.6,
        help="Leiden clustering resolution (default: 0.6)",
    )
    parser.add_argument(
        "--n_top_genes", type=int, default=2000,
        help="Number of highly variable genes for PCA (default: 2000)",
    )
    args = parser.parse_args()

    run_pipeline(
        dataset_path=args.dataset,
        metadata_path=args.metadata,
        patient_id=args.patient_id,
        n_pcs=args.n_pcs,
        leiden_resolution=args.leiden_resolution,
        n_top_genes=args.n_top_genes,
    )
