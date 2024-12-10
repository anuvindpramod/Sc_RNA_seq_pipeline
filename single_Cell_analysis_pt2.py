import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import leidenalg 
def analyze_dataset(preprocessed_h5ad_path, metadata_path):
        bdata_patient=sc.read_h5ad(preprocessed_h5ad_path)
        sc.pp.normalize_total(bdata_patient, target_sum=1e4, inplace=True)
        sc.pp.log1p(bdata_patient)
        
        # Identify highly variable genes
        sc.pp.highly_variable_genes(bdata_patient, flavor='seurat', n_top_genes=2000)
        sc.pl.highly_variable_genes(bdata_patient)
        
        # Perform PCA and UMAP
        sc.pp.pca(bdata_patient, svd_solver='arpack')
        sc.pp.neighbors(bdata_patient, n_pcs=22)
        sc.tl.umap(bdata_patient)
        sc.pl.umap(bdata_patient)
        
        # Cluster the cells
        sc.tl.leiden(bdata_patient, resolution=0.6, key_added='leiden_clusters')
        sc.pl.umap(bdata_patient, color='leiden_clusters', title=f'Leiden Clusters of {patient_id}')
        
        # Rank genes for each cluster
        sc.tl.rank_genes_groups(bdata_patient, groupby='leiden_clusters', method='t-test')
        sc.pl.rank_genes_groups(bdata_patient, n_genes=10, title='Top Genes per Cluster')
        
        # Map cell types from metadata
        cell_type_dict = metadata.set_index('NAME')['celltype_major'].to_dict()
        bdata_patient.obs['celltype_major'] = bdata_patient.obs.index.map(cell_type_dict)
        
        # Plot UMAP with cell types
        sc.pl.umap(bdata_patient, color='celltype_major', title=f'Cell Type in {patient_id}')
        
        # Create a cross-tabulation of cell types and expression groups
        cross_tab = pd.crosstab(bdata_patient.obs['celltype_major'], bdata_patient.obs['expression_group'])
        cross_tab['Total'] = cross_tab.sum(axis=1)
        cross_tab.loc['Total'] = cross_tab.sum()
        cross_tab_normalized = (cross_tab / len(bdata_patient.obs)) * 100
        
        # Plot heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(cross_tab_normalized, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Percentage of Cells'})
        plt.title(f'Distribution of Expression Groups across Cell Types for {patient_id}')
        plt.xlabel('Expression Group')
        plt.ylabel('Cell Type')
        plt.tight_layout()
        plt.show()
