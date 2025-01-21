import os
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def run_pipeline(dataset_path, metadata_path):
    dataset_name = os.path.basename(dataset_path).split('.')[0]
    output_dir = f"./{dataset_name}_output"
    os.makedirs(output_dir, exist_ok=True)
    #data loading
    adata = sc.read_h5ad(dataset_path)
    metadata = pd.read_csv(metadata_path)
    group = 'CID4290A'
    metadata_group = metadata[metadata['Patient'] == group]

    # Subset data
    adata_group = adata[adata.obs_names.str.startswith(group), :].copy()

    # processing
    sc.pp.normalize_total(adata_group, target_sum=1e4)
    sc.pp.log1p(adata_group)
    sc.pp.highly_variable_genes(adata_group, flavor='seurat', n_top_genes=2000)
    sc.pp.pca(adata_group, svd_solver='arpack',mask_var='highly_variable')
    sc.pp.neighbors(adata_group, n_pcs=20)
    sc.tl.umap(adata_group,n_components=2)
    sc.tl.leiden(adata_group, resolution=0.6, key_added='leiden_clusters')
    yap1_exp_sig=adata_group[:,'YAP1'].X.toarray().flatten()
    wwtr1_exp_sig=adata_group[:,'WWTR1'].X.toarray().flatten()
    adata_group.obs['expression_group']='Double_negative'
    adata_group.obs.loc[(yap1_exp_sig>0)&(wwtr1_exp_sig==0),'expression_group']="YAP1_only"
    adata_group.obs.loc[(yap1_exp_sig==0)&(wwtr1_exp_sig>0),'expression_group']="WWTR1_only"
    adata_group.obs.loc[(yap1_exp_sig>0)&(wwtr1_exp_sig>0),'expression_group']="Double_postive"
    group_percentages = adata_group.obs['expression_group'].value_counts(normalize=True) * 100
    group_percentages.to_csv('{dataset_name}_pct_metrics.txt')

    # UMAP plot
    colors={'Double_negative': 'none', #invisible
    'YAP1_only': '#1f77b4',       #blue
    'WWTR1_only': '#2ca02c',      #green
    'Double_postive': '#d62728'    #red
    }
    sc.pl.umap(adata_group,save=f'/{dataset_name}raw_cluster')
    sc.pl.umap(adata_group,color='expression_group',title=f'Expression group of {group}',palette=colors,save=f'/{dataset_name}_expression_clusters.png')
    sc.pl.umap(adata_group, color='leiden_clusters', title='Leiden Clusters', save=f'/{dataset_name}_leiden_clusters.png')

    # Ranking the genes
    sc.tl.rank_genes_groups(adata_group, groupby='leiden_clusters', method='Wilcoxon')
    '''sc.pl.rank_genes_groups(adata_group, n_genes=10, title='Ranked Genes', save=f'/{dataset_name}_ranked_genes.png')
    '''
    result = adata_group.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    top_features = {}
    n_top_genes = 10  # desired number of top genes per cluster
    for group in groups:
        top_features[group] = result["names"][group][:n_top_genes] 
    # extracting top features for each cluster
    with open(f"{output_dir}/{dataset_name}_top_features.txt", "w") as file:
        for group, features in top_features.items():
            file.write(f"Cluster {group} top features:\n")
            for feature in features:
                file.write(f"{feature}\n")
            file.write("\n")
    # Adding the cell type information
    cell_type_dict = metadata_group.set_index('NAME')['celltype_major'].to_dict()
    adata_group.obs['celltype_major'] = adata_group.obs.index.map(cell_type_dict)

    # Ploting cell types
    sc.pl.umap(adata_group, color='celltype_major', title='Cell Type', save=f'/{dataset_name}_cell_type.png')

    #heatmap
    cross_tab = pd.crosstab(adata_group.obs['celltype_major'], adata_group.obs['leiden_clusters'])
    cross_tab_normalized = (cross_tab / len(adata_group)) * 100
    plt.figure(figsize=(12, 8))
    sns.heatmap(cross_tab_normalized, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Percentage of Cells'})
    plt.title('Distribution of Clusters across Cell Types')
    plt.xlabel('Leiden Clusters')
    plt.ylabel('Cell Type')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{dataset_name}_heatmap.png")
    plt.close()

    print(f"Pipeline completed. Outputs saved in {output_dir}")


run_pipeline('/path/to/dataset.h5ad', '/path/to/metadata.csv')
