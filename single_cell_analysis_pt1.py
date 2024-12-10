import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def process_dataset(metadata_path, matrix_path, barcodes_path, features_path, output_h5ad_path):
    
    metadata = pd.read_csv(metadata_path)
    
    # Extract unique patient IDs
    patient_ids = metadata['Patient'].unique()
    
    for patient_id in patient_ids:
        # Filter metadata for the current patient
        metadata_patient = metadata[metadata['Patient'] == patient_id]
        
        # Load expression data
        ##This is done in the case where both the barcodes and features are placed in separate files
        ##These files are then put into a single anndata folder
        adata = sc.read_mtx(matrix_path)
        barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')
        features = pd.read_csv(features_path, header=None, sep='\t')
        
        
        # Set observation and variable names
        adata.obs_names = features[0]
        adata.var_names = barcodes[0]
        adata = adata.T
        
        # Filter data for the current patient
        adata_patient = adata[adata.obs_names.str.startswith(patient_id), :]
        
        # Save and reload the filtered data
        adata_patient.write(output_h5ad_path.format(patient_id=patient_id))
        adata_patient = sc.read_h5ad(output_h5ad_path.format(patient_id=patient_id))
        
        # Calculate mitochondrial gene metrics
        adata_patient.var['mt'] = adata_patient.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata_patient, qc_vars=['mt'], inplace=True)
        
        # Plot QC metrics
        sc.pl.violin(adata_patient, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
        
        # Filter cells based on QC metrics
        #The examples here are 5,800,6800
        bdata_patient = adata_patient.copy()
        bdata_patient = bdata_patient[bdata_patient.obs['pct_counts_mt'] < 5]
        bdata_patient = bdata_patient[bdata_patient.obs['total_counts'] > 800]
        bdata_patient = bdata_patient[bdata_patient.obs['n_genes_by_counts'] < 6800]
        
