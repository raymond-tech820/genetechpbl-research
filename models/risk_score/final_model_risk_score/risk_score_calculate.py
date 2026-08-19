import pandas as pd
import scanpy as sc
import os

# 1. Load data
adata_path = "GSE176206_msc_sokm.h5ad"
risk_genes_path = "finalized_risk_gene_list.csv"

adata = sc.read_h5ad(adata_path)
df_risk = pd.read_csv(risk_genes_path)

# Uncomment if normalization and log1p are required for the loaded raw data
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. Define risk genes
pluripotency_genes = df_risk['Pluripotency'].dropna().tolist()
genomic_stress_genes = df_risk['Genomic_Stress'].dropna().tolist()

# Filter for genes present in the dataset
pluripotency_genes_in_data = [g for g in pluripotency_genes if g in adata.var_names]
genomic_stress_genes_in_data = [g for g in genomic_stress_genes if g in adata.var_names]

print(f"Pluripotency genes in data: {len(pluripotency_genes_in_data)} / {len(pluripotency_genes)}")
print(f"Genomic Stress genes in data: {len(genomic_stress_genes_in_data)} / {len(genomic_stress_genes)}")

# 3. Calculate Risk Scores
print("\nCalculating Risk Scores...")

# Pluripotency Risk Score
sc.tl.score_genes(adata, gene_list=pluripotency_genes_in_data, ctrl_size=50, n_bins=25, score_name='Pluripotency_Risk_Score', use_raw=False, random_state=42)

# Genomic Stress Risk Score
sc.tl.score_genes(adata, gene_list=genomic_stress_genes_in_data, ctrl_size=50, n_bins=25, score_name='Genomic_Stress_Risk_Score', use_raw=False, random_state=42)

# 4. Reverse Pseudotime for trajectory interpretation
adata.obs['Reprogramming_Trajectory'] = 1.0 - adata.obs['velocity_pseudotime']

# 5. Save the scored AnnData object
output_path = "risk_score_calulated.h5ad"
adata.write_h5ad(output_path)

print("Scoring completed successfully.")