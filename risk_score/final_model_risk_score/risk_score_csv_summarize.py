import os
import numpy as np
import pandas as pd
import scanpy as sc

# ==========================================
# 1. Load Data
# ==========================================
# Replace these with your actual file names/paths
adata_path = "GSE176206_msc_sokm.h5ad"
risk_genes_path = "finalized_risk_gene_list.csv"
output_h5ad_path = "cleared_risk_score.h5ad"
output_csv_path = "GSE176206_msc_sokm_risk_donor_summaries.csv"

adata = sc.read_h5ad(adata_path)
df_risk = pd.read_csv(risk_genes_path)

# Normalization and log1p transformation
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ==========================================
# 2. Define and Filter Risk Genes
# ==========================================
pluripotency_genes = df_risk['Pluripotency'].dropna().tolist()
genomic_stress_genes = df_risk['Genomic_Stress'].dropna().tolist()

# Retain only genes present in the dataset
pluripotency_genes_in_data = [g for g in pluripotency_genes if g in adata.var_names]
genomic_stress_genes_in_data = [g for g in genomic_stress_genes if g in adata.var_names]

print(f"Pluripotency genes in data: {len(pluripotency_genes_in_data)} / {len(pluripotency_genes)}")
print(f"Genomic Stress genes in data: {len(genomic_stress_genes_in_data)} / {len(genomic_stress_genes)}")

# ==========================================
# 3. Calculate Single-Cell Risk Scores
# ==========================================
print("\nCalculating Risk Scores...")

# Pluripotency Risk Score
sc.tl.score_genes(adata, gene_list=pluripotency_genes_in_data, ctrl_size=50, n_bins=25, score_name='Pluripotency_Risk_Score', use_raw=False, random_state=42)

# Genomic Stress Risk Score
sc.tl.score_genes(adata, gene_list=genomic_stress_genes_in_data, ctrl_size=50, n_bins=25, score_name='Genomic_Stress_Risk_Score', use_raw=False, random_state=42)

# Reverse Pseudotime for intuitive trajectory interpretation (0 = Somatic, 1 = Pluripotent)
if 'velocity_pseudotime' in adata.obs:
    adata.obs['Reprogramming_Trajectory'] = 1.0 - adata.obs['velocity_pseudotime']

# Save the scored AnnData object
adata.write_h5ad(output_h5ad_path)
print("Scoring completed and AnnData saved successfully.")

# ==========================================
# 4. Donor-Level Aggregation (CSV Output)
# ==========================================
# Calculate gene coverage
coverage_pluri = len(pluripotency_genes_in_data) / len(pluripotency_genes) if pluripotency_genes else 0
coverage_tumor = len(genomic_stress_genes_in_data) / len(genomic_stress_genes) if genomic_stress_genes else 0

df_obs = adata.obs.copy()

# Metadata column names (Adjust these if your metadata columns differ)
col_age = 'age'            
col_treatment = 'treatment'
col_animal = 'animal'   

risk_axes = {
    'pluripotency': 'Pluripotency_Risk_Score',
    'genomic_stress': 'Genomic_Stress_Risk_Score'
}

# Calculate baseline thresholds (90th, 95th, 99th percentiles) from NegCtrl for each age group
thresholds = {}
for age_group in df_obs[col_age].unique():
    ctrl_df = df_obs[(df_obs[col_treatment] == 'NegCtrl') & (df_obs[col_age] == age_group)]
    
    thresholds[age_group] = {}
    for axis_name, score_col in risk_axes.items():
        scores = ctrl_df[score_col].dropna()
        if len(scores) > 0:
            thresholds[age_group][axis_name] = {
                'p90': np.percentile(scores, 90),
                'p95': np.percentile(scores, 95),
                'p99': np.percentile(scores, 99)
            }

results = []

# Group by age, treatment condition, and animal ID
for (age_val, trt_val, animal_val), group in df_obs.groupby([col_age, col_treatment, col_animal]):
    
    age_group = str(age_val).lower()
    
    # Create canonical biological unit ID
    unit_id = f"GSE176206_MSC|{age_group}|{trt_val}|animal_{animal_val}"
    
    # Define analysis role
    role = 'control' if trt_val == 'NegCtrl' else 'SOKM'
    
    row_data = {
        'biological_unit_id': unit_id,
        'dataset_id': 'GSE176206_msc_sokm',
        'age_group': age_group,
        'exact_treatment_arm': trt_val,
        'analysis_role': role,
        'animal_label': animal_val,
        'n_cells': len(group),
        'qc_status': 'pass' if len(group) >= 50 else 'fail_low_cells',
        'pluripotency_gene_coverage': coverage_pluri,
        'genomic_stress_gene_coverage': coverage_tumor
    }
    
    # Calculate metrics for each risk axis
    for axis_name, score_col in risk_axes.items():
        scores = group[score_col].dropna()
        
        if len(scores) == 0:
            continue
            
        # Central tendency & Secondary metrics
        row_data[f'{axis_name}_median'] = np.median(scores)
        row_data[f'{axis_name}_mean'] = np.mean(scores)
        
        # Upper-tail risk
        row_data[f'{axis_name}_p95'] = np.percentile(scores, 95)
        
        # High-risk cell fraction (based on age-matched control thresholds)
        if age_val in thresholds and axis_name in thresholds[age_val]:
            t_90 = thresholds[age_val][axis_name]['p90']
            t_95 = thresholds[age_val][axis_name]['p95']
            t_99 = thresholds[age_val][axis_name]['p99']
            
            row_data[f'{axis_name}_high_risk_frac_p90'] = (scores > t_90).mean()
            row_data[f'{axis_name}_high_risk_frac_p95'] = (scores > t_95).mean() # Primary metric
            row_data[f'{axis_name}_high_risk_frac_p99'] = (scores > t_99).mean()
            
    results.append(row_data)

# Save the aggregated results to CSV
df_summary = pd.DataFrame(results)
df_summary.to_csv(output_csv_path, index=False)

print(f"\nDonor-level Risk Output Table saved to: {output_csv_path}")
print(f"Total rows exported: {len(df_summary)}")