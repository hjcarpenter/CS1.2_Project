# -*- coding: utf-8 -*-
"""
AMR Exploratory Analysis
This script:
Imports the cleaned hAMRonization csv.
Identifies Acquired ARGs (resfinder hits) and Protein variants 
Generates exploratory plots (counts/ proportions)
Identifies genes with low coverage/ % id (possible non functional) and duplicates
Generates a resistance summary dataframe (linking gene to resistance)
Extracts unique antimicrobial drug classes and agents. 
Generates a drug class heatmap plot per contig/ resistance type

Possible Ouputs:
Data:
AMR_agent_class_summary.csv – Contains summary statistics (total hits, unique agents, unique classes) overall and per contig.
Total_AMR_determinants_per_contig_summary.csv – A breakdown of determinant counts per contig, including the percentages of acquired ARGs and protein variants.
Genes_max_coverage_lt100.csv – A subset of genes where the maximum coverage percentage is less than 100%.
Genes_max_id_lt90.csv – A subset of genes where the maximum identity percentage is less than 90%.
clean_AMR_data_subset.csv – The cleaned, filtered subset of the original hAMRonization data containing only the specified columns of interest.
Resistance_df.csv – A deduplicated dataframe focusing specifically on resistance mechanisms, including boolean flags for acquired, protein variant, and intrinsic resistance.
all_antimicrobials_table.csv – A summary frequency table of individual antimicrobial agents and their associated unique genes.
all_drug_classes_table.csv – A summary frequency table of individual drug classes and their associated unique genes.
Table_Acquired_ARGs.csv – A specific sub-table detailing only the acquired antibiotic resistance genes.
Table_Protein_Variants.csv – A specific sub-table detailing only the protein variants conferring resistance.
Table_Intrinsic_Genes.csv – A specific sub-table detailing only the intrinsic (core chromosomal) resistance genes.
Master_AMR_Determinants_Supplemental.csv – A merged master table combining the acquired, protein variant, and intrinsic lists, with an added Resistance_Type column.
Visualizations:
AMR_proportions_by_contig_piechart.png – A pie chart showing the proportion of total AMR determinants found on each contig.
AMR_determinants_per_contig.svg – A grouped bar chart displaying the absolute counts of intrinsic, acquired, and protein variant resistance per contig.
AMR_determinants_per_contig_proportional.svg – A stacked bar chart showing the proportional breakdown of resistance types per contig.
AMR_Functional_Heatmap_Comparison.svg – Test version of a 3-panel heatmap showing the presence/absence of drug classes across contigs separated by resistance type (Intrinsic, Acquired, Protein Variant).
AMR_Heatmap_With_Legend.svg – A highly refined version of the 3-panel heatmap featuring custom Okabe-Ito colormaps, row totals, and a custom legend.

Created on Fri Nov 28 03:52:00 2025
@author: Heath
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

# ==========================================
# --- HELPER FUNCTIONS ---
# ==========================================

def count_true_uniques(series):
    """Safely splits delimited strings and counts the unique individual items."""
    return (
        series
        .dropna()               
        .astype(str)            
        .str.split(';')         
        .explode()              
        .str.strip()            
        .replace('', pd.NA)     
        .nunique()              
    )

def get_summary_stats(df, prefix):
    overall_stats = pd.DataFrame({
        f'{prefix}_total_hits': [len(df)],
        f'{prefix}_unique_agents': [count_true_uniques(df['antimicrobial_agent'])],
        f'{prefix}_unique_classes': [count_true_uniques(df['drug_class'])]
    }, index=['OVERALL_TOTAL'])
    
    contig_stats = df.groupby('contig').agg(
        total_hits=('contig', 'size'),
        unique_agents=('antimicrobial_agent', count_true_uniques),
        unique_classes=('drug_class', count_true_uniques)
    ).rename(columns={
        'total_hits': f'{prefix}_total_hits',
        'unique_agents': f'{prefix}_unique_agents',
        'unique_classes': f'{prefix}_unique_classes'
    })
    
    return overall_stats, contig_stats

def max_from_semicolon(val):
    """Helper function for semicolon numbers to extract maximum float."""
    if pd.isna(val): return np.nan
    if isinstance(val, (int, float, np.number)): return float(val)
    parts = str(val).split(';')
    floats = []
    for p in parts:
        p = p.strip()
        if p:
            try:
                floats.append(float(p))
            except ValueError:
                pass 
    return max(floats) if floats else np.nan

def get_unique_items(series):
    """Extracts a sorted list of individual unique items from a delimited series."""
    unique_array = (
        series
        .dropna()
        .astype(str)
        .str.split(';')
        .explode()
        .str.strip()
        .replace('', pd.NA)
        .dropna()               
        .unique()               
    )
    return sorted(unique_array) 

def create_summary_table(df, column):
    """Creates a frequency table of individual items and lists their associated unique genes."""
    # 1. Create a temporary subset of just the column we need and 'gene'
    temp_df = df[[column, 'gene']].dropna(subset=[column]).copy()
    
    # 2. Explode the target column 
    temp_df[column] = temp_df[column].astype(str).str.split(';')
    temp_df = temp_df.explode(column)
    temp_df[column] = temp_df[column].str.strip()
    
    # Drop empties that resulted from the split
    temp_df = temp_df.replace('', pd.NA).dropna(subset=[column])
    # drop duplicates
    temp_df = temp_df.drop_duplicates()
    
    # 3. Group by the exploded items and aggregate both counts and unique genes
    summary_df = (
        temp_df.groupby(column)
        .agg(
            count=(column, 'size'), # Count total occurrences
            # Grab unique genes, sort them alphabetically, and join with a semicolon
            gene=('gene', lambda x: '; '.join(sorted(x.dropna().astype(str).unique())))
        )
        .reset_index()
        .sort_values(by='count', ascending=False) # Keep it sorted by highest count
        .reset_index(drop=True)
    )
    
    return summary_df

 
def create_presence_matrix(df, flag_col, drug_col='drug_class'):
    """Creates a binary 1/0 matrix for contig vs drug_class, excluding NaNs
    to use for the heatmaps."""
    # Filter for the specific resistance type
    filtered = df[df[flag_col] == True].copy()
    
    # 1. Remove true nulls before processing
    filtered = filtered.dropna(subset=[drug_col])
    
    # 2. Split by semicolon and explode
    filtered[drug_col] = filtered[drug_col].astype(str).str.split(';')
    exploded = filtered.explode(drug_col)
    
    # 3. Clean up whitespace and remove "nan" or empty strings
    exploded[drug_col] = exploded[drug_col].str.strip()
    exploded = exploded[~exploded[drug_col].isin(['nan', 'NaN', 'None', '', 'NA'])]
    
    if exploded.empty:
        return pd.DataFrame()
    
    # Create the pivot table (1 if present, 0 if absent)
    matrix = pd.crosstab(exploded['contig'], exploded[drug_col])
    matrix = (matrix > 0).astype(int)
    return matrix

# ==========================================
# --- 1. IMPORTS AND LOADING ---
# ==========================================
fp = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/hAMRonize/hAMR_Result_Clean.csv'

keep = ['contig', 'antimicrobial_agent', 'drug_class', 'resistance_mechanism',
        'start', 'end', 'gene', 'gene_name', 'genetic_variation_type', 
        'coverage_percentage', 'amino_acid_mutation', 'identity', 'strand', 
        'analysis_software_name', 'reference_database_name']

df = pd.read_csv(fp)

# Check all keep columns exist
missing_cols = [c for c in keep if c not in df.columns]
if missing_cols:
    print(f"Warning: The following columns are missing from CSV: {missing_cols}")
    keep = [c for c in keep if c in df.columns]

subset = df[keep].copy()

# Rename 'plasmid1' for plotting
subset['contig'] = subset['contig'].replace('plasmid1', 'pCS1.2IncF-NDM')

# Simplify variation type
subset['genetic_variation_type'] = (
    subset['genetic_variation_type']
    .replace({'gene_presence_detected': 'resistance gene', 
              'protein_variant_detected': 'protein variant'})
)

# ---- Flag Acquired ARGs (ResFinder hits + known mobile passengers)
subset["acq_arg"] = False

# identify hits from the ResFinder database (definite acq)
m_resfinder = subset["reference_database_name"].astype("string").str.contains(
    "resfinder", case=False, na=False
)
subset.loc[m_resfinder, "acq_arg"] = True

# add known mobile 'passenger' genes missed by the ResFinder flag
# mobile_passengers = ['qacEdelta1', 'ble', 'mrx',
#                      'klebsiella pneumoniae kpnE','klebsiella pneumoniae kpnF' ]
# subset.loc[subset['gene'].isin(mobile_passengers), "acq_arg"] = True
mobile_passengers = ['qacedelta1', 'ble', 'mrx', 'kpn'] # lowercase
subset.loc[subset['gene'].str.lower().str.contains('|'.join(mobile_passengers)), "acq_arg"] = True

# Create the hits dataframe (acquired args) for downstream analysis
resfinder_hits = subset.loc[subset['acq_arg']].copy()

# --- Flag protein variants
subset["is_protein_variant"] = (subset["genetic_variation_type"] == "protein variant")
protein_vars = subset.loc[subset['is_protein_variant']].copy()

# ---- check what genes are falling into the 'intrinsic category'
# edit acquired ARG flagging as required
intrinsic_list = subset[(subset['acq_arg'] == False) & (subset['is_protein_variant'] == False)]
print(intrinsic_list['gene'].unique())

# ==========================================
# --- 2. EXPLORATORY ANALYSIS ---
# ==========================================
print(f'total Resfinder hits = {len(resfinder_hits)}')
print(f'unique resfinder hits = {resfinder_hits.nunique()}')
print(f'total protein variants = {len(protein_vars)}')
print(f'unique protein vars = {protein_vars["gene"].nunique()}')

# Generate Summary Stats Dataframes
args_overall, args_contig = get_summary_stats(resfinder_hits, 'args')
prot_overall, prot_contig = get_summary_stats(protein_vars, 'prot')
all_overall, all_contig = get_summary_stats(subset, 'all_amr')

df_overall = pd.concat([args_overall, prot_overall, all_overall], axis=1)
df_contig_summary = pd.concat([args_contig, prot_contig, all_contig], axis=1).fillna(0).astype(int)
df_summary_table = pd.concat([df_overall, df_contig_summary])
df_summary_table.to_csv('AMR_agent_class_summary.csv', index = True)

# Per-Contig Summary
contig_counts = subset['contig'].value_counts().rename('det_count')
total_dets = contig_counts.sum()
print(f'\ntotal AMR determinants (incl dupes): {total_dets}')

contig_prop = (contig_counts / total_dets).rename('proportion')
contig_summary = pd.concat(
    [contig_counts, contig_prop], axis=1).reset_index(
        ).rename(columns={'index': 'contig'})


# ==========================================
# --- 3. PLOTTING ---
# ==========================================

# A. Quick pie - total determinants/contig 
colors = sns.color_palette('viridis', n_colors=len(contig_counts))
plt.figure(figsize=(7,7))
plt.pie(
    contig_counts,
    autopct='%1.1f%%',
    startangle=90,
    colors=colors,
    labels=None 
)
plt.title("Total AMR Determinants per Contig")
plt.legend(contig_counts.index, bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("AMR_proportions_by_contig_piechart.png", dpi=300, bbox_inches='tight') 
plt.show()

# base table for bar charts
contig_summary = (
    subset.groupby("contig")
    .agg(
        det_count=("contig", "size"),
        acq_arg_count=("acq_arg", "sum"),
        protein_variant_count=("is_protein_variant", "sum")
    )
).sort_index()

contigs = contig_summary.index.astype(str)
x = np.arange(len(contigs))
width = 0.2

# B. Grouped bar chart absolute Ns
plt.figure(figsize=(6, 6))
plt.bar(x - width, (contig_summary["det_count"]-contig_summary["acq_arg_count"]-contig_summary["protein_variant_count"]), width, label="Intrinsic Resistance", color="#0072B2")
plt.bar(x, contig_summary["acq_arg_count"], width, label="Acquired Antibiotic Resistance Genes", color="#E69F00")
plt.bar(x + width, contig_summary["protein_variant_count"], width, label="Protein variants conferring Resistance", color="#009E73")

plt.ylabel("N AMR determinants")
plt.xticks(x, contigs, rotation=0)
ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.set_ylim(top=45)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(direction="out")
plt.legend(fontsize=8, frameon=False)
plt.tight_layout()
plt.savefig("AMR_determinants_per_contig.svg", bbox_inches='tight')
plt.show()

# C. Proportional breakdown per contig
prop = contig_summary.copy()
prop["acq_only"] = prop["acq_arg_count"]
prop["protein_variant"] = prop["protein_variant_count"]
prop["other"] = prop["det_count"] - prop["acq_arg_count"] - prop["protein_variant_count"]
prop_frac = prop[["acq_only", "protein_variant", "other"]].div(prop["det_count"], axis=0)

label_map = {"acq_only": "Acquired ARG", "protein_variant": "Protein variant", "other": "Intrinsic Resistance"}
colors_map = {"acq_only": "#E69F00", "protein_variant": "#009E73", "other": "#0072B2"}

plt.figure(figsize=(6, 6))
bottom = None
for col in ["acq_only", "protein_variant", "other"]:
    plt.bar(x, prop_frac[col], bottom=bottom, label=label_map[col], color=colors_map[col])
    bottom = prop_frac[col] if bottom is None else bottom + prop_frac[col]

plt.ylabel("Proportion of AMR determinants")
plt.xticks(x, contigs, rotation=0)
ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.set_ylim(top=1.0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(direction="out")
plt.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
plt.tight_layout()
plt.savefig("AMR_determinants_per_contig_proportional.svg", bbox_inches='tight') # Added bbox_inches
plt.show()

# --- Save contig summary table
contig_summary["pct_acquired_ARG"] = (contig_summary["acq_arg_count"] / contig_summary["det_count"] * 100)
contig_summary["pct_protein_var"] = (contig_summary["protein_variant_count"] / contig_summary["det_count"] * 100)
contig_summary.to_csv('Total_AMR_determinants_per_contig_summary.csv', index=False)


# ==========================================
# --- 4. FUNCTIONAL CHECKS ---
# ==========================================
gene_contig_counts = subset.groupby("gene")["contig"].nunique().sort_values(ascending=False)
genes_multi_contig = gene_contig_counts[gene_contig_counts > 1]
print("\n--- Genes found on >1 contig ---")
print(genes_multi_contig)

subset['max_coverage'] = subset['coverage_percentage'].apply(max_from_semicolon)
subset['max_identity'] = subset['identity'].apply(max_from_semicolon)

subset_cov_lt100 = subset[subset['max_coverage'] < 100].copy()
subset_id_lt90 = subset[subset['max_identity'] < 90].copy()
subset_cov_lt100.to_csv('Genes_max_coverage_lt100.csv', index = False)
subset_id_lt90.to_csv('Genes_max_id_lt90.csv', index = False)

gene_counts = subset['gene'].value_counts()
multi_copy_genes = gene_counts[gene_counts > 1] 

cov_multi = subset_cov_lt100[subset_cov_lt100['gene'].isin(multi_copy_genes.index)]
id_multi = subset_id_lt90[subset_id_lt90['gene'].isin(multi_copy_genes.index)]

print("\nGenes <100% coverage that appear in duplicate in full dataset:")
print(cov_multi[['gene','contig','max_coverage']])
print("\nGenes <90% identity that appear in duplicate in full dataset:")
print(id_multi[['gene','contig','max_identity']])


# ==========================================
# --- 5. RESISTANCE DFS & SUMMARIES ---
# ==========================================
RDF_cols = ['gene','antimicrobial_agent', 'drug_class',
            'resistance_mechanism', 'gene_name', 'contig', 
            'genetic_variation_type', 'acq_arg']
res_df = subset[RDF_cols].drop_duplicates().reset_index(drop=True)
# Add a 'protein_var' flag
# 'genetic_variation_type' == 'protein variant'
res_df['protein_var'] = res_df['genetic_variation_type'] == 'protein variant'

# Add 'intrinsic' flag
# A gene is intrinsic if it is NEITHER an acquired ARG NOR a protein variant
res_df['intrinsic'] = (res_df['acq_arg'] == False) & (res_df['protein_var'] == False)

# Value Lists (ALL)
unique_classes = get_unique_items(subset['drug_class'])
print(f"\nNumber of unique drug_class values: {len(unique_classes)}")
for item in unique_classes: 
    print(" •", item)

unique_agents = get_unique_items(subset['antimicrobial_agent'])
print(f"\nNumber of unique antimicrobial agent values: {len(unique_agents)}")
for item in unique_agents: 
    print(" •", item)

# Summary Tables (ALL)
drug_class_table = create_summary_table(subset, 'drug_class') # class, gene count, genes, contig
antimicrobial_table = create_summary_table(subset, 'antimicrobial_agent') 
# SAVE 
subset.to_csv('clean_AMR_data_subset.csv', index= False)
res_df.to_csv('Resistance_df.csv', index=False)
antimicrobial_table.to_csv('all_antimicrobials_table.csv', index = False)
drug_class_table.to_csv('all_drug_classes_table.csv', index = False)

# ---- SUB-DATAFRAMES: Acquired, protein vars and intrinsics
# Define the columns for the specific sub-DFs
display_cols = ['gene', 'contig', 'drug_class', 'antimicrobial_agent']

# Create the Acquired ARGs df
# Includes ResFinder hits + manual mobile additions (qacEdelta1, ble, etc.)
df_acquired = res_df[res_df['acq_arg'] == True][display_cols].copy()

# Create protein variants df
# Includes genes with specific resistance-conferring mutations
df_protein_variants = res_df[res_df['protein_var'] == True][display_cols].copy()

# Create "intrinsic" resistance
# Includes core chromosomal genes (e.g., acrAB, mdfA) not flagged above
df_intrinsic = res_df[res_df['intrinsic'] == True][display_cols].copy()

# print summary
print("\nBreakdown of unique AMR determinants:")
print(f" • Acquired ARGs: {len(df_acquired)}")
print(f" • Protein Variants: {len(df_protein_variants)}")
print(f" • Intrinsic Genes: {len(df_intrinsic)}")

# Save subdfs to CSV 
df_acquired.to_csv('Table_Acquired_ARGs.csv', index=False)
df_protein_variants.to_csv('Table_Protein_Variants.csv', index=False)
df_intrinsic.to_csv('Table_Intrinsic_Genes.csv', index=False)

# Create a Master Supplemental table: gene, contig, class, agent, res type
# add a 'Resistance_Type' label to each sub-dataframe
df_acquired['Resistance_Type'] = 'Acquired ARG'
df_protein_variants['Resistance_Type'] = 'Protein Variant'
df_intrinsic['Resistance_Type'] = 'Intrinsic Resistance'

# merge them
master_supplemental_table = pd.concat([
    df_acquired, 
    df_protein_variants, 
    df_intrinsic
], axis=0).reset_index(drop=True)

# clean strings (alphabetical sort)
master_supplemental_table['drug_class'] = master_supplemental_table['drug_class'].apply(
    lambda x: '; '.join(sorted(set(str(x).split(';')))) if pd.notna(x) else x
)
# clean spaces 
for col in ['gene', 'contig', 'drug_class', 'antimicrobial_agent']:
    master_supplemental_table[col] = master_supplemental_table[col].astype(str).str.strip()

# Replace 'nan' strings with NaN
master_supplemental_table = master_supplemental_table.replace('nan', pd.NA)

# save 
master_supplemental_table.to_csv('Master_AMR_Determinants_Supplemental.csv', index=False)
print(f"\nMaster Supplemental Table generated with {len(master_supplemental_table)} total entries.")

# ==========================================
# --- 6. DATA INTEGRITY & CATEGORY CHECK ---
# ==========================================

print("\n" + "="*40)
print("FINAL DATA INTEGRITY REPORT")
print("="*40)

# 1. Check for Overlapping Categories...
# A gene should only belong to ONE category (Acquired, Protein Var, or Intrinsic)
overlap = res_df[
    (res_df['acq_arg'] & res_df['protein_var']) | 
    (res_df['acq_arg'] & res_df['intrinsic']) | 
    (res_df['protein_var'] & res_df['intrinsic'])
]

if not overlap.empty:
    print(f"WARNING: {len(overlap)} entries have overlapping categories!")
    print(overlap[['gene', 'contig', 'acq_arg', 'protein_var', 'intrinsic']])
else:
    print("All genes have mutually exclusive categories.")

# 2. Check for Missing Categories...
# Every entry in res_df must be assigned to at least one category
unassigned = res_df[~(res_df['acq_arg'] | res_df['protein_var'] | res_df['intrinsic'])]

if not unassigned.empty:
    print(f"WARNING: {len(unassigned)} entries are unassigned!")
    print(unassigned[['gene', 'contig']])
else:
    print("No unassigned genes found.")

# 3. Contig-Specific Summary (verification for Heatmap)
print("\nFinal Determinant Count by Contig and Type:")
summary_pivot = res_df.groupby('contig')[['acq_arg', 'protein_var', 'intrinsic']].sum()
print(summary_pivot)

# 4. Check for 'Kpn' or 'qacE' remaining in the Intrinsic list
# This confirms that the string-matching fix worked
intrinsic_check = res_df[res_df['intrinsic'] == True]['gene'].unique()
suspect_mobile = [g for g in intrinsic_check if any(x in g.lower() for x in ['qace', 'ble', 'mrx', 'kpn'])]

if suspect_mobile:
    print(f"WARNING: Potential mobile genes found in intrinsic list: {suspect_mobile}")
else:
    print("intrinsic list looks OK.")


##############################################################
# --- 7. Plotting: Heatmaps for resistance to drug class by contig and res type
##############################################################

# Map the 'intrinsic' flag from res_df to 'subset' 

intrinsic_genes = res_df.loc[res_df['intrinsic'] == True, 'gene'].unique()
subset['is_intrinsic'] = subset['gene'].isin(intrinsic_genes)

# Generate presence/absence matrices for each res type
matrix_intrinsic = create_presence_matrix(subset, 'is_intrinsic')
matrix_acquired  = create_presence_matrix(subset, 'acq_arg')
matrix_prot      = create_presence_matrix(subset, 'is_protein_variant')

# get a master list of all drug classes found across all three categories 
# (excluding 'nan')
all_valid_classes = sorted(set(matrix_intrinsic.columns) | 
                           set(matrix_acquired.columns) | 
                           set(matrix_prot.columns))

# reindex all so they share the same X-axis
matrix_intrinsic = matrix_intrinsic.reindex(columns=all_valid_classes, fill_value=0)
matrix_acquired  = matrix_acquired.reindex(columns=all_valid_classes, fill_value=0)
matrix_prot      = matrix_prot.reindex(columns=all_valid_classes, fill_value=0)

# Test plot (optional):
fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
# Heatmap A: Intrinsic Resistance
sns.heatmap(matrix_intrinsic, ax=axes[0], cmap="Blues", 
            cbar=False, linewidths=.5, linecolor='lightgrey')
axes[0].set_title("Core Chromosomal Determinants", fontsize=14, pad=10)
axes[0].set_ylabel("")
axes[0].set_xlabel("")

# Heatmap B: Acquired ARGs (Mobilome)
sns.heatmap(matrix_acquired, ax=axes[1], cmap="Oranges", 
            cbar=False, linewidths=.5, linecolor='lightgrey')
axes[1].set_title("Acquired Resistance Genes", 
                  fontsize=14, pad=10)
axes[1].set_ylabel("")
axes[1].set_xlabel("")

# Heatmap C: Protein Variants (Evolutionary Adaptations)
sns.heatmap(matrix_prot, ax=axes[2], cmap="Greens", 
            cbar=False, linewidths=0.5, linecolor='lightgrey')
axes[2].set_title("Protein Variants", fontsize=14, pad=10)
axes[2].set_ylabel("")
axes[2].set_xlabel("Antimicrobial Drug Class")
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.tight_layout()
plt.savefig("AMR_Functional_Heatmap_Comparison_TEST.svg", bbox_inches='tight')
plt.show()

#### REFINED PLOT:

# Custom colormaps based on hex codes (Okabe ito)
cmap_intrinsic = LinearSegmentedColormap.from_list("int", ["#f8f9fa", "#0072B2"])
cmap_acquired = LinearSegmentedColormap.from_list("acq", ["#f8f9fa", "#E69F00"])
cmap_protein = LinearSegmentedColormap.from_list("pro", ["#f8f9fa", "#009E73"])

fig, axes = plt.subplots(
    3, 1, 
    figsize=(15, 10), 
    sharex=True, 
    gridspec_kw={'height_ratios': [1, 2, 2]} 
)

matrices = [matrix_intrinsic, matrix_acquired, matrix_prot]
custom_cmaps = [cmap_intrinsic, cmap_acquired, cmap_protein]
text_x_pos = matrices[0].shape[1] + 0.5

for i, (ax, matrix, cmap) in enumerate(zip(axes, matrices, custom_cmaps)):
    sns.heatmap(matrix, ax=ax, cmap=cmap, cbar=False, linewidths=0.5, linecolor='lightgrey')
    
    # Remove labels
    ax.set_ylabel("")
    ax.set_xlabel("")
    
    # header for the total n classes column
    if i == 0:
        ax.text(text_x_pos, -0.15, "N Classes", va='bottom', ha='left', 
                fontsize=11, fontweight='bold')

    # counts per row
    row_totals = matrix.sum(axis=1)
    for row_idx, total in enumerate(row_totals):
        ax.text(text_x_pos, row_idx + 0.5, 
                f"      {int(total)}", va='center', ha='left', 
                fontsize=10, fontweight='bold')

    # formatting: make sure top 2 plots don't have x labels
    if i < 2:
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=10)

# --- CREATE THE LEGEND ---
# Define the legend handles using the chosen colors
legend_handles = [
    mpatches.Patch(color='#0072B2', label='Intrinsic'),
    mpatches.Patch(color='#E69F00', label='Acquired'),
    mpatches.Patch(color='#009E73', label='Protein Variant')
]

fig.legend(
    handles=legend_handles, 
    title='Resistance Type\n', 
    title_fontsize='11', 
    fontsize='10', 
    loc='upper left', 
    bbox_to_anchor=(-0.025, 0.65), # edit to move the legend
    frameon=False
)

# adjust layout to accommodate legend and total col
plt.tight_layout()
plt.subplots_adjust(hspace=0.02, right=0.95, left=0.15) 

plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel("Antimicrobial Drug Class", fontsize=12, fontweight='bold', labelpad=15)
# save it!
plt.savefig("AMR_Heatmap_With_Legend.svg", bbox_inches='tight')
plt.show()

