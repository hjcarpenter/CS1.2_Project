# -*- coding: utf-8 -*-
"""
Heatmap of AMR genotype, predicted phenotype, actual phenotype:
Script compares the experimentally observed AST data to the resfinder predicted phenotype and
the full genotype.
Generates a resistance heatmap with the associated genes 

Created on Sat Nov 29 04:59:57 2025
@author: Heath
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches


# ==========================================
# --- 1. HELPER FUNCTIONS ---
# ==========================================

def get_unique_items(series):
    """Splits, explodes, and cleans a column to find unique values."""
    # Define case-insensitive junk values
    junk = {'<na>', 'nan', 'none', '', 'null'}
    
    unique_array = (
        series.astype(str)
        .dropna()
        .str.split(';')
        .explode()
        .str.strip()
    )
    
    # Filter out anything that is the junk list (converted to lowercase for checking)
    clean_items = [item for item in unique_array if item.lower() not in junk]
    
    # Return as a deduplicated, sorted list
    return sorted(list(set(clean_items)))

def sort_gene_string(gene_str):
    """Splits a semicolon string, sorts alphabetically ignoring symbols, and rejoins."""
    if pd.isna(gene_str) or str(gene_str).strip() in ('', '<NA>', 'nan'):
        return ''
    
    # Split the string and clean up extra spaces
    genes_list = [g.strip() for g in str(gene_str).split(';') if g.strip()]
    
    # Custom sorting key: ignore case, delta, asterisks, and superscripts
    def clean_sort_key(g):
        return g.replace('∆', '').replace('*', '').replace('⁻', '').lower()
        
    # Deduplicate with set(), then sort using the custom key
    sorted_genes = sorted(list(set(genes_list)), key=clean_sort_key)
    
    return '; '.join(sorted_genes)

# ==========================================
# --- 2. IMPORTS ---
# ==========================================
fp1 = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/staramr_detailed_summary.tabular'
star = pd.read_csv(fp1, sep='\t')

fp2 = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/resistance_df.csv' 
res = pd.read_csv(fp2)

fp3 = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/CS1_2_AST_phenotype.csv' 
pheno = pd.read_csv(fp3, keep_default_na=False, na_values=['', 'NaN', 'null', 'N/A', '<NA>'])

fp4 = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/all_antimicrobials_table.csv' 
agent_gene = pd.read_csv(fp4)

fp5 = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/all_drug_classes_table.csv' 
dc_tbl = pd.read_csv(fp5)


# ==========================================
# --- 3. STARAMR DATA PREP ---
# ==========================================
star.columns = star.columns.str.lower().str.replace(' ', '_')

# cols to keep
keep = ['data', 'data_type', 'cge_predicted_phenotype']
star = star[keep].copy()

# drop rows where type is not Resistance
star = star[star['data_type'] == 'Resistance'].reset_index(drop=True)

# Rename 'data' --> gene
star = star.rename(columns={'data' : 'gene'})

# normalize `cge_predicted_phenotype` 
star['cge_predicted_phenotype'] = (
    star['cge_predicted_phenotype']
    .str.lower()
    .str.replace(',', ';', regex=False)
)

# normalize the gene column using exact match dict replacement
star['gene'] = star['gene'].replace({
    'parC (S80I)': 'parC',
    'parE (S458A)': 'parE',
    "aac(6')-Ib-cr": "aac(6')-Ib-cr5"  
})
star = star.drop_duplicates()


# ==========================================
# --- 4. RES DF & AGENT PREP ---
# ==========================================
res['gene'] = res['gene'].replace({
    'Escherichia coli parC conferring resistance to fluoroquinolones': 'parC'
})

agent_gene = agent_gene.rename(columns={'antimicrobial_agent': 'antibiotic'})


# ==========================================
# --- 5. MERGING & PHENOTYPE MAPPING ---
# ==========================================

# outer merge res onto star on gene 
stargenes = star.merge(res, on='gene', how='outer')

# Get unified set of unique agents
unique_agents = sorted(set(
    get_unique_items(stargenes['cge_predicted_phenotype']) +
    get_unique_items(res['antimicrobial_agent']) +
    get_unique_items(pheno['antibiotic'])
))

# Convert to sets for faster 'in' lookups
unique_cge = set(get_unique_items(stargenes['cge_predicted_phenotype']))
unique_hAMR = set(get_unique_items(stargenes['antimicrobial_agent']))

# create df_agents for merge with pheno
df_agents = pd.DataFrame({
    'antibiotic': unique_agents,
    'cge_pheno': [2 if a in unique_cge else 0 for a in unique_agents],
    'hAMR_pheno': [2 if a in unique_hAMR else 0 for a in unique_agents]
})

# merge df_agents onto pheno on antibiotic 
pheno_merge = pheno.merge(df_agents, on='antibiotic', how='left')

# Map the AST phenotypes
mapping = {'R': 2, 'S': 0, 'I': 1}
pheno_merge['phenotype'] = pheno_merge['phenotype'].map(mapping)

# bring in gene/ contig data by agent
pheno_merge = pd.merge(pheno_merge, agent_gene, on='antibiotic', how='left')

# ==========================================
# --- CLEAN GENE NAMES FOR PLOTTING ---
# ==========================================

# Define mapping dictionary (with Unicode superscript minus ⁻ for p variants)
gene_mapping = {
    'Escherichia coli AcrAB-TolC with AcrR mutation conferring resistance to ciprofloxacin, tetracycline, and ceftazidime': 'acrR⁻', 
    'Escherichia coli AcrAB-TolC with MarR mutations conferring resistance to ciprofloxacin and tetracycline': 'marR⁻',
    'Escherichia coli acrA': 'acrA',
    'Escherichia coli soxR with mutation conferring antibiotic resistance': 'soxR⁻',
    'Escherichia coli soxS with mutation conferring antibiotic resistance': 'soxS⁻',
    'Klebsiella pneumoniae KpnE': 'kpnE',
    'Klebsiella pneumoniae KpnF': 'kpnF', 
    'Escherichia coli gyrA conferring resistance to fluoroquinolones': 'gyrA⁻',
    'Escherichia coli parC conferring resistance to fluoroquinolones': 'parC⁻',
    'Escherichia coli mdfA': 'mdfA',
    'tetR': 'tetR⁻',
    'nfsA': 'nfsA⁻',
    'ftsI': 'ftsI⁻',
}

# Create new gene column
pheno_merge['gene_clean'] = pheno_merge['gene'].astype(str)

# Loop through the dictionary and replace substrings
for long_name, short_name in gene_mapping.items():
    pheno_merge['gene_clean'] = pheno_merge['gene_clean'].str.replace(long_name, short_name, regex=False)

# parE in case
pheno_merge['gene_clean'] = pheno_merge['gene_clean'].str.replace('parE (S458A)', 'parE⁻', regex=False)

# Clean up 'nan' strings
pheno_merge['gene_clean'] = pheno_merge['gene_clean'].replace('nan', pd.NA)

# ==========================================
# --- MANUAL OVERRIDES (FROM OUTSIDE ANALYSIS) ---
# ==========================================

# Update partial genes to include the Delta symbol (∆)
pheno_merge['gene_clean'] = (
    pheno_merge['gene_clean']
    .str.replace('catB3', '∆catB3', regex=False)
    # Catching both exact and lowercase just in case
    .str.replace('dfrA17', '∆dfrA17', regex=False) 
    .str.replace('dfra17', '∆dfrA17', regex=False) 
)

# Helper function to safely add new genes to the semicolon-separated strings
def append_new_genes(current_str, new_genes_list):
    if pd.isna(current_str) or str(current_str).strip() in ('', '<NA>', 'nan'):
        current_genes = []
    else:
        current_genes = [g.strip() for g in str(current_str).split(';')]
        
    # Combine, deduplicate (using set), and alphabetize
    updated_genes = sorted(list(set(current_genes + new_genes_list)))
    return '; '.join(updated_genes)


# Update Cefpodoxime / CPD
mask_cpd = pheno_merge['antibiotic'].isin(['cefpodoxime', 'CPD'])

# Force genotype (hAMR) to Resistant (2)
pheno_merge.loc[mask_cpd, 'hAMR_pheno'] = 2
# Add blaCTX-M-15 to the gene list
pheno_merge.loc[mask_cpd, 'gene_clean'] = pheno_merge.loc[mask_cpd, 'gene_clean'].apply(
    lambda x: append_new_genes(x, ['blaCTX-M-15*'])
)

# Update Trimethoprim+Sulfamethoxazole / SXT
mask_sxt = pheno_merge['antibiotic'].isin(['trimethoprim+sulfamethoxazole', 'SXT'])

# Force genotype (hAMR) to Resistant (2)
pheno_merge.loc[mask_sxt, 'hAMR_pheno'] = 2
# Add sul1, dfrA12, and ∆dfrA17 to the gene list
pheno_merge.loc[mask_sxt, 'gene_clean'] = pheno_merge.loc[mask_sxt, 'gene_clean'].apply(
    lambda x: append_new_genes(x, ['sul1*', 'dfrA12*', '∆dfrA17*'])
)

# ---------- Generate Resistance Heat map

# 1. SORT AND EXTRACT LABELS FIRST
pheno_sort = pheno_merge.sort_values('antibiotic').reset_index(drop=True)

# Using abbreviations
antibiotics = pheno_sort['abbreviation'].tolist() 

# Apply sorting function to every row in the gene_clean column
raw_genes = pheno_sort['gene_clean'].tolist()
genes = [sort_gene_string(g) for g in raw_genes]

# Wrap the long gene strings. 
# width=N means it will drop to a new line after ~N characters.
# adjust the number up or down to make the text block wider or narrower
wrapped_genes = [textwrap.fill(g, width=50) for g in genes]

# 2. PREP DATA FOR HEATMAP
heat_cols = ['phenotype', 'cge_pheno', 'hAMR_pheno']
heat_df = pheno_sort.set_index('abbreviation')[heat_cols]

heat_df = heat_df.rename(
    columns={'phenotype': 'Phenotype',
             'cge_pheno': 'ResFinder 4.0',
             'hAMR_pheno': 'Genotype'}
)

data = heat_df.to_numpy(dtype=float)
data_masked = np.ma.masked_invalid(data)

# Setup figure 
fig, ax = plt.subplots(figsize=(9, 0.5 * len(heat_df) + 2))

# 3. CUSTOM COLORS & BOUNDARIES
color_list = ['mediumseagreen', 'orange', 'firebrick']
cmap = mcolors.ListedColormap(color_list)
bounds = [-0.5, 0.5, 1.5, 2.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

im = ax.imshow(data_masked, aspect='auto', cmap=cmap, norm=norm)

# 4. LEFT AXIS (Antibiotics)
ax.set_yticks(np.arange(len(heat_df)))
ax.set_yticklabels(antibiotics)
ax.set_xticks(np.arange(len(heat_df.columns)))
ax.set_xticklabels(heat_df.columns, rotation=45, ha='center')

# 5. RIGHT AXIS (Genes)
ax2 = ax.twinx()  
ax2.set_ylim(ax.get_ylim())  
ax2.set_yticks(np.arange(len(heat_df)))
# Pass the 'wrapped_genes' list here instead of 'genes'
ax2.set_yticklabels(wrapped_genes, ha='left', fontsize=8, fontstyle='italic')

ax2.tick_params(axis='y', length=0, pad=10) 
ax2.spines[:].set_visible(False)

# ADD 'Associated Genes' Header
# x=1.00 = right edge
# y=1.00 = top edge
ax.text(1.04, 1.01, 'Associated Genes', 
        transform=ax.transAxes, 
        ha='left', va='bottom', 
        fontweight='bold', fontsize=10)

# 6. LEGEND (Anchored below the main axis)
legend_patches = [
    mpatches.Patch(color='mediumseagreen', label='Sensitive'),
    mpatches.Patch(color='orange', label='Intermediate'),
    mpatches.Patch(color='firebrick', label='Resistant')
]

# Attach to the main axis (ax), place below, and spread out horizontally (ncol=3)
ax.legend(handles=legend_patches, 
          bbox_to_anchor=(0.5, -0.15), 
          loc='upper center', 
          ncol=3,
          borderaxespad=0.,
          frameon=False)

# Footnote (Bottom Right) <-- optional. Leave off if too busy/ ref in legend instead. 
# y=-0.02 places it just below the bottom edge of the heatmap
# ax.text(1.1, -0.03, '⁻ Protein variant\n* Inferred from literature', 
#         transform=ax.transAxes, 
#         ha='left', va='top', 
#         fontsize=9, color='dimgrey')

# 7. AESTHETICS
ax.set_xticks(np.arange(len(heat_df.columns)) - 0.5, minor=True)
ax.set_yticks(np.arange(len(heat_df)) - 0.5, minor=True)
ax.grid(which="minor", color="white", linestyle='-', linewidth=2)
ax.tick_params(which="minor", bottom=False, left=False)
ax.spines[:].set_visible(False)

plt.tight_layout()
# Note: tight_layout might sometimes clip bottom legends. 
# bbox_inches='tight' in savefig will rescue it
plt.savefig('pheno_geno_heatmap.svg', bbox_inches='tight') 
plt.show()