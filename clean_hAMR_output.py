# -*- coding: utf-8 -*-
"""
### Clean the hAMRonization output ####

Imports the hAMRonization TSV output, saves a raw copy as csv.
IMPORTANT: Also imports the Resfinder tabular output from starAMR to backfill predictions lost by hAMRonization!!
Cleans gene symbols (e.g. adding 'bla' to beta-lactamases, fixing specific symbols, and normalising RGI capitalisation). 
Normalises gene coordinates, groups overlapping hits on the same contig and gene-prefix into clusters, applying a special rule to avoid merging the KpnE and KpnF calls.
Within each locus (contig + gene-prefix + cluster), it merges results from
all tools using tool-specific priority rules for key annotation fields (e.g., 'drug class':RGI > AMRFinder > staramr) and concatenates all non-null values for the remaining columns to preserve data. 
Then adds summary information on tool support, renames some columns for downstream compatibility, and saves the final table as a csv. 

Created on Thu Nov 27 04:40:18 2025
@author: Heath
"""

import pandas as pd
import numpy as np

# ==========================================
# --- CONSTANTS & CONFIGURATION ---
# ==========================================

# to ensure everything is in CARD ontology format
CLASS_MAPPING = {
    'trimethoprim': 'diaminopyrimidine antibiotic',
    'quinolone': 'fluoroquinolone antibiotic',
    'phenicol': 'phenicol antibiotic',
    'beta-lactam': 'penicillin beta-lactam',
    'bleomycin': 'glycopeptide antibiotic', 
    'aminoglycoside': 'aminoglycoside antibiotic', 
    'macrolide': 'macrolide antibiotic',
    'sulfonamide': 'sulfonamide antibiotic',
    'tetracycline': 'tetracycline antibiotic'
}

# priority order/s to apply later
PRIORITY_RGI_FIRST = ['rgi', 'amrfinderplus', 'staramr']
PRIORITY_AMRF_FIRST = ['amrfinderplus', 'rgi', 'staramr']

LOCUS_COLS = ['input_sequence_id', 'gene_prefix', 'cluster_id']
RGI_PRIORITY_COLS = ['resistance_mechanism']
AMRF_PRIORITY_COLS = [
    'start_corrected', 'stop_corrected', 'input_gene_start', 
    'input_gene_stop', 'gene_symbol', 'genetic_variation_type', 'gene_name'
]
PRESERVE_ALL_COLS = ['antimicrobial_agent', 'drug_class']

SPECIAL_COLS = (
    set(LOCUS_COLS + ['analysis_software_name']) 
    | set(RGI_PRIORITY_COLS) 
    | set(AMRF_PRIORITY_COLS) 
    | set(PRESERVE_ALL_COLS)
)

# ==========================================
# --- HELPER FUNCTIONS ---
# ==========================================

def lower_first_letter(s):
    """Lowercase first letter only."""
    if pd.isna(s) or len(s) == 0:
        return s
    return s[0].lower() + s[1:]

def update_agents(row, reference_dict):
    """Check and update agents in df against RFdf reference."""
    gene = row['gene_symbol']
    current_val = row['antimicrobial_agent']
    
    if gene not in reference_dict:
        return current_val
        
    rf_agents = reference_dict[gene]
    
    if pd.isna(current_val) or current_val.strip() == '':
        current_list, current_set = [], set()
    else:
        current_list = current_val.split('; ')
        current_set = set(current_list)
        
    added_new = False
    for agent in rf_agents:
        if agent not in current_set and agent != '':
            current_list.append(agent)
            added_new = True
            
    return '; '.join(current_list) if added_new else current_val

def choose_by_priority(group, col, priority_order):
    if col not in group.columns:
        return pd.NA

    for tool in priority_order:
        subset = group.loc[group['analysis_software_name'] == tool, col].dropna()
        if not subset.empty:
            return subset.iloc[0]

    non_null = group[col].dropna()
    return non_null.iloc[0] if not non_null.empty else pd.NA

def concat_all_non_null(group, col):
    if col not in group.columns:
        return pd.NA

    vals = (
        group[col]
        .dropna()
        .astype(str)
        .str.strip()
        .loc[lambda s: s.ne('')]
        .drop_duplicates()
        .tolist()
    )
    return pd.NA if not vals else ';'.join(vals)

def concat_all_non_null_by_tool_priority(group, col, priority_order):
    """Concatenate all non-null values, deduped, ordered by tool priority."""
    if col not in group.columns:
        return pd.NA

    seen = set()
    out = []

    # 1) values in priority tool order
    for tool in priority_order:
        s = group.loc[group['analysis_software_name'] == tool, col].dropna().astype(str).str.strip()
        for v in s:
            if v and v not in seen:
                seen.add(v)
                out.append(v)

    # 2) remaining values
    s = group[col].dropna().astype(str).str.strip()
    for v in s:
        if v and v not in seen:
            seen.add(v)
            out.append(v)

    return pd.NA if not out else ';'.join(out)

def merge_locus(group):
    out = {}
    seq, prefix, cid = group.name
    out['input_sequence_id'] = seq
    out['gene_prefix'] = prefix
    out['cluster_id'] = cid

    out['analysis_software_name'] = concat_all_non_null(group, 'analysis_software_name')

    for col in PRESERVE_ALL_COLS:
        if col in group.columns:
            out[col] = concat_all_non_null_by_tool_priority(group, col, PRIORITY_RGI_FIRST)

    for col in RGI_PRIORITY_COLS:
        if col in group.columns:
            out[col] = choose_by_priority(group, col, PRIORITY_RGI_FIRST)

    for col in AMRF_PRIORITY_COLS:
        if col in group.columns:
            out[col] = choose_by_priority(group, col, PRIORITY_AMRF_FIRST)

    other_cols = [c for c in group.columns if c not in SPECIAL_COLS]
    for col in other_cols:
        out[col] = concat_all_non_null(group, col)

    return pd.Series(out)

def clean_and_dedupe(val, mapping_dict=None, remove_terms=None):
    """Splits, maps, removes exclusions, and deduplicates semicolon strings."""
    if pd.isna(val):
        return val

    ignore_list = {'<na>', 'nan', 'na', 'none', 'null'}
    if remove_terms:
        ignore_list.update([t.lower() for t in remove_terms])

    parts = [p.strip() for p in str(val).split(';') if p.strip()]
    seen, deduped = set(), []
    
    for p in parts:
        if p.lower() in ignore_list:
            continue
            
        if mapping_dict:
            p = mapping_dict.get(p, p)
            
        if p not in seen:
            seen.add(p)
            deduped.append(p)

    return np.nan if not deduped else '; '.join(sorted(deduped))

def add_ceph(row):
    """Safely append cephalosporin to drug_class."""
    current_class = str(row['drug_class']) if pd.notna(row['drug_class']) else ""
    if 'cephalosporin' not in current_class:
        if current_class in ("", "nan"):
            return 'cephalosporin'
        else:
            return current_class + '; cephalosporin'
    return current_class


# ==========================================
# --- MAIN SCRIPT ---
# ==========================================

# --- 1. Load Data ---
fp = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/hAMRonize/combined_AMR_report.tsv' # FILEPATH to output from hAMRonization
df = pd.read_csv(fp, sep='\t')
df.to_csv('hAMRonize_Raw_Results.csv', index=False)

RFdffp = "D:/CRE_PROJECT_DATA_ANALYSIS/AMR/staramr/resfinder.tabular" # FILEPATH starAMR resfinder output
RFdf = pd.read_csv(RFdffp, sep='\t')


# --- 2. Clean hAMR df ---
# Beta-lactamases: prepend "bla" 
mask_bla = (df['analysis_software_name'] == 'rgi') & (df['gene_name'].str.contains('beta-lactamase', case=False, na=False))
df.loc[mask_bla & ~df['gene_symbol'].str.startswith('bla', na=False), 'gene_symbol'] = 'bla' + df.loc[mask_bla, 'gene_symbol'].astype(str)

# specific-case fixes 
df.loc[df['gene_symbol'] == 'mphA', 'gene_symbol'] = 'mph(A)'
mask_aac = df['gene_symbol'].str.startswith('AAC', na=False)
df.loc[mask_aac, 'gene_symbol'] = 'aac' + df.loc[mask_aac, 'gene_symbol'].str[3:]

# rename: gyrA and parC to match RGI calls
df['gene_symbol'] = df['gene_symbol'].replace({
    'gyrA': 'Escherichia coli gyrA conferring resistance to fluoroquinolones',
    'parC': 'Escherichia coli parC conferring resistance to fluoroquinolones'
})

# Apply first-letter-lowercasing for RGI genes
mask_rgi = (df['analysis_software_name'] == 'rgi')
exclude = (
    df['gene_symbol'].isin(['CRP', 'H-NS']) |
    df['gene_symbol'].str.startswith('Escherichia', na=False) |
    df['gene_symbol'].str.startswith('Haemophilus', na=False) |
    df['gene_symbol'].str.startswith('Klebsiella', na=False)
)
mask_modify = mask_rgi & ~exclude
df.loc[mask_modify, 'gene_symbol'] = df.loc[mask_modify, 'gene_symbol'].astype(str).apply(lower_first_letter)

# StarAMR: Move drug_class to antimicrobial_agent
mask_star = df['analysis_software_name'] == 'staramr'
df.loc[mask_star, 'antimicrobial_agent'] = df.loc[mask_star, 'drug_class']
df.loc[mask_star, 'drug_class'] = pd.NA

# Clean up StarAMR entries
df.loc[mask_star, 'antimicrobial_agent'] = (
    df.loc[mask_star, 'antimicrobial_agent']
      .str.replace('I/R', '', regex=False)
      .str.replace('/', '+', regex=False)
      .str.strip()
)

# Lowercase and clean separators
df['antimicrobial_agent'] = df['antimicrobial_agent'].astype(str).str.lower().replace('nan', np.nan)
df['drug_class'] = df['drug_class'].astype(str).str.lower().replace('nan', np.nan)

df['antimicrobial_agent'] = df['antimicrobial_agent'].str.replace('/', '; ').str.replace(', ', '; ').str.replace(r'\s*;\s*', '; ', regex=True).str.strip()
df['drug_class'] = df['drug_class'].str.replace('/', '; ').str.replace(r'\s*;\s*', '; ', regex=True).str.strip()


# --- 3. Clean RFdf ---
RFdf = (
    RFdf
    .rename(columns=lambda x: x.lower().replace(' ', '_'))
    .rename(columns={'cge_predicted_phenotype': 'antimicrobial_agent', 'gene':'gene_symbol'})
    .assign(
        antimicrobial_agent=lambda df_: df_['antimicrobial_agent'].str.lower().str.replace(', ', '; '),
        gene_symbol=lambda df_: df_['gene_symbol'].replace({"aac(6')-Ib-cr": "aac(6')-Ib10"})
    )
)


# --- 4. Merge Stage (RFdf -> hAMR df) ---
rf_reference = (
    RFdf.dropna(subset=['antimicrobial_agent'])
    .groupby('gene_symbol')['antimicrobial_agent']
    .apply(lambda x: set('; '.join(x).split('; ')))
    .to_dict()
)

# Apply dictionary using lambda to pass reference
df['antimicrobial_agent'] = df.apply(lambda row: update_agents(row, rf_reference), axis=1)


# --- 5. Overlapping Hits Clustering ---
df['start_corrected'] = df[['input_gene_start', 'input_gene_stop']].min(axis=1)
df['stop_corrected']   = df[['input_gene_start', 'input_gene_stop']].max(axis=1)

df2 = df.copy()
df2['gene_prefix'] = df2['gene_symbol'].astype(str).str[:5]

TOL = 0 
clusters = []  
for (seq, prefix), sub in df2.groupby(['input_sequence_id', 'gene_prefix']):
    sub = sub.sort_values('start_corrected')
    cluster_id = 0
    current_cluster_end = None
    
    for idx, row in sub.iterrows():
        start, end = row['start_corrected'], row['stop_corrected']
        if current_cluster_end is None:
            cluster_id += 1
            current_cluster_end = end
        else:
            if start <= current_cluster_end + TOL:
                current_cluster_end = max(current_cluster_end, end)
            else:
                cluster_id += 1
                current_cluster_end = end
        
        clusters.append({
            'index': idx, 'input_sequence_id': seq,
            'gene_prefix': prefix, 'cluster_id': cluster_id
        })

cluster_df = pd.DataFrame(clusters).set_index('index')
df2['cluster_id'] = cluster_df['cluster_id']

# Special overlap case fix: prevent merging of KpnE and KpnF
gene_e, gene_f = "Klebsiella pneumoniae KpnE", "Klebsiella pneumoniae KpnF"
bad_pair = {gene_e, gene_f}
max_cluster = df2['cluster_id'].max()

for (seq, cid), sub in df2.groupby(['input_sequence_id', 'cluster_id']):
    if bad_pair.issubset(set(sub['gene_symbol'])):
        max_cluster += 1
        mask = (df2['input_sequence_id'] == seq) & (df2['cluster_id'] == cid) & (df2['gene_symbol'] == gene_f)
        df2.loc[mask, 'cluster_id'] = max_cluster

# Generate cluster summary
summary2 = (
    df2.groupby(LOCUS_COLS)
      .agg({
          'start_corrected': 'min',
          'stop_corrected': 'max',
          'analysis_software_name': lambda x: sorted(set(x)),
          'gene_symbol': lambda x: sorted(set(x))
      })
      .reset_index()
)

summary2['n_tools'] = summary2['analysis_software_name'].apply(len)
summary2['n_genes'] = summary2['gene_symbol'].apply(len)
summary2['gene_disagreement'] = summary2['n_genes'] > 1


# --- 6. Duplicate Handling & Priority Merging ---
df3 = df2.copy().dropna(axis=1, how='all')

df_locus_rep = (
    df3
    .groupby(LOCUS_COLS)
    .apply(merge_locus, include_groups=False)
    .reset_index(drop=True)
)


# --- 7. Final String Formatting & Filtering ---
# Clean basic drug classes
df_locus_rep['drug_class'] = df_locus_rep['drug_class'].apply(
    lambda x: clean_and_dedupe(x, mapping_dict=CLASS_MAPPING)
)

# Add cephalosporin where required
mask_ceph = df_locus_rep['antimicrobial_agent'].str.contains('cephalosporin', case=False, na=False)
df_locus_rep.loc[mask_ceph, 'drug_class'] = df_locus_rep.loc[mask_ceph].apply(add_ceph, axis=1)

# Clean exclusions from columns using the clean_and_dedupe function
df_locus_rep['antimicrobial_agent'] = df_locus_rep['antimicrobial_agent'].apply(
    lambda x: clean_and_dedupe(x, remove_terms=['cephalosporin', 'quinolone'])
)

df_locus_rep['drug_class'] = df_locus_rep['drug_class'].apply(
    lambda x: clean_and_dedupe(x, remove_terms=['efflux'])
)

# Filter rows
df_locus_rep = df_locus_rep[df_locus_rep['gene_prefix'] != 'Haemo']

# --- 8. Append Summary Info & Save ---
summary_info = summary2[LOCUS_COLS + ['n_tools', 'analysis_software_name']].copy()
summary_info = summary_info.rename(columns={'n_tools': 'n_tools_supporting', 'analysis_software_name': 'supporting_tools'})

df_locus_rep = df_locus_rep.merge(summary_info, on=LOCUS_COLS, how='left')
df_locus_rep = df_locus_rep.dropna(axis=1, how='all')

df_locus_rep = df_locus_rep.rename(columns={
    "start_corrected": "start", "stop_corrected": "end",
    "sequence_identity": "identity", 'strand_orientation': 'strand', 
    'gene_symbol': 'gene', 'input_sequence_id': 'contig'
})

df_locus_rep.to_csv("hAMR_Result_Clean.csv", index=False)