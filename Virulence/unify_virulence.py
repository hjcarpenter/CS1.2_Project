# -*- coding: utf-8 -*-
"""
Virulence Genes: 
unify the results from Virulence Finder, VFDB (hits extracted from Bakta output) and 
ARMFinder Plus Virulence subset. 

The script takes the outputs from THE 3 tools, standardises their coordinates 
and names, outer-joins them by genomic region, builds unified gene, 
product, identity, and tool-support fields, merges overlapping calls for the 
same gene, and outputs a full unified table.

Created on Tue Nov 18 05:49:04 2025

@author: Heath
"""

import pandas as pd
import numpy as np

# Import tool output tsvs
fp1 = 'amrfinder_plus_Ecoli.tsv'
AFP = pd.read_csv(fp1, sep='\t')

fp2 = 'VirulenceFinder/results_tab.tsv'
VF = pd.read_csv(fp2, sep='\t')

# need to provide col headers for results extracted from bakta
fp3 = 'VFDB_hits_bakta/VFDB_hits.tsv'
cols = ['seq_id','type','start','end','strand','locus_tag','gene','product','db_xrefs']
BAK = pd.read_csv(fp3, sep='\t', skiprows=1, names=cols) # skip '#comment line'


# ------------ Prepare Data ------------

# Subset AFP for VIRULENCE hits
AFP = AFP[AFP['Type'].str.upper() == 'VIRULENCE']

# Make sure all dfs share contig, start, end cols 
AFP = AFP.rename(columns={
    "Contig id": "contig",
    "Start": "start",
    "Stop": "end"
})
BAK = BAK.rename(columns={"seq_id": "contig"})
VF = VF.rename(columns={
    "Contig": "contig"
})

# Extract the numeric positions for start and end from VF position in contig 
# extract two groups of digits around '..'
VF[['start', 'end']] = VF['Position in contig'].str.extract(r'(\d+)\.\.(\d+)')
# Convert to nullable integer type (allows NaN)
VF['start'] = pd.to_numeric(VF['start'], errors='coerce').astype('Int64')
VF['end']   = pd.to_numeric(VF['end'],   errors='coerce').astype('Int64')

# make sure all coordinates are integers
AFP['start'] = AFP['start'].astype(int)
AFP['end'] = AFP['end'].astype(int)
VF['start'] = VF['start'].astype(int)
VF['end'] = VF['end'].astype(int)

# clean contig names in VF
# Keep only the first token before the first space: "chromosome", "plasmid1", etc.
VF['contig'] = VF['contig'].str.split(' ', n=1).str[0]

# clean up the cols
# Col names:
def clean_columns(df):
    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(r'[^0-9a-zA-Z]+', '_', regex=True)  # replace ANY non-alphanumeric with underscore
        .str.replace(r'_+', '_', regex=True)             # collapse repeated underscores
        .str.strip('_')                                   # remove leading/trailing underscores
    )
    return df
AFP = clean_columns(AFP)
VF  = clean_columns(VF)
BAK = clean_columns(BAK)

# col content:
def clean_string_values(df):
    obj_cols = df.select_dtypes(include='object').columns
    df[obj_cols] = df[obj_cols].apply(lambda col: col.str.strip())
    return df

VF  = clean_string_values(VF)
AFP = clean_string_values(AFP)
BAK = clean_string_values(BAK)

# add a coverage column to VF
# Split "157 / 213" into two columns
nums = VF['query_template_length'].str.split('/', expand=True)

# Convert to numeric and compute coverage
VF['coverage'] = (
    nums[0].astype(float).div(nums[1].astype(float)) * 100
)

# drop rows in VF with < 90% coverage
VF = VF[VF['coverage'] >= 90].copy()

# add prefixes to all cols for reference
AFP_tagged = AFP.add_prefix('AFP_')
VF_tagged  = VF.add_prefix('VF_')
BAK_tagged = BAK.add_prefix('BAK_')

# rename the join cols to drop prefixes
AFP_tagged = AFP_tagged.rename(columns={
    'AFP_contig': 'contig',
    'AFP_start':  'start',
    'AFP_end':    'end'
})
VF_tagged = VF_tagged.rename(columns={
    'VF_contig': 'contig',
    'VF_start':  'start',
    'VF_end':    'end'
})
BAK_tagged = BAK_tagged.rename(columns={
    'BAK_contig': 'contig',
    'BAK_start':  'start',
    'BAK_end':    'end'
})

# ------------ Outer merge the dfs ------------------

# outer join BAK and AFP
BAK_AFP = BAK_tagged.merge(
    AFP_tagged,
    on=['contig', 'start', 'end'],
    how='outer'
)
# outer join Full with VF
FULL = BAK_AFP.merge(
    VF_tagged,
    on=['contig', 'start', 'end'],
    how='outer'
)
#print(FULL.info())
# infer gene name when all gene cols are nan
FULL['gene_inferred'] = (
    FULL['BAK_gene']
    .fillna(FULL['AFP_element_symbol'])
    .fillna(FULL['VF_virulence_factor'])
)

# if still empty, infer from product-related columns
mask_no_gene = FULL['gene_inferred'].isna()

FULL.loc[mask_no_gene, 'gene_inferred'] = (
    FULL['BAK_product']
    .where(mask_no_gene)
    .fillna(FULL['AFP_element_name'])
    .fillna(FULL['VF_protein_function'])
)
# now write value back to the correct column depending on which product field exists
FULL.loc[mask_no_gene & FULL['BAK_product'].notna(), 'BAK_gene'] = FULL.loc[mask_no_gene, 'BAK_product']
FULL.loc[mask_no_gene & FULL['AFP_element_name'].notna(), 'AFP_element_symbol'] = FULL.loc[mask_no_gene, 'AFP_element_name']
FULL.loc[mask_no_gene & FULL['VF_protein_function'].notna(), 'VF_virulence_factor'] = FULL.loc[mask_no_gene, 'VF_protein_function']

# ---- Add tool column based on which tools hit 
def get_tool(row):
    tools = []
    if pd.notna(row['BAK_gene']):
        tools.append('VirulenceFinderDatabase')
    if pd.notna(row['AFP_element_symbol']):
        tools.append('AMRFinderPlus')
    if pd.notna(row['VF_virulence_factor']):
        tools.append('VirulenceFinder')
    return '; '.join(tools) if tools else np.nan

FULL['tool'] = FULL.apply(get_tool, axis=1)

# Rearrange cols - put keys first for readability
cols = ['contig', 'start', 'end', 'tool',
        'BAK_gene', 'AFP_element_symbol', 'VF_virulence_factor', 'BAK_product',
        'AFP_element_name', 'VF_protein_function', 'BAK_db_xrefs', 
        'AFP_closest_reference_accession', 'VF_accession_number'
        ] + [
            c for c in FULL.columns if c not in ['contig', 'start', 'end', 
                                                 'tool', 'BAK_gene', 
                                                 'AFP_element_symbol', 
                                                 'VF_virulence_factor',
                                                 'BAK_product',
                                                 'AFP_element_name', 
                                                 'VF_protein_function', 
                                                 'BAK_db_xrefs', 
                                                 'AFP_closest_reference_accession', 
                                                 'VF_accession_number']]
FULL = FULL[cols]

# ---------- Unify genes -----------------

# coalesce the first non-null value across the relevant cols
FULL['gene'] = FULL[['BAK_gene',
                     'AFP_element_symbol',
                     'VF_virulence_factor']].bfill(axis=1).iloc[:, 0]
cols = FULL.columns.tolist()
cols.insert(3, cols.pop(cols.index('gene')))   # move 'gene' to position 3
FULL = FULL[cols]

# Drop all NAN cols
FULL = FULL.dropna(axis=1, how='all')

# Reorder cols
first_cols = [
    'contig', 'start', 'end', 'gene',
    'BAK_product', 'AFP_element_name', 'VF_protein_function',
    'AFP_identity_to_reference', 'VF_identity'
]
remaining_cols = [c for c in FULL.columns if c not in first_cols]
new_order = first_cols + remaining_cols
FULL = FULL[new_order]

# add an overlap flag
def find_overlaps(df):
    df = df.sort_values(['contig', 'start', 'end']).reset_index(drop=True)
    df['overlap'] = False
    
    for contig, group in df.groupby('contig'):
        idx = group.index
        starts = group['start'].values
        ends   = group['end'].values
        
        # Compare each interval to the next one
        for i in range(len(group) - 1):
            # If this end overlaps with the next start
            if ends[i] >= starts[i+1]:
                df.loc[idx[i],   'overlap'] = True
                df.loc[idx[i+1], 'overlap'] = True
    
    return df

FULL = find_overlaps(FULL)

# Merge overlapping rows that are the same (contig, gene)
def first_non_null(series):
    s = series.dropna()
    if not s.empty:
        return s.iloc[0]
    return np.nan

def merge_same_gene_overlaps(df):
    # Split into rows with and without gene name
    with_gene = df[df['gene'].notna()].copy()
    no_gene   = df[df['gene'].isna()].copy()
    
    merged_rows = []

    # Work per (contig, gene)
    for (contig, gene), grp in with_gene.groupby(['contig', 'gene']):
        grp = grp.sort_values(['start', 'end'])
        grp = grp.reset_index(drop=True)

        # clusters of overlapping intervals
        clusters = []
        current = [0]  # indices within grp

        for i in range(1, len(grp)):
            prev = grp.loc[current[-1]]
            cur  = grp.loc[i]

            # overlap condition: same contig+gene already guaranteed
            if cur['start'] <= prev['end']:
                # same gene, overlapping -> same cluster
                current.append(i)
            else:
                clusters.append(current)
                current = [i]

        clusters.append(current)  # last cluster

        # aggregate each cluster into a single row
        for cluster in clusters:
            sub = grp.loc[cluster]

            agg = {}
            for col in df.columns:
                if col == 'start':
                    agg['start'] = sub['start'].min()
                elif col == 'end':
                    agg['end'] = sub['end'].max()
                else:
                    agg[col] = first_non_null(sub[col])
            merged_rows.append(agg)

    merged_with_gene = pd.DataFrame(merged_rows, columns=df.columns)

    # Put everything back together
    out = pd.concat([merged_with_gene, no_gene], ignore_index=True)
    # sort nicely
    out = out.sort_values(['contig', 'start', 'end']).reset_index(drop=True)
    return out

FULL_merged = merge_same_gene_overlaps(FULL)

# add %identity col
FULL_merged['identity'] = (
    FULL_merged[['AFP_identity_to_reference', 'VF_identity']]
    .bfill(axis=1)
    .iloc[:, 0]
)

# -------- unify product
# create unified product column using first non-null across the three product/name fields
FULL_merged['product'] = FULL_merged[['BAK_product',
                                      'AFP_element_name',
                                      'VF_protein_function']].bfill(axis=1).iloc[:,0]

#move product column to sit next to gene for readability
cols = FULL_merged.columns.tolist()
cols.insert(cols.index('gene')+1, cols.pop(cols.index('product')))
FULL_merged = FULL_merged[cols]

# Specific cleaning 
# AslA -->
mask = FULL_merged['gene'] == 'AslA'
FULL_merged.loc[mask, 'gene'] = 'aslA'
FULL_merged.loc[mask, 'product'] = 'arylsulfatase AslA'

# Add a flag column
FULL_merged['flag'] = 'virulence'

# save this dataframe
FULL_merged.to_csv('VF_VFDB_AFP_virulence_full.csv', index = False)





