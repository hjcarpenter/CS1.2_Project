# -*- coding: utf-8 -*-
"""
Import TADB finder results: 
    Filter by H value (blastp-based measurement of similarity to seqs in the TADB - higher = better)
    clean coords, name and flag for anno file. 

Created on Fri Jan  2 03:56:49 2026

@author: Heath
"""

import pandas as pd

fp = 'result_output.csv'
df = pd.read_csv(fp)

print(df.info())

# normalise col names: lc/ '_'
df.columns = (
    df.columns
      .str.strip()          # remove leading/trailing whitespace
      .str.lower()          # lowercase
      .str.replace(' ', '_', regex=False)  # spaces to underscore
)
# cols needed: 
    # start, end from coordinates
    # type 
    # + flag col --> 'TA system'

# extract start and end from coords
coords = df['coordinates'].astype(str).str.extract(r'(\d+)\.\.(\d+)')

df['start'] = pd.to_numeric(coords[0], errors='coerce')
df['end']   = pd.to_numeric(coords[1], errors='coerce')

# filter by Ha
threshold = 0.45
# coerce ha value to float
df['ha_value'] = pd.to_numeric(df['ha_value'], errors='coerce')
# Find pairs where the minimum ha_value in the pair meets the cutoff
# make sure: both members present, both numeric, both ha>=0.45 
valid_ta_ids = (
    df.groupby('ta_id')['ha_value']
      .filter(lambda x: x.notna().all() and (x >= threshold).all())
      .index
)

df_filtered = df.loc[valid_ta_ids]

# extract name from blast hit
df_filtered['name'] = (
    df_filtered['blast_hit']
    .astype(str)
    .str.extract(r'\(([^)]+)\)')
)

df_filtered[['blast_hit', 'name']].head()

# lowercase type col
df_filtered['type'] = df_filtered['type'].astype(str).str.lower()

# add a flag column
df_filtered['flag'] = 'TA system'

# save out all 
df_filtered.to_csv('plasmid1_TA_systems_all_cols_filtered.csv', index=False)

# now filter for anno cols
print(df_filtered.info())

cols_to_keep = ['ta_id', 'type', 'strand', 'ta_type', 'start', 'end', 'name','flag']
df_anno = df_filtered[cols_to_keep]

# save
df_anno.to_csv('plasmid1_TA_systems_anno.csv', index=False)
