# -*- coding: utf-8 -*-
"""
This script Cleans and filters the mobile element finder results
1. Import data from tsv
2. Rename columns and normalise naming and values
3. Filter as:
evalue <= 1e-10]
identity>= 30
coverage>= 40
allele_len >= 300 # remove very short alignmnets
4. Drop duplicates keeping the best hit per contig+start+end based on coverage, identity, gaps, substitution, evalue.
5. save

Created on Tue Nov 25 06:34:16 2025
@author: Heath
"""
import pandas as pd
import matplotlib.pyplot as plt

fp = 'mgefinder_result.tsv'
df = pd.read_csv(fp, sep='\t')

# ----  clean data 

# rename columns for downstream work
rename_map = {
    "stop": "end",
    "e_value": "evalue",
    "prediction method": "prediction_method"
}
df = df.rename(columns=rename_map)

# Normalise the contig names
df['contig'] = df['contig'].str.split().str[0] 
# normalise identity and coverage: decimal -> %
df['identity']  = df['identity']  * 100
df['coverage']  = df['coverage']  * 100
# lc type column
df['type'] = df['type'].str.strip().str.lower()
# add a source column 
df["source"] = "mobileelementfinder"

print(df.info())

# --- Filtering
# evalue
mef_filt = df[df["evalue"] <= 1e-10].copy()
# identity and coverage
mef_filt = mef_filt[
    (mef_filt["identity"] >= 30) &
    (mef_filt["coverage"] >= 40)
]
# short alignmnets
mef_filt = mef_filt[mef_filt["allele_len"] >= 300]

# have a look at duplicates to decide on next step
dupe_mask = mef_filt.duplicated(subset=['contig','start','end'], keep=False)
dupes = mef_filt[dupe_mask].sort_values(['contig','start','end'])
non_dupes = mef_filt[~dupe_mask].copy()
dupes.groupby(['contig','start','end']).apply(lambda x: x)

dupes_sorted = dupes.sort_values(
    ['contig', 'start', 'end',
     'coverage', 'identity', 'gaps', 'substitution', 'evalue'],
    ascending=[True,   True,   True,
               False,   False,   True,  True,          True]
)

# keep best per contig+start+end, regardless of name
dupes_best = dupes_sorted.drop_duplicates(
    subset=['contig', 'start', 'end'],
    keep='first'
)

keep_mask     = ~dupes_sorted.duplicated(subset=['contig', 'start', 'end'],
                                         keep='first')
dupes_kept    = dupes_sorted[keep_mask].copy()
dupes_dropped = dupes_sorted[~keep_mask].copy()

df_final = (
    pd.concat([non_dupes, dupes_best], ignore_index=True)
      .sort_values(['contig', 'start', 'end'])
      .reset_index(drop=True)
)

# save
df_final.to_csv('MobileElementFinder_filtered.csv', index=False)



