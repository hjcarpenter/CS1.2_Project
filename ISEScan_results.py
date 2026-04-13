# -*- coding: utf-8 -*-
"""
ISEScan results: 
This script cleans and filters the output from ISEScan 
1. clean tabular file
2. save raw data
3. filter the result:
- drop partial hits
- evalue <= 1e-10
- tnp_orf_len >=300 # filter out very small tnp orfs
4. save

Created on Tue Nov 25 09:38:26 2025

@author: Heath
"""

import pandas as pd

fp = 'tabular results.tabular'
df = pd.read_csv(fp, sep='\t')

# strip out whitespace from col headers
df.columns = df.columns.str.strip()

print(df.info())

# rename columns
df_clean = df.rename(columns={"seqID": "contig", 
                              "family": "is_family", 
                              "isBegin": "start", 
                              "isEnd": "end", 
                              "start1":"ir_start1", 
                              "end1": "ir_end1", 
                              "start2":"ir_start2", 
                              "end2": "ir_end2",
                              "orfBegin": "tnp_orf_start",
                              "orfEnd": "tnp_orf_end", 
                              "orfLen": "tnp_orf_len", 
                              "cluster": "tnp_cluster", 
                              "E-value": "evalue_best",
                              "E-value4copy": "evalue",
                              "type": "hit_type"})

# Add columns
df_clean['name'] = df_clean['is_family']
df_clean['source'] = 'isescan'
df_clean['type'] = 'insertion sequence'
print(df.info())
# make type col more intelligible
df_clean['hit_type'] = df_clean['hit_type'].map({'c': 'complete', 'p': 'partial'})

# save raw as csv
df_clean.to_csv('ISEScan_raw_results.csv', index = False)

# ---- Filtering
# drop partials
filt = df_clean[df_clean["hit_type"] == "complete"].copy()
# keep only high confidence  (HMMER confidence threshold)
filt = filt[filt["evalue"] <= 1e-10]
# filter out tiny tnp ORFs
filt = filt[filt["tnp_orf_len"] >= 300]   # at least 300 bp

# save
filt.to_csv('ISEScan_results_filtered.csv', index=False)