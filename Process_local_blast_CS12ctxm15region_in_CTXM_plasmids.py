# -*- coding: utf-8 -*-
"""
Process the local blastn result for CS1.2 ISEcp1-CTXM15 insert (query) versus
the CTXM15-containing plasmid1 blast hits. 
Merge result with metadata 
Flag presence/absence of NDM5
Filter out sequences >=4000000 length (chromosomes) SAVE
Apply a strong filter for best hits with a date (100% identity) SAVE

Created on Sat Jan 17 11:34:36 2026

@author: Heath
"""

import pandas as pd

####### Imports #########

# NDM-5 containing plasmid metadata
fpndm = "D:/CRE_PROJECT_DATA_ANALYSIS/plasmid1_BLASTn/NDM-5_metadata_table.csv"
ndm = pd.read_csv(fpndm)

# Results from local blastn (ctxm15)
cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

blast = pd.read_csv(
    "query_vs_CTXM_plasmids.tsv",
    sep="\t",
    names=cols
)

# CTXM15 plasmids metadata
meta = pd.read_csv(
    "CTX-M-15_metadata_table.csv") # ctxm metadata

# 1: Local blast hittable merge to metadata
merge = blast.merge(
    meta,
    left_on="sseqid",
    right_on="accession",
    how="left"
)

# 2. flag accessions containing NDM
ndm['has_ndm'] = True 
filt = ndm[['accession', 'has_ndm']] # filter for accession and flag

# LEFT merge filt onto merge on accession
ctxplus = pd.merge(merge, filt, how = 'left', on = 'accession')
# make flag bool:
ctxplus["has_ndm"] = ctxplus["has_ndm"].fillna(False)

#  3. filter out obviously chromosomal hits
initial_rows = len(ctxplus)
# filter
ctxplus = ctxplus[ctxplus['length_y'] < 4000000]
# post-filter count
final_rows = len(ctxplus)
# print summary
rows_dropped = initial_rows - final_rows
print(f"Rows dropped: {rows_dropped}")
print(f"Rows remaining: {final_rows} (out of original {initial_rows})")

ctxplus.to_csv('CTXMplasmids_blast_hit_table_w_metadata.csv', index= False)

# ---- Now further filter ctxplus: 
# with hits pident 100, length >2000 
# these will be hits containing most of the isecp1-ctxm-orf477 unit seq. 
# NOTE: some of them contain a longer and a shorter hit,these will be dropped 
# but the aim is to get good quality hits with a year. 

# ensure numeric types 
ctxplus["pident"] = pd.to_numeric(ctxplus["pident"], errors="coerce")
ctxplus["length_x"] = pd.to_numeric(ctxplus["length_x"], errors="coerce")
ctxplus["collection_year"] = pd.to_numeric(ctxplus["collection_year"], errors="coerce")

ctxplus_isecp1 = ctxplus[
    (ctxplus["collection_year"].notna()) &
    (ctxplus["pident"] == 100) &
    (ctxplus["length_x"] > 2000)
].copy()

# save it
ctxplus_isecp1.to_csv('ctx_plasmids_strong_filter.csv', index=False)

# how many unique accessions?
n_unique = ctxplus_isecp1["accession"].nunique()
print(n_unique)

