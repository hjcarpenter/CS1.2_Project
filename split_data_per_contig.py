# -*- coding: utf-8 -*-
"""
This script splits the results from Bakta annotation, the cleaned combined AMR results and the cleaned combined virulence results by contig,
to create separate CSVs. These form the basis of the annotation CSVs. 

Created on Fri Dec  5 07:54:08 2025

@author: Heath
"""

import pandas as pd

df = pd.read_csv("bakta_CS1_2.tsv", sep="\t", skiprows=5) # skip comment lines

# clean the headers
df.columns = df.columns.str.replace("^#", "", regex=True)
df.columns = df.columns.str.replace(" ", "_")
df.columns = df.columns.str.lower()

# See available contigs
print(df["sequence_id"].unique())

# extract per contig and save out
for seq in df["sequence_id"].unique():
    df_sub = df[df["sequence_id"] == seq]
    df_sub.to_csv(f"bakta_CS1_2_{seq}.tsv", sep="\t", index=False)

# -------------
# split additional annotation files - extract chromosome and plasmid1

# import amr and virulence dfs
amrfp = 'D:/CRE_PROJECT_DATA_ANALYSIS/AMR/AMR_results_analysis/AMR_subset.csv'
virfp = 'D:/CRE_PROJECT_DATA_ANALYSIS/VIRULENCE/VF_VFDB_AFP_virulence_full.csv'

amr = pd.read_csv(amrfp)
vir = pd.read_csv(virfp)

# --- clean up the virulence df
keep_cols = ['contig', 'start', 'end', 'gene', 'product', 'BAK_strand', 'AFP_strand']
vir = vir[keep_cols]
# create a unified strand col
vir["strand"] = vir["BAK_strand"].combine_first(vir["AFP_strand"])
vir = vir.drop(columns = ["BAK_strand", "AFP_strand"])
# subset for plasmid1 only
vir_pl1 = vir[vir["contig"] == "plasmid1"].sort_values(by=["start", "end"])
# subset for chr only
vir_chr = vir[vir["contig"] == "chromosome"].sort_values(by=["start", "end"])
# save
vir_pl1.to_csv("plasmid1_virulence_genes.csv", index=False)
vir_chr.to_csv("chromosome_virulence_genes.csv", index=False)

# ---- amr df
print(amr.info())
cols = ['contig', 'start', 'end', 'gene', 'strand' ]
amr = amr[cols]
# subset for plasmid1
amr_pl1 = amr[amr["contig"] == "plasmid1"].sort_values(by=["start", "end"])
# subset for chromosome
amr_chr = amr[amr["contig"] == "chromosome"].sort_values(by=["start", "end"])
# save
amr_pl1.to_csv("plasmid1_amr_genes.csv", index=False)
amr_chr.to_csv("chromosome_amr_genes.csv", index=False)



