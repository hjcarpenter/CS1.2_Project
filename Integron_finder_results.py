# -*- coding: utf-8 -*-
"""
This script cleans the output of Integron Finder

Created on Fri Nov 21 09:40:17 2025

@author: Heath
"""
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

# -------- Load result file (.integrons)

fpI = 'Results_Integron_Finder_final_assembly/final_assembly.integrons'
df = pd.read_csv(
    fpI,
    sep="\t",
    comment="#",
    na_values=["NA"],             # treat NA as missing
    keep_default_na=True
)

# ------ clean data

# lowercase all column names
df.columns = df.columns.str.lower()
# strip out whitespace from col headers
df.columns = df.columns.str.strip()
# add a source column 
df["source"] = "integronfinder"
# rename columns for downstream work
rename_map = {
    "id_replicon": "contig", 
    "pos_beg": "start",
    "pos_end":"end", 
    "distance_2attc": "distance_to_attc", 
    "type_elt": "type",
    "type": "element_type"
}
df = df.rename(columns=rename_map)

# fix annotation: carry over to new column and strip
df['name'] = df['annotation'].copy()
df['name'] = (
    df['name']
    .str.replace(r'-NCBIFAM$', '', regex=True)   # remove trailing '-NCBIFAM'
    .str.replace(r'^(SMR_|trim_)', '', regex=True)  # remove 'SMR_' or 'trim_' at start
)

# replace complete with a more informative string
df['element_type'] = df['element_type'].replace('complete', 'complete integron')

# add in integron boundary cols
df['integron_start'] = df.groupby('id_integron')['start'].transform('min')
df['integron_end']   = df.groupby('id_integron')['end'].transform('max')

# save cleaned results
df.to_csv('Integron_Finder_results.csv', index = False)

