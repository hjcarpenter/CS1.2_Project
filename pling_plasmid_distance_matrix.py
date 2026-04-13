# -*- coding: utf-8 -*-
"""
This script creates a cluster heatmap from the Pling distance matrix
Applies UPGMA clustering for the tree

Created on Fri Mar 20 05:29:04 2026

@author: Heath
"""

import pandas as pd
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

# 1. Load and prepare square matrix 

fp = 'D:/CRE_PROJECT_DATA_ANALYSIS/pling/output_dir_cd_03/all_plasmids_distances.tsv' # <-- fp to distance tsv
df = pd.read_csv(fp, sep='\t')
df_reverse = df.rename(columns={'plasmid_1': 'plasmid_2', 'plasmid_2': 'plasmid_1'})
df_symmetric = pd.concat([df, df_reverse], ignore_index=True)
matrix = df_symmetric.pivot(index='plasmid_1', columns='plasmid_2', values='distance')

# Fill diagonals with 0 (distance to self) and clean up missing connections
for col in matrix.columns:
    if col in matrix.index:
        matrix.loc[col, col] = 0
matrix = matrix.fillna(0)

# 2. Convert to a "Condensed Distance Matrix"
# This flattens the square grid into the specific format the maths requires
condensed_distances = squareform(matrix.values)

# 3. Calculate the Linkage (build the Tree)
# Can change 'average' to 'complete' or 'single' to see how the groupings shift.
print("calculating hierarchical linkage...")
my_linkage = linkage(condensed_distances, method='average') # use the pling distances for linkage

# 4. Draw the Heatmap
print("generating the figure...")
cg = sns.clustermap(
    matrix, 
    row_linkage=my_linkage,    
    col_linkage=my_linkage,     
    cmap='Reds_r', 
    annot=True, 
    fmt=".0f", 
    figsize=(12, 10),
    cbar_kws={'label': 'Pling Distance'},
    linewidths=.5
)
cg.cax.set_position([0.05, 0.175, 0.03, 0.2]) # move color map
# Clean up
#cg.ax_heatmap.set_title("Plasmid Pairwise Distances (UPGMA Clustering)", pad=20, fontsize=14) # title if wanted
cg.ax_heatmap.set_xlabel('')
cg.ax_heatmap.set_ylabel('')

# Save figure
cg.savefig('plasmid_clustermap_strict.svg',bbox_inches='tight')
print("Complete! Image saved as plasmid_clustermap_strict.png.")