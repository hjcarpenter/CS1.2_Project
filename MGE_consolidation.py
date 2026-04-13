# -*- coding: utf-8 -*-
"""
Mobile genetic elements consolidation file
Combines the cleaned, filtered results from ISFinder, ISEScan, Mobile Element Finder and integronfinder

Merges data from exact positional duplicates.
Generates consolidated intervals for overlapping and near overlapping (+/-5bp) alignments of a given [contig, subject] and collapses the results on these intervals to remove remaining redundancy while preserving fragmented elements.
Creates contig-specific subdataframes and saves the full and contig specific results to csv. 

Created on Sat Nov 22 09:24:04 2025

@author: Heath
"""

import pandas as pd
import numpy as np

# Import csvs
fp1 = 'ISEScan/ISEScan_results_filtered.csv'     	# ISEScan
fp2 = 'ISFinder/ISFinder_BLAST_deduped.csv' 		# ISFinder
fp3 = 'MobileElementFinder/MobileElementFinder_filtered.csv' 	# Mobile element finder
fp4 = 'IntegronFinder/Integron_Finder_results.csv' 		# Integron finder

ISCAN = pd.read_csv(fp1)
ISF = pd.read_csv(fp2)
MEF = pd.read_csv(fp3)
IntF = pd.read_csv(fp4)

# add name col to isfinder
ISF['name'] = ISF['subject_id']

# just want full integron span for Integron Finder plotting on map. 
# subset and collapse dupes
print(IntF.info())
cols = ["id_integron", "contig", "element_type", "source", "integron_start", "integron_end"]
new_intf = IntF[cols].drop_duplicates().copy()
new_intf = new_intf.rename(columns={"element_type": "type", 
                                    "integron_start": "start",
                                    "integron_end": "end", 
                                    "id_integron": "name"})

# choose columns to keep in the combined df
keep_cols = [
    "contig", "name", "type", "start", "end",
    "score", "evalue", "identity", "coverage", "length",
    "gaps", "synonyms", "prediction method",
    "source", "is_family", "tnp_orf_start", "tnp_orf_end", "ir_start1", "ir_end1",
    "ir_start2", "ir_end2"
    
]

def clean_df(df, keep):
    return df[[c for c in keep if c in df.columns]].copy()

ISCAN_clean = clean_df(ISCAN, keep_cols)
ISF_clean   = clean_df(ISF, keep_cols)
MEF_clean   = clean_df(MEF, keep_cols)
IntF_clean  = clean_df(new_intf, keep_cols)

# concatenate : 
combined = pd.concat(
    [ISCAN_clean, ISF_clean, MEF_clean, IntF_clean],
    ignore_index=True,
    sort=False
)

# recast floats to ints for coordinate cols
coord_cols = [
    "start", "end",
    "tnp_orf_start", "tnp_orf_end",
    "ir_start1", "ir_end1",
    "ir_start2", "ir_end2",
    "integron_start", "integron_end", 
    "start_consolidated",  "end_consolidated"
]

for col in coord_cols:
    if col in combined.columns:
        combined[col] = (
            pd.to_numeric(combined[col], errors="coerce")  # parse floats / strings / NaN
              .astype("Int64")                             # nullable integer dtype
        )


# --- Exact positional duplicate handling

group_cols = ["contig", "start", "end"]

def collapse_positional_group(g: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse rows with identical (contig, start, end).

    - If an 'isescan' row exists, use it as the base.
    - If a non-isescan row exists, use its 'name'.
    - Combine sources.
    - For other columns, prefer non-NA values (from isescan row first, then others).
    """
    # single row --> nothing to do
    if len(g) == 1:
        return g

    # choose base row: prefer isescan, else first row
    if (g["source"] == "isescan").any():
        base = g[g["source"] == "isescan"].iloc[0].copy()
    else:
        base = g.iloc[0].copy()

    # --- set name to the non-isescan name if present ---
    other_names = (
        g.loc[g["source"] != "isescan", "name"]
        .dropna()
        .unique()
    )
    if len(other_names) > 0:
        # if multiple different names, you could join them instead
        base["name"] = other_names[0]

    # --- combine sources into a semicolon-separated list ---
    base["source"] = ";".join(sorted(g["source"].dropna().unique()))

    # --- fill missing values in other columns from any row in group ---
    for col in g.columns:
        if col in ["name", "source", "contig", "start", "end"]:
            continue

        if pd.isna(base[col]) or (isinstance(base[col], float) and np.isnan(base[col])):
            # take first non-null value in this group for that column
            non_null_vals = g[col].dropna()
            if len(non_null_vals) > 0:
                base[col] = non_null_vals.iloc[0]

    # return as a one-row DataFrame
    return base.to_frame().T

# collapse isescan associated positional dupes
collapsed = (
    combined
    .groupby(group_cols, as_index=False, group_keys=False)
    .apply(collapse_positional_group)
    .reset_index(drop=True)
)

# check for remaining positional dupes
pos_cols = ["contig", "start", "end"]

# Check: All remaining positional duplicates 
pos_dupes = (
    collapsed[collapsed.duplicated(subset=pos_cols, keep=False)]
    .sort_values(pos_cols)
    .reset_index(drop=True)
)

# --- find overlapping alignmnets of the SAME CONTIG, NAME -------

def flag_overlaps_by_contig_name(df):
    """
    Add a boolean 'overlap' column:
    True if an interval overlaps at least one other interval
    with the same (contig, subject).
    """
    df = df.sort_values(['contig', 'name', 'start', 'end']).copy()
    df['overlap'] = False

    for (contig, subject), grp in df.groupby(['contig', 'name'], sort=False):
        if len(grp) < 2:
            continue

        idx = grp.index.to_numpy()
        starts = grp['start'].to_numpy()
        ends   = grp['end'].to_numpy()

        # compare each interval with the previous one
        overlap_with_prev = starts[1:] <= ends[:-1]  # use < if you *don’t* want touching intervals counted

        # mark both sides of each overlapping pair as True
        df.loc[idx[1:][overlap_with_prev],  'overlap'] = True
        df.loc[idx[:-1][overlap_with_prev], 'overlap'] = True

    return df

flagged = flag_overlaps_by_contig_name(collapsed)
# re-sort
flagged = flagged.sort_values(['contig','name', 'start', 'end'])

# create combined intervals for the overlaps
# Add new columns for overall start stop: copy start and end
flagged['start_consolidated'] = flagged['start']
flagged['end_consolidated'] = flagged['end']

# use fuzzy window and interval clustering logic to merge overlapping pairs but keep fragments separate

def consolidate_jitter(flagged, window_bp=5):
    """
    For each (contig, subject) group where overlap == True, cluster hits whose
    intervals are overlapping or within a small gap (window_bp).
    
    For each cluster, set:
        start_consolidated = min(start in cluster)
        end_consolidated   = max(end in cluster)

    Non-overlapping rows keep their original start/end.
    """
    df = flagged.copy()

    # initialise consolidated to original coords if not already present
    if 'start_consolidated' not in df.columns:
        df['start_consolidated'] = df['start']
    if 'end_consolidated' not in df.columns:
        df['end_consolidated'] = df['end']

    # work only on overlapping rows
    mask_overlap = df['overlap'] == True
    df_overlap = df[mask_overlap]

    # process per contig+subject
    for (contig, subject), grp in df_overlap.groupby(['contig', 'name'], sort=False):
        # sort by start so we can sweep
        grp = grp.sort_values(['start', 'end'])
        idxs = grp.index.tolist()
        if not idxs:
            continue

        # start first cluster
        current_cluster = [idxs[0]]
        cluster_start = grp.loc[idxs[0], 'start']
        cluster_end   = grp.loc[idxs[0], 'end']

        for i in idxs[1:]:
            s = grp.loc[i, 'start']
            e = grp.loc[i, 'end']

            # interval-overlap / near-overlap condition:
            # if the start of this interval is inside or very close to the current cluster
            # (allowing a small gap of window_bp), treat as same fragment
            if s <= cluster_end + window_bp:
                # same cluster – expand bounds
                current_cluster.append(i)
                cluster_start = min(cluster_start, s)
                cluster_end   = max(cluster_end,   e)
            else:
                # close previous cluster
                df.loc[current_cluster, 'start_consolidated'] = cluster_start
                df.loc[current_cluster, 'end_consolidated']   = cluster_end

                # start new cluster
                current_cluster = [i]
                cluster_start = s
                cluster_end   = e

        # final cluster in this group
        df.loc[current_cluster, 'start_consolidated'] = cluster_start
        df.loc[current_cluster, 'end_consolidated']   = cluster_end

    return df

# merge overlapping pairs in contig, subject groups. 
final = consolidate_jitter(flagged, window_bp=5)

# now collapse overlapping pairs on the intervals 

def collapse_overlaps_with_isescan(df):
    """
    Collapse overlapping hits with the same (contig, name, start_consolidated, end_consolidated).

    Rules:
    - Work only on rows where overlap == True.
    - Within each group, if any row has source containing 'isescan',
      use that row as the base.
    - Fill any remaining NaNs in the base row using non-null values
      from the other rows in the group.
    - 'source' becomes the concatenation of all unique sources in the group.
    """
    df = df.copy()

    key_cols = ['contig', 'name', 'start_consolidated', 'end_consolidated']

    # Split overlapping / non-overlapping
    overlap_df = df[df['overlap']]
    non_overlap_df = df[~df['overlap']]

    def collapse_group(g: pd.DataFrame) -> pd.Series:
        # if only one row, no need to do anything fancy
        if len(g) == 1:
            base = g.iloc[0].copy()
        else:
            # choose base row: prefer source containing 'isescan'
            isescan_mask = g['source'].str.contains('isescan', case=False, na=False)
            if isescan_mask.any():
                base = g[isescan_mask].iloc[0].copy()
            else:
                base = g.iloc[0].copy()

            # fill NaNs/empty strings in base from other rows
            for col in g.columns:
                if col == 'source':
                    continue  # handle separately below

                val = base[col]

                # treat NaN or empty string as "missing"
                if (isinstance(val, str) and val == '') or pd.isna(val):
                    non_null = g[col].dropna()
                    if len(non_null) > 0:
                        base[col] = non_null.iloc[0]

        # combine all sources in the group
        sources = g['source'].dropna().astype(str).unique()
        base['source'] = ';'.join(sources)

        return base

    # Collapse only overlapping groups (by consolidated interval)
    collapsed_overlap = (
        overlap_df
        .groupby(key_cols, as_index=False, group_keys=False)
        .apply(collapse_group)
    )

    # Put back non-overlap rows unchanged
    result = pd.concat([non_overlap_df, collapsed_overlap], ignore_index=True)

    # Optional: sort again for readability
    result = result.sort_values(['contig', 'name', 'start_consolidated', 'end_consolidated']).reset_index(drop=True)

    return result

collapsed_final = collapse_overlaps_with_isescan(final)

# create 2 subdfs - one for each contig
plasmid1_df = collapsed_final[collapsed_final["contig"] == "plasmid1"].copy()
# sort by start end
plasmid1_df = plasmid1_df.sort_values(["start", "end"]).reset_index(drop=True)
print(plasmid1_df.info())
chr_df = collapsed_final[collapsed_final["contig"] == "chromosome"].copy()
# sort by start end
chr_df = chr_df.sort_values(["start", "end"]).reset_index(drop=True)

# uniques
unique = plasmid1_df['name'].nunique()

# save
collapsed_final.to_csv('MGE_combined_results_chr_p1.csv', index=False) #all
chr_df.to_csv('MGE_combined_results_chromosome.csv', index=False) # chromosome only
plasmid1_df.to_csv('MGE_combined_results_plasmid1.csv', index=False) # plasmid only



