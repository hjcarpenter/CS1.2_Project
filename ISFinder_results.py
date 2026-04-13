# -*- coding: utf-8 -*-
"""
This script cleans and filters the IS Finder BLASTn results
1. Import ISFinder BLAST results as TSV and add headers
2. Filter exact coordinate duplicate start/end hits
3. Filter the remaining results for:
evalue <= 1e-10 (should all be below this anyway)
identity >=35
alignment length >=250
bitscore >=50
4. collapse overlapping and nearby (<=50bp) hits into non-redundant loci: keep best hit per locus ["bitscore"
>"identity">"alignment_length"]
5. save

Created on Tue Nov 25 04:54:55 2025

@author: Heath
"""

import pandas as pd

# ADD PATH TO TSV
fp = 'ISFinder_Blast_result.tsv'

# column headers
cols = [
    "query id", "subject id", "% identity", "alignment length",
    "mismatches", "gap opens", "q. start", "q. end", "s. start", 
    "s. end", "evalue", "bit score"
]

# read in tsv with the standard blast headers
df = pd.read_csv(fp, sep='\t', names=cols)

# optional: rename columns for downstream work
rename_map = {
    "query id": "contig",
    "subject id":"subject_id",
    "% identity": "identity",
    "q. start": "start",
    "q. end": "end",
    "s. start":"s_start",
    "s. end":"s_end",
    "gap opens":"gap_opens",
    "bit score": "bitscore",
    "alignment length": "alignment_length"
}

df = df.rename(columns=rename_map)

# add a source column 
df["source"] = "isfinder"
# add a type column
df["type"] = "insertion sequence"

# save the raw results dataframe to csv 
df.to_csv('ISFinder_BLAST_Raw.csv', index=False)
print(df.info())

# ---------- Drop EXACT coordinate duplicates

# drop duplicate start/end hits using a heirarchical tie-break
group_cols = ['contig', 'start', 'end']   

# priority order: max bitscore > max identity > min mismatch
blast_filt = (
    df.sort_values(
        by=['bitscore', 'identity', 'mismatches'],
        ascending=[False, False, True]     
    )
    .drop_duplicates(subset=group_cols, keep='first')
)

# Now filter: evalue <= 1e-10 (should all be below this anyway)
# identity >=35
# alignment length >=250
# bitscore >=50
blast_filt = blast_filt[
    (blast_filt["evalue"] <= 1e-10) &
    (blast_filt["identity"] >= 30) &
    (blast_filt["alignment_length"] >= 250) &
     (blast_filt["bitscore"] >=50)
].copy()


# sort by contig, start, end
blast_filt = blast_filt.sort_values(["contig", "start", "end"]).reset_index(drop=True)

# collapse hits into non redundant loci keeping the best hit per locus
def collapse_IS_hits(df, max_gap=50):
    """
    Collapse overlapping/nearby ISFinder BLAST hits into non-redundant loci.
    Keep the best hit per locus (best = highest bitscore, then identity, then length).

    Adds:
        - locus_hits: string listing all hits (subject_id + key stats)
        - locus_n_hits: number of hits in the locus
    """
    # we assume columns:
    # ['contig', 'subject_id', 'identity', 'alignment_length',
    #  'start', 'end', 'evalue', 'bitscore']

    df = df.sort_values(["contig", "start", "end"]).reset_index(drop=True)
    keep = []
    locus_hits_list = {}
    locus_n_hits_list = {}

    for contig, sub in df.groupby("contig", sort=False):
        sub = sub.sort_values(["start", "end"]).copy()

        current_cluster = [sub.index[0]]
        current_end = sub.iloc[0]["end"]

        for idx, row in sub.iloc[1:].iterrows():
            s, e = row["start"], row["end"]

            # overlap or small gap → same cluster
            if s <= current_end + max_gap:
                current_cluster.append(idx)
                current_end = max(current_end, e)
            else:
                # finalize previous cluster
                cluster = df.loc[current_cluster]

                # pick best hit in this cluster
                best_idx = (
                    cluster.sort_values(
                        ["bitscore", "identity", "alignment_length"],
                        ascending=[False, False, False]
                    )
                    .index[0]
                )
                keep.append(best_idx)

                # build locus_hits string for this cluster
                hits_str = "; ".join(
                    f"{r.subject_id}"
                    f" (id={r.identity:.1f}%, len={int(r.alignment_length)}, "
                    f"bitscore={r.bitscore:.1f})"
                    for _, r in cluster.iterrows()
                )
                locus_hits_list[best_idx] = hits_str
                locus_n_hits_list[best_idx] = len(cluster)

                # start new cluster
                current_cluster = [idx]
                current_end = e

        # finalize last cluster for this contig
        cluster = df.loc[current_cluster]
        best_idx = (
            cluster.sort_values(
                ["bitscore", "identity", "alignment_length"],
                ascending=[False, False, False]
            )
            .index[0]
        )
        keep.append(best_idx)

        hits_str = "; ".join(
            f"{r.subject_id}"
            f" (id={r.identity:.1f}%, len={int(r.alignment_length)}, "
            f"bitscore={r.bitscore:.1f})"
            for _, r in cluster.iterrows()
        )
        locus_hits_list[best_idx] = hits_str
        locus_n_hits_list[best_idx] = len(cluster)

    collapsed = df.loc[keep].sort_values(["contig", "start"]).reset_index(drop=True)

    # map the locus info onto the collapsed df
    collapsed["locus_hits"] = collapsed.index.map(
        lambda i: locus_hits_list.get(collapsed.index[i])
    )
    collapsed["locus_n_hits"] = collapsed.index.map(
        lambda i: locus_n_hits_list.get(collapsed.index[i])
    )

    return collapsed

blast_nonredundant = collapse_IS_hits(blast_filt, max_gap=50)

blast_nonredundant[[
    "contig", "start", "end", "subject_id",
    "identity", "alignment_length", "bitscore",
    "locus_n_hits", "locus_hits"
]].head()

# save
blast_nonredundant.to_csv('ISFinder_BLAST_deduped.csv', index=False)



