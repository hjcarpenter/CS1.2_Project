# -*- coding: utf-8 -*-
"""
Pycirclize mapping of chromosome from flagged annotation csv
Colour features by FLAG
Treat each row of the df as a feature
Import Genbank file - Add: GC content and GC skew

@author: Heath
"""

import pandas as pd
import numpy as np
from pycirclize import Circos
from pycirclize.parser import Genbank
from Bio.SeqFeature import SeqFeature, FeatureLocation
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42  # TrueType fonts
plt.rcParams["ps.fonttype"] = 42

# ----------------------------
# INPUTS
# ----------------------------
# IMPORT flagged chromosome annotation
fpp = 'D:/CRE_PROJECT_DATA_ANALYSIS/All_Anno/chromosome/chromosome_flagged_all.csv'
df = pd.read_csv(fpp)

# IMPORT Genbank file (for GC)
fpg = 'D:/CRE_PROJECT_DATA_ANALYSIS/Bakta_compliant/bakta_out/per_contig/CS1_2_chromosome.gbff'
gbk = Genbank(fpg)

# Get seqs/sizes from GenBank
seqid2seq = gbk.get_seqid2seq()
seqid2size_gbk = gbk.get_seqid2size()

# Optional: restrict dataframe to seqids found in GenBank (prevents KeyError later)
df = df[df["sequence_id"].isin(seqid2size_gbk.keys())].copy()

# ---- normalise annotation flag column
df["flag_norm"] = df["flag"].astype(str).str.strip().str.lower()

# ----------------------------
# Color mapping EDIT THESE
# ----------------------------
# GC colors
skew_pos = "#EECC66"
skew_neg = "#997700"

content_pos = "#000000"
content_neg = "grey"

FEATURE_COLORS = {
    # ---- Mobile / recombination / integration

    "is": "#009E73",          # Okabe-Ito teal  **
    "integron": "#E69F00",    # Okabe-Ito orange **
    "cc1i": "#33BBEE",  # Tol cyan **
    "calin": "#332288",         # indigo
    "integrase": "#56B4E9",   # Okabe-Ito light blue **
    "attc": "#56B4E9",        # Okabe-Ito light blue **
    "atti": "#56B4E9",        # Okabe-Ito light blue **
    "orit": "#882255",          # TOL WINE

    # Inferred mobile regions
    "cn": "#CCDDAA",    # Tol pale green
    "tn": "#F0E442",    # Okabe-Ito yellow **
    "arr": "#cc79a7",    # Okabe-Ito pink  **
    
    # ---- Resistance & virulence
    "amr": "#D55E00",         # Okabe-Ito vermillion ***
    "virulence": "#0072B2",   #  Okabe-Ito blue*** 

    # ---- Defence systems
    "ta system": "#EE3377",       # Tol magenta *
    "crispr": "#FFCCCC",          # Tol pale red
    "crispr-repeat": "#FFCCCC",   # Tol pale red
    "crispr-spacer": "#FFCCCC",   # Tol pale red

    # ---- Coding sequences
    "cds": "#8f8f8f",         # Neutral gray
    "sorf": "#DDDDDD",        # Very light gray

    # ---- RNA features
    "ncrna": "#999933",           # Olive
    "ncrna-region": "#DDCC77",    # Light yellow/gold
    "trna": "#BBCC33",            # Lime/yellow-green
    "rrna": "#CCBB44",            # Gold
    "tmrna": "#AAAA00",           # Dark yellow

    # ---- Core & regulation
    "oriv": "#000000",        # Black (replication origin) ***
    "promoter": "#bbcc33",    # Tol PEAR
    "repcds": "#DDDDDD",
    "repregion": "#ddcc77",
    "traregion": "#dcbeff"
}

DEFAULT_COLOR = "#DDDDDD" # light grey

def color_for_flag(flag_norm):
    f = str(flag_norm).strip().lower() if pd.notna(flag_norm) else ""
    return FEATURE_COLORS.get(f, DEFAULT_COLOR)

# ---- Track routing rules 
TRACK_FLAGS = {
    "middle": {"virulence", "amr", "cds"},
    "outer": set(),
}

# flags NOT to plot (in any track)
EXCLUDE_FLAGS = {
    "promoter", "attc", "atti", "integrase", "sorf", "oric", "orit",
    "crispr", "crispr-repeat", "crispr-spacer", "inverted repeat",
    "is_ir", "ncrna", "ncrna-region", "rrna", "tmrna", "trna"
}

DEFAULT_TRACK = "skip"
LABEL_TRACKS = {"outer"}

TRACK_ORDER = ["middle", "outer"]
TRACK_LABELS = {"middle": "Middle track", "outer": "Outer track"}

# Labels to use in legend
FLAG_DISPLAY = {
    "is": "Insertion Sequence",
    "amr": "AMR",
    "ta system": "Toxin-Antitoxin system",
    "gi": "Genomic Island (predicted)",
    "cn": "Composite Tn (inferred)",
    "integron": "Integron",
    "cds": "CDS",
    "ncrna": "ncRNA",
    "oric": "oriC",
    "orit": "oriT",
    "virulence": "Virulence",
}

def route_track(flag: str) -> str:
    """Return which track a feature goes to: inner/middle/outer/skip."""
    if flag is None:
        return DEFAULT_TRACK
    f = str(flag).strip().lower()
    if f in EXCLUDE_FLAGS:
        return "skip"
    for track, flags in TRACK_FLAGS.items():
        if f in {x.lower() for x in flags}:
            return track
    return DEFAULT_TRACK

# Flags that will actually be drawn (i.e. not routed to "skip")
plotted_flags = sorted(
    df.loc[df["flag_norm"].apply(route_track) != "skip", "flag_norm"]
      .dropna()
      .unique()
)
print("Plotted flags:", len(plotted_flags), plotted_flags)

# ----------------------------
# SECTORS (use GenBank sizes)
# ----------------------------
seqid2size = seqid2size_gbk
space = 0 if len(seqid2size) == 1 else 2
circos = Circos(sectors=seqid2size, space=space)
circos.text("CS1.2 contig-1\n4.89 kb", size=17)

# ----------------------------
# helpers
# ----------------------------
def df_to_features(sub_df):
    features = []
    for _, row in sub_df.iterrows():
        start = int(row["start"])
        end = int(row["end"])

        s = row.get("strand", 0)
        if s in ["+", 1, "1", "+1"]:
            strand = 1
        elif s in ["-", -1, "-1"]:
            strand = -1
        else:
            strand = 0

        qualifiers = {}
        for k in ["gene", "name", "product", "flag"]:
            v = row.get(k, None)
            if pd.notna(v):
                qualifiers[k] = [str(v)]

        feat_type = str(row.get("type", "CDS"))

        feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type=feat_type,
            qualifiers=qualifiers,
        )
        features.append(feature)
    return features

def _split_rlim(track_range, strand):
    """Split a track range into two sub-bands by strand (+ on outer half)."""
    r0, r1 = track_range
    mid = (r0 + r1) / 2
    if strand == 1:
        return (mid, r1)
    else:
        return (r0, mid)

def draw_features(track, feat_df, do_labels=False, track_range=None):
    features = df_to_features(feat_df)
    for feature in features:
        flag = feature.qualifiers.get("flag", [None])[0]
        flag_norm = str(flag).strip().lower() if flag is not None else None

        color = color_for_flag(flag_norm)
        strand = feature.location.strand

        # 1. Determine Plot Style
        if flag_norm in {"orit", "oric", "cn", "is", "integron", "gi"}:
            plotstyle = "box"
            alpha = 0.5
        else:
            if flag_norm in {"cds"}:
                plotstyle = "arrow"
                alpha = 0.6
            else:
                plotstyle = "arrow"
                alpha = 1.0

        # 2. Determine Height (r_lim)
        #    FIX: Only split the track height if it's an ARROW (CDS). 
        #    If it's a BOX (Transposon/IS), use the full track range.
        if track_range is None:
            # Fallback if track_range not provided (uses track defaults)
            if plotstyle == "box":
                r_lim = track.r_lim
            else:
                r_lim = _split_rlim(track.r_lim, strand)
        else:
            # Use the provided track_range
            if plotstyle == "box":
                r_lim = track_range  # <--- Use FULL height for boxes
            else:
                r_lim = _split_rlim(track_range, strand) # <--- Keep split for arrows

        track.genomic_features(
            feature,
            plotstyle=plotstyle,
            r_lim=r_lim,
            fc=color,
            alpha=alpha
        )

        if not do_labels:
            continue

        label = feature.qualifiers.get("name", [None])[0]
        if not label or str(label).startswith("hypothetical"):
            continue

        label_pos = (int(feature.location.start) + int(feature.location.end)) / 2
        track.annotate(label_pos, label, label_size=8)

# ----------------------------
# PLOT
# ----------------------------
# Baseline GC for delta-GC content (whole genome)
genome_gc = gbk.calc_genome_gc_content(seq=gbk.full_genome_seq)

for sector in circos.sectors:
    seqid = sector.name
    seq = seqid2seq.get(seqid, None)
    if seq is None:
        # skip sectors not present in genbank (shouldn't happen if filtered)
        continue

    sub_df = df[df["sequence_id"] == seqid].copy()
    sub_df["flag_norm"] = sub_df["flag"].astype(str).str.strip().str.lower()
    sub_df["track"] = sub_df["flag_norm"].apply(route_track)

    #inner_df  = sub_df[sub_df["track"] == "inner"] # commented out, not plotting
    middle_df = sub_df[sub_df["track"] == "middle"]
    outer_df  = sub_df[sub_df["track"] == "outer"]

    # ---- GC skew track (innermost)
    gc_skew_track = sector.add_track((60, 70))
    gc_skew_track.axis(fc="#FFFFFF", ec="none")

    pos_list, gc_skews = gbk.calc_gc_skew(seq=seq)
    pos_list = np.asarray(pos_list)
    gc_skews = np.asarray(gc_skews)

    pos_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
    neg_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
    abs_max_gc_skew = np.max(np.abs(gc_skews)) if len(gc_skews) else 1.0
    vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew

    gc_skew_track.fill_between(pos_list, pos_gc_skews, 0, 
                               vmin=vmin, vmax=vmax, color=skew_pos)
    gc_skew_track.fill_between(pos_list, neg_gc_skews, 0, 
                               vmin=vmin, vmax=vmax, color=skew_neg)

    # ---- GC content track (delta from genome mean)
    gc_content_track = sector.add_track((70, 80))
    gc_content_track.axis(fc="#FFFFFF", ec="none")

    pos_list, gc_contents = gbk.calc_gc_content(seq=seq)
    pos_list = np.asarray(pos_list)
    gc_contents = np.asarray(gc_contents)

    # delta from genome mean GC
    gc_contents = gc_contents - genome_gc

    pos_gc = np.where(gc_contents > 0, gc_contents, 0)
    neg_gc = np.where(gc_contents < 0, gc_contents, 0)
    abs_max_gc = np.max(np.abs(gc_contents)) if len(gc_contents) else 1.0
    vmin, vmax = -abs_max_gc, abs_max_gc

    gc_content_track.fill_between(pos_list, pos_gc, 0, vmin=vmin, 
                                  vmax=vmax, color=content_pos)
    gc_content_track.fill_between(pos_list, neg_gc, 0, vmin=vmin, 
                                  vmax=vmax, color=content_neg)

    # 2) Inner feature track
    # inner_track = sector.add_track((68, 80))
    # inner_track.axis(fc="#F5F5F5", ec="none")

    # 3) Middle feature track
    middle_track = sector.add_track((80, 88))
    middle_track.axis(fc="#FFFFFF", ec="none")

    # 4) Outer feature track

    # A. Define the Feature Track FIRST
    #    RANGE: 88 to 98.5 (extend it by 0.5 to overlap into the grey ring area)
    feat_outer_track = sector.add_track((88, 98.5), r_pad_ratio=0.0)
    feat_outer_track.axis(fc="#FFFFFF", ec="none")

    # B. Define the Coordinate (Grey) Track SECOND
    #    RANGE: 98 to 100. Because this is added second, it draws ON TOP.
    #    It will cover the extra 0.5 overlap of the genes, creating a perfect seal.
    outer_track = sector.add_track((90, 90.5), r_pad_ratio=0.0)
    outer_track.axis(fc="lightgrey", ec="none")

    # C. Draw the Ticks on the Grey Track
    outer_track.xticks_by_interval(
        interval=250000,
        outer=True,                     # Ticks point OUT from the ring
        tick_length=2,
        label_margin=0.5,
        label_size=10,
        label_orientation="horizontal",   # Radial text to prevent "floating" labels
        label_formatter=lambda v: f"{v / 1000:.0f} Kb",
        line_kws=dict(ec="black", lw=0.8),
    )

    # D. Draw the Features
    #    IMPORTANT: Use the extended range (88, 98.5) for the calculation
    draw_features(middle_track,     middle_df, do_labels=("middle" in LABEL_TRACKS), track_range=(80, 88))
    draw_features(feat_outer_track, outer_df,  do_labels=("outer"  in LABEL_TRACKS), track_range=(88, 98.5))
# set up figure size
fig = circos.plotfig(figsize=(15, 15))

# ---- legend
legend_handles = []
for track in TRACK_ORDER:
    flags_in_track = sorted(set(
        df.loc[df["flag_norm"].apply(route_track) == track, "flag_norm"].dropna()
    ).intersection(plotted_flags))

    if len(flags_in_track) == 0:
        continue

    for flag in sorted(flags_in_track):
        legend_handles.append(
            Patch(facecolor=color_for_flag(flag), edgecolor="none", label=FLAG_DISPLAY.get(flag, flag))
        )

# ----------------------------
# Add GC Tracks to Legend
# ----------------------------
# Create a separator or header (optional), or just comment out and append directly
#legend_handles.append(Patch(color="none", label=" "))  # Spacer
#legend_handles.append(Patch(color="none", label="GC Tracks:")) # Header

# GC Skew Handles
legend_handles.append(Patch(facecolor= skew_pos, edgecolor="none", label="GC Skew (+)"))
legend_handles.append(Patch(facecolor= skew_neg, edgecolor="none", label="GC Skew (-)"))

# GC Content Handles
legend_handles.append(Patch(facecolor=content_pos, edgecolor="none", label="GC Content (+)"))
legend_handles.append(Patch(facecolor=content_neg, edgecolor="none", label="GC Content (-)"))

fig.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1.05, 0.5),
    fontsize=12,
    frameon=False,
    #handlelength=1.0,
    #handletextpad=0.6,
)

#fig.tight_layout()
fig.savefig("cs1.2_map_GC.png", dpi=300, bbox_inches="tight")
fig.savefig("cs1.2_map_GC.pdf", bbox_inches="tight")
fig.savefig("cs1.2_map_GC.svg", bbox_inches="tight")
