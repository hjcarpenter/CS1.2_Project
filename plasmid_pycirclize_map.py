# -*- coding: utf-8 -*-
"""
Pycirclize mapping of plasmid1 from flagged annotation csv
Colour features by FLAG
Treat each row of the df as a feature

Created on Sat Dec  6 05:40:53 2025

@author: Heath
"""

import pandas as pd
from pycirclize import Circos
from Bio.SeqFeature import SeqFeature, FeatureLocation
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = 42  # TrueType fonts
plt.rcParams["ps.fonttype"] = 42

# IMPORT
fpp = 'D:/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all_fullmapUTF8.csv'
df = pd.read_csv(fpp, encoding='utf-8-sig')

# ---- normalize flag column
df["flag_norm"] = df["flag"].astype(str).str.strip().str.lower()

# ----------------------------
# Color mapping: EDIT THESE 
# ----------------------------
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

DEFAULT_COLOR = "#DDDDDD"

def color_for_flag(flag_norm):
    f = str(flag_norm).strip().lower() if pd.notna(flag_norm) else ""
    return FEATURE_COLORS.get(f, DEFAULT_COLOR)

# ----------------------------
# ---- Track routing rules (EDIT THESE)
# ----------------------------
TRACK_FLAGS = {
    "innermost": {"arr", "repregion", "traregion"},                    
    "inner": {"tn",  "cc1i"},
    "middle": {"integron", "is", "calin" },
    "outer": { "cds", "ta system", "virulence", "amr", "oriv", "orit"}
}

# flags NOT to plot (in any track)
EXCLUDE_FLAGS = {"promoter", "attc", "atti", "integrase", "sorf", "is_ir",
                 "ncrna", "ncrna-region", "is_dr", "iscr1", "repcds", "cn"}

# default behaviour for flags not listed
DEFAULT_TRACK = "inner"

# label behaviour per track
LABEL_TRACKS = {"outer"}  # only outer gets labels

# ---- Legend configuration
#TRACK_ORDER = ["innermost", "inner", "middle", "outer"] 
TRACK_ORDER = ["outer","middle","inner","innermost"  ]
TRACK_LABELS = {
    "innermost": "Innermost track", 
    "inner": "Inner track",
    "middle": "Middle track",
    "outer": "Outer track",
}

# pretty display names (optional but recommended --> edit these)
FLAG_DISPLAY = {
    "is": "IS Element",
    "amr": "Antimicrobial Resistance Gene",
    "ta system": "Toxin-Antitoxin System",
    "tn": "Transposon",
    "tu": "Putative mobile unit",
    "cn": "Composite Transposon",
    "integron": "Integron",
    "calin": "Cassette Array Lacking Integron",
    "cds": "CDS",
    "ncrna": "ncRNA",
    "oriv": "OriV",
    "orit": "OriT",
    "virulence": "Virulence Factor",
    "arr": "Resistance Region",
    "cc1i": "Complex Class I Integron", 
    "repregion": "Replicon Region", 
    #"repcds": "Replication Initiator Protein",
    "traregion": "Transfer Region"
}

# -----

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

# ---- PLASMID sectors
seqid2size = df.groupby("sequence_id")["end"].max().to_dict()
space = 0 if len(seqid2size) == 1 else 2
circos = Circos(sectors=seqid2size, space=space)
circos.text("pCS1.2IncF-NDM\n161.4 kb", size=20)

def df_to_features(sub_df):
    features = []
    for _, row in sub_df.iterrows():
        start = int(row["start"])
        end = int(row["end"])

        # strand
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

# ---- plot
for sector in circos.sectors:
    seqid = sector.name
    sub_df = df[df["sequence_id"] == seqid].copy()
    sub_df["flag"] = sub_df["flag"].astype(str)

# 1) Coordinate/tick track (INSIDE innermost)
    tick_track = sector.add_track((70, 72))
    tick_track.axis(fc="#FFFFFF", ec="none")
    tick_track.xticks_by_interval(
        interval=10000,
        outer=False,
        label_formatter=lambda v: f"{v / 1000:.1f} Kb",
        label_orientation="vertical",
        line_kws=dict(ec="grey"),
        text_kws=dict(size=12, color="black")
    )

    # 2) Innermost feature track
    innermost_track = sector.add_track((72, 75))
    innermost_track.axis(fc="#FFFFFF", ec="lightgrey")

    # 3) Inner feature track 
    inner_track = sector.add_track((75, 80))
    inner_track.axis(fc="#FFFFFF", ec="lightgrey")

    # 4) Middle feature track 
    middle_track = sector.add_track((80, 85))
    middle_track.axis(fc="#FFFFFF", ec="lightgrey")

    # 5) Outer feature track 
    outer_track = sector.add_track((85, 100))
    outer_track.axis(fc="#F5F5F5", ec="none")

    # split df into groups by flag
    sub_df = df[df["sequence_id"] == seqid].copy()
    sub_df["flag_norm"] = sub_df["flag"].astype(str).str.strip().str.lower()
    sub_df["track"] = sub_df["flag_norm"].apply(route_track)

    innermost_df = sub_df[sub_df["track"] == "innermost"] 
    inner_df     = sub_df[sub_df["track"] == "inner"]
    middle_df    = sub_df[sub_df["track"] == "middle"]
    outer_df     = sub_df[sub_df["track"] == "outer"]
    # skipped rows are simply not drawn

    # helper to draw a dataframe on a given track
    def draw_features(track, feat_df, do_labels=False, track_range=None):
        features = df_to_features(feat_df)
        for feature in features:
            flag = feature.qualifiers.get("flag", [None])[0]
            flag_norm = str(flag).strip().lower() if flag is not None else None

            color = color_for_flag(flag_norm)
            strand = feature.location.strand

            # 1. Determine Plot Style: BOX, or ARROW/ alpha
            if flag_norm in {"orit", "oriv", "cn", "is", "integron", "tu", "tn",
                             "cc1i", "calin", "arr", "repregion", "traregion"
                             }:
                plotstyle = "box"
                alpha = 1.0

            else:
                if flag_norm in {"cds"}:
                    plotstyle = "arrow"
                    alpha = 0.5
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
            track.annotate(label_pos, label, label_size=10)
            label_pos = (int(feature.location.start) + int(feature.location.end)) / 2

            # stagger label option
            # label_pos = (int(feature.location.start) + int(feature.location.end)) / 2
            
            # # Create a toggle to alternate label heights
            # if not hasattr(draw_features, "toggle"):
            #     draw_features.toggle = False
            # draw_features.toggle = not draw_features.toggle
            
            # # Alternate the maximum radius between 120 and 140
            # stagger_max = 106 if draw_features.toggle else 102
            
            # track.annotate(
            #     label_pos, 
            #     label, 
            #     label_size=8, 
            #     min_r=100.5, 
            #     max_r=stagger_max, 
            #     line_kws=dict(ec="grey", lw=0.5, alpha=0.8)
            # )
            
    # ---> EXECUTE THE DRAWING FOR EACH TRACK <---
    draw_features(innermost_track, innermost_df, do_labels=False)
    draw_features(inner_track, inner_df, do_labels=False)
    draw_features(middle_track, middle_df, do_labels=False)
    draw_features(outer_track, outer_df, do_labels=True)
        
# set up figure size
fig = circos.plotfig(figsize=(16, 16))

# set up legend
legend_handles = []

for track in TRACK_ORDER:

    flags_in_track = sorted(set(
        df.loc[df["flag_norm"].apply(route_track) == track, "flag_norm"].dropna()
    ).intersection(plotted_flags))

    if len(flags_in_track) == 0:
        continue

    # --- track label (option)
    legend_handles.append(
        Line2D([], [], linestyle="none", label=TRACK_LABELS.get(track, track))
    )

    # --- flags sorted alphabetically
    for flag in sorted(flags_in_track):
        legend_handles.append(
            Patch(
                facecolor=color_for_flag(flag),  # <-- legend uses same fixed mapping
                edgecolor="none",
                label=FLAG_DISPLAY.get(flag, flag),
            )
        )

    # spacer line
    legend_handles.append(Line2D([], [], linestyle="none", label=" "))

fig.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1.05, 0.5),
    fontsize=12,
    frameon=False,
    handlelength=1.0,
    handletextpad=0.6,
)

fig.tight_layout()
fig.savefig("IncF_plasmid_map_edit5.pdf", bbox_inches="tight")
fig.savefig("IncF_plasmid_map_edit5.svg", format="svg", bbox_inches="tight")
