# -*- coding: utf-8 -*-
"""
Linear DNAFeatures viewer plot for a specified region: 
Plots up to 3 stacked tracks with configurable routing by flag
Modifiable legend labels
Modifiable feature colours
Option to show specified features as boxes rather than arrows
Input: flagged annotation file (csv) with at least: start, end, flag and name cols
Output: png, pdf and svg files. 
"""

import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

# ----------------------------
# --- Settings ---
# ----------------------------

# contig to display:
contig = 'plasmid'

# flagged annotation csv file path
csv_file = r'D:/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all_fullmapUTF8.csv'

# ----------------------------
# region to plot: coordinates (bp): EDIT THESE
# Note: 
# plasmid arr1: 12931, 16760
# plasmid arr2: 30410, 57860

start = 12931
end   = 16760

# ----------------------------
# --- Track routing rules: EDIT THESE.. what do you want to plot?
# ----------------------------
# Use lower-case keys/values for flags bc we normalise flags to flag_norm (lc)

TRACK_FLAGS = {
    "top":    {"amr", "virulence", "ta system", "cds",  
               "integrase", "attc", "is_ir", "calin", "is"},       
    "middle": {"is",  "tn", "calin", "cc1i", "iscr1",
               "attc", "atti", "promoter", "integron", "is_ir"},              
    "bottom": {"arr"}
}

EXCLUDE_FLAGS = {"ncrna-region", "arr", "is_dr", "ncrna", "sorf","cn" }  # DON'T plot these flags

DEFAULT_TRACK = "top"  # track to send flags not listed above

# which tracks should get labels? 
LABEL_TRACKS = {"top", "middle"}     # or {"bottom"} etc.

# ordering of tracks in the figure (top->bottom): 
# ALSO decides HOW MANY tracks to plot: i.e., n tracks = len(TRACK_ORDER)
TRACK_ORDER = ["top"]

# thickness per track (optional, defaults to 12: adjust as needed)
TRACK_THICKNESS = {"top": 16, "middle": 16, "bottom": 10}

# label position per track (optional)
TRACK_LABEL_POS = {"top": "top", "middle": "top", "bottom": "center"}

# ----------------------------
# --- Manual Legend Labels 
# ----------------------------
# Map the normalised flag to the text you want in the legend

LEGEND_LABELS = {
    "amr": "Antimicrobial Resistance Gene",
    "virulence": "Virulence Factor",
    "is": "IS element",
    "integron": "Integron",
    "ta system": "Toxin-Antitoxin System",
    "cds": "CDS",
    "cc1i": "Complex Class I Integron",
    "integrase": "Integrase",
    "tn": "Transposon",
    "promoter": "Promoter",
    "iscr1": "ISCR Element",
    "cn": "Composite Transposon (hypothetical)",
    "calin": "Cassette Array Lacking Integron",
    "spacer": "Spacer region",
    "is_ir": "Inverted Repeat",
    "is_dr": "Direct Repeat",
    "attc": "attC", 
    "atti":"attI"

}

# ----------------------------
# Color mapping to hex codes (edit as needed)
# ----------------------------
FEATURE_COLORS = {
    
    # ---- Mobile / recombination / integration

    "is": "#009E73",          # Okabe-Ito teal  **
    "is_ir": "#F78FE6",
    "is_dr": "#66CCEE",
    "spacer": "#E69F00",

    "integron": "#E69F00",    # Okabe-Ito orange **
    "cc1i": "#33BBEE",  # Tol cyan **
    "calin": "#332288",         # indigo
    "integrase": "#FFCCCC",   
    "attc": "#E69F00",        
    "atti": "#E69F00",        
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

DEFAULT_COLOR = "#DDDDDD" # if no colour specified

def color_for_flag(flag_norm):
    f = str(flag_norm).strip().lower() if pd.notna(flag_norm) else ""
    return FEATURE_COLORS.get(f, DEFAULT_COLOR)

# ----------------------------
# ----- Feature display
# ----------------------------

# elements to show as BOXES rather than arrows
BOX_FLAGS = {"orit", "oric", "cn", "is", "integron", "tn", "iscr1", "cc1i", 
             "calin", "arr","spacer", "is_ir", "is_dr", "attc"}  # shown as boxes (no direction)

# --- read + normalise columns ---

df = pd.read_csv(csv_file, encoding='utf-8-sig')
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

rename_map = {
    "begin": "start",
    "from": "start",
    "stop": "end",
    "to": "end",
}

df = df.rename(columns={c: rename_map[c] for c in df.columns if c in rename_map})

required = {"start", "end", "flag", "name"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

df["start"] = pd.to_numeric(df["start"], errors="coerce")
df["end"]   = pd.to_numeric(df["end"], errors="coerce")
df = df.dropna(subset=["start", "end"]).copy()

swap = df["start"] > df["end"]
df.loc[swap, ["start", "end"]] = df.loc[swap, ["end", "start"]].values

# normalise flag once
df["flag_norm"] = df["flag"].astype(str).str.strip().str.lower()

# --- filter to region (keep ALL overlapping features) ---
sub = df[(df["end"] >= start) & (df["start"] <= end)].copy()
sub["length"] = sub["end"] - sub["start"]

# --- Routing part ---

# normalise config sets 
TRACK_FLAGS_NORM = {t: {x.strip().lower() for x in s} for t, s in TRACK_FLAGS.items()}
EXCLUDE_FLAGS_NORM = {x.strip().lower() for x in EXCLUDE_FLAGS}

def route_track(flag_norm: str) -> str:
    if pd.isna(flag_norm):
        return DEFAULT_TRACK
    f = str(flag_norm).strip().lower()
    if f in EXCLUDE_FLAGS_NORM:
        return "skip"
    for track, flags in TRACK_FLAGS_NORM.items():
        if f in flags:
            return track
    return DEFAULT_TRACK

sub["track"] = sub["flag_norm"].apply(route_track)
sub = sub[sub["track"] != "skip"].copy()


def parse_strand(v):
    s = str(v).strip().lower() if pd.notna(v) else "+"
    if s in ("+", "1", "forward"):
        return 1
    if s in ("-", "-1", "reverse"):
        return -1
    return 0


# --- optional: source of label (name > gene >.. product, etc) ---
# (if you only have name, it just uses that)
def pick_label(row):
    for col in ("name", "gene"):
        val = row.get(col, None)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None

sub["label_for_plot"] = sub.apply(pick_label, axis=1)

# --- Build features per track ---
region_length = end - start

def build_features(df_track, thickness=12, label_position="center", show_labels=True):
    feats = []
    df_track = df_track.sort_values(["length", "start"], ascending=[False, True])

    for _, r in df_track.iterrows():
        a = max(int(r["start"]), start)
        b = min(int(r["end"]), end)
        if b <= a:
            continue

        f_start = a - start
        f_end   = b - start

        label = r["label_for_plot"] if show_labels else None

        # ---- style decisions 
        flag_norm = r.get("flag_norm", None)
        f = str(flag_norm).strip().lower() if pd.notna(flag_norm) else ""

        feat_color = color_for_flag(flag_norm)  # pass hex directly

        strand_val = parse_strand(r.get("strand", "+"))
        
        if f in BOX_FLAGS:
            strand_val = 0  # box (no arrow)

        feats.append(
            GraphicFeature(
                start=f_start,
                end=f_end,
                strand=strand_val,
                color=feat_color,
                label=label,
                label_position=label_position,
                thickness=thickness,
            )
        )

    return feats

records = {}
for track in TRACK_ORDER:
    tdf = sub[sub["track"] == track].copy()
    feats = build_features(
        tdf,
        thickness=TRACK_THICKNESS.get(track, 12),
        label_position=TRACK_LABEL_POS.get(track, "center"),
        show_labels=(track in LABEL_TRACKS),
    )
    records[track] = GraphicRecord(sequence_length=region_length, features=feats)

# --- plot stacked tracks (N tracks) ---

n = len(TRACK_ORDER)
height_ratios = [1] * n
fig, axes = plt.subplots(
    n, 1, figsize=(18, 4 * n), sharex=True, # adjust dimensions as needed #
    gridspec_kw={"height_ratios": height_ratios}
)
if n == 1:
    axes = [axes]

for i, track in enumerate(TRACK_ORDER):
    ax = axes[i]
    with_ruler = (i == n - 1)  # ruler only on bottom track
    records[track].plot(ax=ax, with_ruler=with_ruler, strand_in_label_threshold=7)
    ax.set_ylabel(track)

#axes[0].set_title(f"{start:,}–{end:,}bp") # settitle (or comment out)

# lots of x ticks...
#axes[-1].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{int(x + start):,}"))
# ... OR...
# Set ticks at the beginning (0) and the end (region_length) ONLY <- less messy
axes[-1].set_xticks([0, region_length])

# Label ticks with the actual genomic coordinates
axes[-1].set_xticklabels([f"{start:,}", f"{end:,}"])


# --- Build Custom Legend ---

# Identify which flags are in the plotted subset
used_flags = sorted(sub["flag_norm"].dropna().unique())

legend_handles = []
for flag in used_flags:
    color = color_for_flag(flag)
    
    # Check the dictionary for a manual label 
    # If not found, fall back to capitalized/title-case text
    if flag in LEGEND_LABELS:
        label_str = LEGEND_LABELS[flag]
    elif len(flag) <= 4:
        label_str = flag.upper()
    else:
        label_str = flag.title()
        
    patch = mpatches.Patch(color=color, label=label_str)
    legend_handles.append(patch)

# Add the legend to the overall figure at the bottom
# reduced ncol to 4 to give longer manual labels a bit more room
fig.legend(handles=legend_handles, 
           loc='upper center', 
           bbox_to_anchor=(0.5, 0.0), 
           ncol=min(len(legend_handles), 4), 
           frameon=False)

# --- Layout Adjustments ---
plt.tight_layout()

# NOTE: You may need to increase the bottom margin slightly (e.g., 0.20) 
# if manual labels wrap to multiple rows
fig.subplots_adjust(hspace=0.5, bottom=0.2) 

# --- Save options ----
plt.savefig(f"ARR1{contig}_{start}-{end}_linear_map.pdf", bbox_inches="tight")
# Save as SVG 
fig.savefig(f"ARR1{contig}_{start}-{end}_linear_map.svg", format="svg", bbox_inches="tight")

plt.show()
