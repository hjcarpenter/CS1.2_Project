# -*- coding: utf-8 -*-
"""
This script generates and plots a MUMmer alignment for two genbank files (reference and query) using PyCirclize
Plots arrows for AMR and virulence genes on the reference 
# (requires an annotation CSV with a 'flag' column containing 'AMR' and 'virulence' flagged entries, 
# 'start'/ 'end' coordinate columns and 'strand' (+ or -) for direction   
Allows you to:
Edit filepaths for different alignments. 
Edit Display names
Edit minimum alignment length to display (higher to show higher relatedness)
MODIFY filepaths/ params as needed...
THEN
RUN ON COMMAND LINE:
1. Navigate to directory holding this script    
2. conda activate circos-mummer
3. python mummer_plot_tune.py  # run script
4. conda deactivate

Created on Sun Feb  1 03:19:27 2026

@author: Heath
"""

from pathlib import Path
from pycirclize import Circos
from pygenomeviz.parser import Genbank
from pygenomeviz.align import MUMmer
import pandas as pd
from matplotlib.patches import Patch

# -----------------------
# Inputs: Two filepaths for the genbank files to use as ref and query
# -----------------------
fp_1 = "/mnt/d/CRE_PROJECT_DATA_ANALYSIS/Bakta_compliant/bakta_out/per_contig/CS1_2_plasmid1.gbff"  # REFERENCE
fp_2 = "/mnt/d/CRE_PROJECT_DATA_ANALYSIS/Other_strains/zurfluh_strain_675SK2/plasmid/bakta_out/bakta_675SK2.gbff" # QUERY

DISPLAY_NAME = {
    fp_1: "pCS1.2IncF-NDM", # REFERENCE
    fp_2: "p675SK2_B", # QUERY
}

# --- ANNOTATION CSV file to show AMR/Virulence arrows on the reference ---
ANNOT_CSV = "/mnt/d/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all_fullmapUTF8.csv"  # <-- set filepath to annotation csv
ARROW_TRACK = (95.0, 99.0)    # near the outside edge
ARROW_R_PAD = 0.05
ARROW_EC = None               # no outline

# ---- LEGEND Colours and Labels ----
legend_handles = [
    Patch(fc="#E69F00", label="Forward synteny"),
    Patch(fc="#cc79a7", label="Inverted synteny"),
    Patch(fc="#D55E00", label="AMR gene"),
    Patch(fc="#0072B2", label="Virulence gene"),
]

# -----------------------
# Tuning parameters
# -----------------------
TICKS_INTERVAL = 10_000

MIN_ALIGN_LEN = 3000   # Smallest alignment length to display
MAX_LINKS = None       # Maximum number of alignments to display (hard cap if needed)

FIGSIZE = (11, 11)
START_DEG = -358
END_DEG = 2
SPACE_DEG = 6         # breathing room between seqs

LINK_R = 94
AXIS_TRACK = (99.6, 100)

# Synteny Colors <-- same as legend colors
FWD_COLOR = "#E69F00"
INV_COLOR = "#cc79a7"

# -----------------------
# Load GenBank
# -----------------------
for fp in (fp_1, fp_2):
    if not Path(fp).exists():
        raise FileNotFoundError(fp)

ref_gbk = Genbank(fp_1)
query_gbk = Genbank(fp_2)

ref_name = DISPLAY_NAME.get(fp_1, ref_gbk.name)
query_name = DISPLAY_NAME.get(fp_2, query_gbk.name)

# -----------------------
# Circos layout (two sectors)
# -----------------------
ref_seqid2size = ref_gbk.get_seqid2size()
query_seqid2size = query_gbk.get_seqid2size()

sectors = dict(**ref_seqid2size, **dict(reversed(list(query_seqid2size.items()))))
sector2clockwise = {seqid: False for seqid in query_seqid2size.keys()}

circos = Circos(
    sectors=sectors,
    start=START_DEG,
    end=END_DEG,
    space=SPACE_DEG,
    sector2clockwise=sector2clockwise,
)

ref_kb = round(ref_gbk.full_genome_length / 1000, 1)
query_kb = round(query_gbk.full_genome_length / 1000, 1)

circos.text(f"{ref_name}\n({ref_kb} kb)", r=132, deg=35, size=13)
circos.text(f"{query_name}\n({query_kb} kb)", r=132, deg=-35, size=13)

# Axis + ticks
for sector in circos.sectors:
    track = sector.add_track(AXIS_TRACK)
    track.axis(fc="black")
    if sector.size >= TICKS_INTERVAL:
        track.xticks_by_interval(
            TICKS_INTERVAL,
            label_formatter=lambda v: f"{int(v/1000)} kb",
            label_orientation="vertical",
        )

# --- Load annotation CSV ---
df = pd.read_csv(ANNOT_CSV)

# Filter AMR + virulence
flag_s = df["flag"].astype(str)
mask_amr = flag_s.str.contains("amr", case=False, na=False)
mask_vir = flag_s.str.contains("virulence", case=False, na=False)

feat = df[(mask_amr | mask_vir)].copy()
feat["is_amr"] = feat["flag"].astype(str).str.contains("amr", case=False, na=False)

# IMPORTANT: map CSV sequence_id to the reference sector name(s)
ref_seqids = set(ref_gbk.get_seqid2size().keys())
feat = feat[feat["sequence_id"].isin(ref_seqids)].copy()

# Draw arrows ONLY on reference sector(s)
for sector in circos.sectors:
    if sector.name not in ref_seqids:
        continue  # skip query sector(s)

    track = sector.add_track(ARROW_TRACK, r_pad_ratio=ARROW_R_PAD)

    sdf = feat[feat["sequence_id"] == sector.name].sort_values("start")
    if sdf.empty:
        continue

    for _, r in sdf.iterrows():
        start = int(r["start"])
        end = int(r["end"])
        strand = str(r.get("strand", "+"))

        color = "#D55E00" if bool(r["is_amr"]) else "#0072B2"

        # --- Draw arrow ---
        if hasattr(track, "arrow"):
            # direction by strand
            forward = not strand.startswith("-")
            if forward:
                track.arrow(start, end, fc=color, ec=ARROW_EC)
            else:
                track.arrow(end, start, fc=color, ec=ARROW_EC)
        else:
            # Fallback: draw a rectangle if arrow() isn't available
            track.rect(start, end, color=color, ec=ARROW_EC)


# -----------------------
# MUMmer alignment + filtering
# -----------------------
align_coords = MUMmer([query_gbk, ref_gbk]).run()

def alen(ac):
    # robust length definition
    return max(abs(ac.query_end - ac.query_start), abs(ac.ref_end - ac.ref_start))

align_coords = [ac for ac in align_coords if alen(ac) >= MIN_ALIGN_LEN]

# Optional: sort by length and cap number of links (helps dense plots)
align_coords.sort(key=alen, reverse=True)
if MAX_LINKS is not None:
    align_coords = align_coords[:MAX_LINKS]

# -----------------------
# Links
# -----------------------
for ac in align_coords:
    region_query = (ac.query_name, min(ac.query_start, ac.query_end), max(ac.query_start, ac.query_end))
    region_ref = (ac.ref_name,   min(ac.ref_start, ac.ref_end),     max(ac.ref_start, ac.ref_end))

    color = INV_COLOR if ac.is_inverted else FWD_COLOR

    # alpha/lw make it prettier
    circos.link(region_query, region_ref, color=color, r1=LINK_R, r2=LINK_R, alpha=0.4, lw=0.6)

fig = circos.plotfig(figsize=FIGSIZE)

fig.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(-0.15, 0.5),  # move legend outside, left
    frameon=False,
)

# 2 save options: png and svg
fig.savefig(f"mummer2_{ref_name}_vs_{query_name}_min_align_len_{MIN_ALIGN_LEN}bp.png", dpi=300, bbox_inches="tight")
fig.savefig(f"mummer2_{ref_name}_vs_{query_name}_min_align_len_{MIN_ALIGN_LEN}bp.svg", bbox_inches="tight")
print("Saved: mummer tuned")
