# -*- coding: utf-8 -*-
"""
This script is an extension of mummer_plot_tune.py for a batch run on a folder containing multiple query genbank files.
Requires reference filepath, query folder filepath and filepath to annotation csv (flagged/ start/end coords/ strand required) for
amr and virulence gene arrow annotations on the reference. 
Generates individual Circos MUMmer plots for each reference/query pairing with PyCirclize
Also extracts metadata from Genbank files for plot annotation: NOTE! these can be messy!
All outputs save as png and svg (for downstream edit) to a plots folder inside the query folder

MODIFY filepaths/ params as needed 
THEN
RUN ON COMMAND LINE:
1. Navigate to directory holding this script E.G., cd /mnt/d/CRE_PROJECT_DATA_ANALYSIS/Other_strains/plasmid_comparison_plots/plasmid_homologues_mummer    
2. conda activate circos-mummer
3. python mummer_plot_batch.py # run it!
4. conda deactivate
"""

import os
import textwrap
from pathlib import Path
import matplotlib.pyplot as plt
from pycirclize import Circos
from pygenomeviz.parser import Genbank
from pygenomeviz.align import MUMmer
import pandas as pd
from matplotlib.patches import Patch

# -----------------------
# Inputs - EDIT THESE (WSL paths!)
# -----------------------
# Reference
REF_FP = Path("/mnt/d/CRE_PROJECT_DATA_ANALYSIS/Bakta_compliant/bakta_out/per_contig/CS1_2_plasmid1.gbff") # REFERENCE GBK FILEPATH HERE
REF_NAME = "pCS1.2IncF-NDM"

# Queries Directory
QUERY_DIR = Path("/mnt/d/CRE_PROJECT_DATA_ANALYSIS/Other_strains/plasmid_comparison_plots/plasmid_homologues_mummer") # QUERY FOLDER FILEPATH HERE

# Output Directory Setup
OUTPUT_DIR = QUERY_DIR / "plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True) # Creates the 'plots' folder if it doesn't exist

# Annotations
ANNOT_CSV = "/mnt/d/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all_fullmapUTF8.csv" # ANNOTATION CSV FILEPATH HERE

# -----------------------
# Tuning parameters
# -----------------------
TICKS_INTERVAL = 10_000
MIN_ALIGN_LEN = 3000   # minimum alignment length to plot
MAX_LINKS = None       # max number of alignments to plot (hard cap if needed)
FIGSIZE = (12, 12)     # Make larger to accommodate longer labels
START_DEG = -358
END_DEG = 2
SPACE_DEG = 6         
LINK_R = 94
AXIS_TRACK = (99.6, 100)
ARROW_TRACK = (95.0, 99.0)    
ARROW_R_PAD = 0.05
ARROW_EC = None               

# Plot synteny Colors
FWD_COLOR = "#E69F00"
INV_COLOR = "#cc79a7"

# ---- LEGEND colours and labels  ----
legend_handles = [
    Patch(fc="#E69F00", label="Forward synteny"), 
    Patch(fc="#cc79a7", label="Inverted synteny"),
    Patch(fc="#D55E00", label="AMR gene"),
    Patch(fc="#0072B2", label="Virulence gene"),
]

# -----------------------
# Pre-load Reference & Annotations (Outside Loop)
# -----------------------
if not REF_FP.exists():
    raise FileNotFoundError(f"Reference not found: {REF_FP}")

ref_gbk = Genbank(REF_FP)
ref_seqid2size = ref_gbk.get_seqid2size()
ref_seqids = set(ref_seqid2size.keys())
ref_kb = round(ref_gbk.full_genome_length / 1000, 1)

# Load CSV
df = pd.read_csv(ANNOT_CSV)
flag_s = df["flag"].astype(str)
mask_amr = flag_s.str.contains("amr", case=False, na=False)
mask_vir = flag_s.str.contains("virulence", case=False, na=False)

feat = df[(mask_amr | mask_vir)].copy()
feat["is_amr"] = feat["flag"].astype(str).str.contains("amr", case=False, na=False)
feat = feat[feat["sequence_id"].isin(ref_seqids)].copy()

# -----------------------
# Batch Processing Loop
# -----------------------
# Grab all .gb, .gbk, and .gbff files in the directory
query_files = list(QUERY_DIR.glob("*.gb*"))

if not query_files:
    print(f"No GenBank files found in {QUERY_DIR}")
else:
    print(f"Found {len(query_files)} query files. Starting batch processing...")

for query_fp in query_files:
    print(f"\nProcessing: {query_fp.name}...")
    
    # Load Query
    query_gbk = Genbank(query_fp)
    
    # --- Metadata Extraction ---
    # pygenomeviz parses into Biopython records under the hood
    record = query_gbk.records[0] 
    definition_line = record.description
    
    # Search for date and location in the 'source' features
    collection_date = "N/A"
    geo_loc_name = "N/A"
    
    for f in record.features:
        if f.type == "source":
            if "collection_date" in f.qualifiers:
                collection_date = f.qualifiers["collection_date"][0]
            if "geo_loc_name" in f.qualifiers:
                geo_loc_name = f.qualifiers["geo_loc_name"][0]
            break # Usually only one source feature per record
            
    # Wrap text so long definition lines don't run off the edge
    wrapped_def = "\n".join(textwrap.wrap(definition_line, width=40))
    query_display_label = f"{wrapped_def}\nLocation: {geo_loc_name} | Date: {collection_date}"

    # --- Circos Layout ---
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

    query_kb = round(query_gbk.full_genome_length / 1000, 1)

    circos.text(f"{REF_NAME}\n({ref_kb} kb)", r=132, deg=35, size=13)
    # Use the extracted metadata for the query text
    circos.text(f"{query_display_label}\n({query_kb} kb)", r=132, deg=-35, size=10)

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

    # --- Draw arrows on reference sector(s) ---
    for sector in circos.sectors:
        if sector.name not in ref_seqids:
            continue  

        track = sector.add_track(ARROW_TRACK, r_pad_ratio=ARROW_R_PAD)
        sdf = feat[feat["sequence_id"] == sector.name].sort_values("start")
        if sdf.empty:
            continue

        for _, r in sdf.iterrows():
            start, end = int(r["start"]), int(r["end"])
            strand = str(r.get("strand", "+"))
            color = "#D55E00" if bool(r["is_amr"]) else "#0072B2"

            if hasattr(track, "arrow"):
                forward = not strand.startswith("-")
                if forward:
                    track.arrow(start, end, fc=color, ec=ARROW_EC)
                else:
                    track.arrow(end, start, fc=color, ec=ARROW_EC)
            else:
                track.rect(start, end, color=color, ec=ARROW_EC)

    # --- MUMmer alignment ---
    align_coords = MUMmer([query_gbk, ref_gbk]).run()

    def alen(ac):
        return max(abs(ac.query_end - ac.query_start), abs(ac.ref_end - ac.ref_start))

    align_coords = [ac for ac in align_coords if alen(ac) >= MIN_ALIGN_LEN]
    align_coords.sort(key=alen, reverse=True)
    
    if MAX_LINKS is not None:
        align_coords = align_coords[:MAX_LINKS]

    # Links
    for ac in align_coords:
        region_query = (ac.query_name, min(ac.query_start, ac.query_end), max(ac.query_start, ac.query_end))
        region_ref = (ac.ref_name,   min(ac.ref_start, ac.ref_end),     max(ac.ref_start, ac.ref_end))
        color = INV_COLOR if ac.is_inverted else FWD_COLOR
        circos.link(region_query, region_ref, color=color, r1=LINK_R, r2=LINK_R, alpha=0.4, lw=0.6)

    # --- Plot & Save ---
    fig = circos.plotfig(figsize=FIGSIZE)

    fig.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(-0.15, 0.5), 
        frameon=False,
    )

    # Clean filename using the query file's original stem (no extension) to avoid messy definition lines breaking the save path
    safe_name = query_fp.stem
    
    png_path = OUTPUT_DIR / f"mummer_{REF_NAME}_vs_{safe_name}.png"
    svg_path = OUTPUT_DIR / f"mummer_{REF_NAME}_vs_{safe_name}.svg"
    
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    
    print(f" -> Links plotted: {len(align_coords)}")
    print(f" -> Saved to: {png_path.name}")
    
    # CRITICAL: Close figure to free up memory during the loop
    plt.close(fig) 

print("\nBatch processing complete!")