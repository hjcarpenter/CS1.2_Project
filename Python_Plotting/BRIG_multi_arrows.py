# -*- coding: utf-8 -*-
"""
Generates a BRIG-like plot using pCS1.2IncF-NDM as a reference. 
Modified to dynamically load query fastas and parse headers.
"""

from pathlib import Path
import pandas as pd
import re

from pycirclize import Circos
from pygenomeviz.parser import Fasta
from pygenomeviz.utils import interpolate_color
from pygenomeviz.align import AlignCoord, Blast
from matplotlib.patches import Patch

# --- DIRECTORIES & FILES ---
# Path to the reference FASTA
fpREF = 'D:/CRE_PROJECT_DATA_ANALYSIS/Plasmid1_Blastn/plasmid_fastas/plasmid1.fasta'

# Path to the annotation CSV
annot_csv = 'D:/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all_fullmapUTF8.csv'

# Directory containing ALL the query fastas you want to map against the reference
QUERY_DIR = 'D:/CRE_PROJECT_DATA_ANALYSIS/Plasmid1_Blastn/plasmid_fastas/p_for_analysis_multi.fasta.split/BRIG/subcommunitiesABCD' 

REF_NAME = "pCS1.2IncF-NDM"  # for the center

# --- HELPER FUNCTION: PARSE FASTA HEADERS ---
def get_fasta_info(fasta_path):
    """Extracts accession and plasmid name from a FASTA header."""
    with open(fasta_path, 'r', encoding='utf-8') as f:
        header = f.readline().strip()
    
    # 1. Accession: everything after '>' up to the first space
    acc_match = re.search(r'^>?([^\s]+)', header)
    accession = acc_match.group(1) if acc_match else Path(fasta_path).stem
    
    # 2. Plasmid name: everything after 'plasmid ' up to the first comma or space
    plas_match = re.search(r'plasmid\s+([^,\s]+)', header, re.IGNORECASE)
    plasmid_name = plas_match.group(1) if plas_match else "Unknown_Plasmid"
    
    return plasmid_name, accession

# --- DYNAMICALLY LOAD COMPARISONS ---
comparisons = []
LEGEND_TEXT = {}

# Iterate through all fasta files in the QUERY_DIR
for fp in Path(QUERY_DIR).glob("*.fa*"): # Catches .fasta, .fa, .fna
    p_name, acc = get_fasta_info(fp)
    
    # Create a unique short ID for the script's internal dictionary
    disp_name = f"{p_name}_{acc}" 
    
    # Create the nice, readable format for the plot legend
    legend_label = f"{p_name} [{acc}]"
    
    comparisons.append((disp_name, Fasta(fp)))
    LEGEND_TEXT[disp_name] = legend_label

# --- EXTENDED COLORBLIND-SAFE PALETTE ---
# Combined Okabe-Ito and Paul Tol's palettes for distinct colors
CB_SAFE_PALETTE = [
    "#E69F00", # Orange (Okabe-Ito)
    "#56B4E9", # Sky Blue (Okabe-Ito)
    "#009E73", # Bluish Green (Okabe-Ito)
    "#D55E00", # Vermillion (Okabe-Ito)
    "#F0E442", # Yellow (Okabe-Ito)
    "#0072B2", # Blue (Okabe-Ito)
    "#CC79A7", # Reddish Purple (Okabe-Ito)
    "#332288", # Indigo (Tol)
    "#88CCEE", # Cyan (Tol)
    "#44AA99", # Teal (Tol)
    "#117733", # Green (Tol)
    "#999933", # Olive (Tol)
    "#DDCC77", # Sand (Tol)
    "#CC6677", # Rose (Tol)
    "#882255", # Wine (Tol)
    "#AA4499",  # Purple (Tol)
    "#888888"  # grey
]

# --- Plot + BLAST settings ---
QUERY_TRACK_SIZE = 7        # FOR LOTS of rings this needs to be smaller! 
MIN_IDENTITY = 95            # min blast alignment threshold to display (%)
TICKS_INTERVAL = 20_000

# Legend fade identities: from min identity to 100.. edit as needed
LEGEND_IDENTITIES = [100, 97, 95] 

# Label/tick settings
LABEL_TRACK = (97.0, 100.0)   
TICK_W_BP = 120               
TEXT_R = 99.8                 

# BRIG outer ring position
BRIG_START_R = 95.0           
AXIS_PAD = 0.7                

# --- Load Target FASTA & Annotations ---
if not Path(fpREF).exists():
    raise FileNotFoundError(f"Reference FASTA not found: {fpREF}")
if not Path(annot_csv).exists():
    raise FileNotFoundError(f"Annotation CSV not found: {annot_csv}")

target_fasta = Fasta(fpREF)
df = pd.read_csv(annot_csv, encoding='utf-8-sig')

# Filter AMR + virulence rows in annotation df by flag for labels
flag_s = df["flag"].astype(str)
mask_amr = flag_s.str.contains("amr", case=False, na=False)
mask_vir = flag_s.str.contains("vir", case=False, na=False)  

feat = df[(mask_amr | mask_vir) ].copy()
feat["is_amr"] = feat["flag"].astype(str).str.contains("amr", case=False, na=False)

# Pick label text with fallbacks
def pick_label(row):
    for col in ["name", "name_auto", "gene", "product"]:
        v = row.get(col)
        if pd.notna(v) and str(v).strip():
            return str(v)
    return "feature"

feat["label"] = feat.apply(pick_label, axis=1)

# --- Initialise circos ---
seqid2size = target_fasta.get_seqid2size()
circos = Circos(sectors=seqid2size, space=0 if len(seqid2size) == 1 else 2)
genome_kb = target_fasta.full_genome_length / 1000
circos.text(f"{REF_NAME}\n({genome_kb:.1f} kb)", size=13)

# Ensure CSV seqids match FASTA sector names
ref_seqids = set(seqid2size.keys())
feat = feat[feat["sequence_id"].isin(ref_seqids)].copy()

# --- Add Arrows and Labels track ---
ARROW_TRACK = (95.0, 99.0) 
ARROW_EC = "black"         # Optional: add a border to arrows for visibility

for sector in circos.sectors:
    sdf = feat[feat["sequence_id"] == sector.name].sort_values("start")
    if sdf.empty:
        continue
        
    # Create the track for the arrows
    track = sector.add_track(ARROW_TRACK, r_pad_ratio=0.05)

    for _, row in sdf.iterrows():
        start, end = int(row["start"]), int(row["end"])
        mid = (start + end) // 2
        strand = str(row.get("strand", "+"))
        
        # Color logic
        color = "#D55E00" if bool(row["is_amr"]) else "#0072B2"

        # 1. Draw the Arrow
        forward = not strand.startswith("-")
        if forward:
            track.arrow(start, end, fc=color, ec=ARROW_EC)
        else:
            track.arrow(end, start, fc=color, ec=ARROW_EC)

        # 2. Draw the Line & Text Label (pointing outward from the arrow)
        sector.line(r=(99, 102), start=mid, end=mid, color="grey", lw=0.5)
        sector.text(
            row["label"],
            x=mid,
            r=102.5,               
            size=10,                  
            orientation="vertical",
            ha="center",         
            va="bottom"          
        )

# --- BRIG rings (inside) ---
min_r_pos = BRIG_START_R
comp_name2color = {}

for i, (disp_name, comp_fasta) in enumerate(comparisons):
    align_coords = Blast([target_fasta, comp_fasta]).run()
    align_coords = AlignCoord.filter(align_coords, identity_thr=MIN_IDENTITY)

    # Pick color from extended palette (loops back if > 16 comparisons)
    color = CB_SAFE_PALETTE[i % len(CB_SAFE_PALETTE)]
    comp_name2color[disp_name] = color
    
    min_r_pos -= QUERY_TRACK_SIZE
    for sector in circos.sectors:
        sector.add_track((min_r_pos, min_r_pos + QUERY_TRACK_SIZE), r_pad_ratio=0.1)

    for ac in align_coords:
        track = circos.get_sector(ac.query_name).tracks[-1]
        rect_color = interpolate_color(color, v=ac.identity, vmin=MIN_IDENTITY)
        track.rect(ac.query_start, ac.query_end, color=rect_color)

# --- Reference axis + ticks ---
for sector in circos.sectors:
    track = sector.add_track((min_r_pos - AXIS_PAD, min_r_pos - 0.2))
    track.axis(fc="black")
    if sector.size >= TICKS_INTERVAL:
        track.xticks_by_interval(
            TICKS_INTERVAL,
            outer=False,
            label_formatter=lambda v: f"{int(v/1000)} kb",
            label_orientation="vertical",
        )

fig = circos.plotfig(figsize=(15,15))

# --- Build a BRIG-style Gradient Legend ---
legend_handles = []

for disp_name, color in comp_name2color.items():
    # Strain Name (Title)
    full_name = LEGEND_TEXT.get(disp_name, disp_name)
    legend_handles.append(Patch(facecolor='none', edgecolor='none', label=full_name))
    
    # Gradient Boxes
    for ident in LEGEND_IDENTITIES:
        shade = interpolate_color(color, v=ident, vmin=MIN_IDENTITY)
        legend_handles.append(Patch(facecolor=shade, edgecolor='none', label=f"    {ident}% identity"))

# Draw the legend
_ = fig.legend(
    handles=legend_handles, 
    loc="upper left",
    bbox_to_anchor=(-0.25, 1.0),   
    fontsize=9,                    
    frameon=False,                 
    labelspacing=0.4               
)

# save: 2 formats
fig.savefig(f"ABCDbrig_ARROWS_blast_thresh_{MIN_IDENTITY}.svg", bbox_inches="tight")
fig.savefig(f"ABCDbrig_ARROWS_blast_thresh_{MIN_IDENTITY}.pdf", bbox_inches="tight")
