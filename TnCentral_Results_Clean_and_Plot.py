# -*- coding: utf-8 -*-
"""
TnCentral BLASTn processing + alignment visualisation (zoom plotter)

1) Read per-contig TnCentral CSV summary tables -> merged_summary
2) Parse BLAST .fa.out -> coords (HSP-level hits)
3) De-duplicate exact positional duplicates -> coords_nodup
4) Merge coords_nodup with merged_summary -> combined
5) Subset to a target contig (e.g., plasmid1) -> hits_contig
6) Collapse overlapping hits (within contig+subject) by best per-hit bitscore -> hits_collapsed
7) Build gene_features from annotation CSV (flags filtered)
8) Reusable function make_zoom_plot() to plot any [zoom_start, zoom_end] window
   with gene map + hits colored by per-hit bitscore

"""
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from dna_features_viewer import GraphicFeature, GraphicRecord
from typing import Optional, List

# ============================================================
# ===================== EDITABLE PARAMETERS ==================
# ============================================================

# --- Input: TnCentral per-contig csv "alignments-table" files ---
FILES = {
    "chromosome": "chromosome_alignments-table.csv",
    "plasmid1": "plasmid1_alignments-table.csv",
}

# --- Input: BLAST text output (.fa.out) ---
FA_OUT_PATH = "JOB_C7F7BD.fa.out"

# --- Input: annotation CSV for the contig you want to plot ---
ANNOTATION_CSV = "D:/CRE_PROJECT_DATA_ANALYSIS/All_Anno/plasmid1/plasmid1_flagged_all.csv"

# --- Which contig are we plotting? ---
TARGET_CONTIG = "plasmid1"

# --- Reference/plot region for the full contig (used when building features) ---
CONTIG_START = 0
CONTIG_END   = 161_365   # plasmid length

# --- Flags to show on the gene map track ---
FLAGS_TO_KEEP = {f.lower() for f in {
    "IS", "IS_ir", "Promoter", "amr", "attC", "attI", "integrase", "integron", "virulence"
}}

# --- Colors for flags (keys must be lowercase) ---
FLAG_COLOR_MAP = {
    "is": "#009E73",
    "is_ir": "#009E73",
    "integron": "#E69F00",
    "integrase": "#D55E00",
    "attc": "#56B4E9",
    "atti": "#0072B2",
    "promoter": "#F0E442",
    "amr": "#EE6677",
    "virulence": "#CC79A7",
}
DEFAULT_FEATURE_COLOR = "#9b5de5"

# --- Overlap collapsing behavior ---
# If True: intervals that "touch" (start == previous end) are treated as overlap.
# If False: touching intervals are NOT overlaps (must actually overlap).
COUNT_TOUCHING_AS_OVERLAP = True

# --- Output file names ---
MERGED_SUMMARY_OUTCSV = "TnCentral_BLASTn_merged_alignments-table.csv"
HITS_COLLAPSED_OUTCSV = "plasmid_TnCentral_results_deduped_filtered.csv"

# --- styling ---
FIGSIZE = (22, 10)
CMAP_NAME = "cividis" # choose colour map
TICK_INTERVAL_DEFAULT = 5_000

# --- Zoom windows to plot:  (start, end, 'plot title', tick interval) ---
ZOOMS_TO_PLOT = [
     (38_000, 50_000, "Tn_alignments_gene_map_p1_zoom_38000_50000", 2_000),
     (48_000, 60_000, "Tn_alignments_gene_map_p1_zoom_48000_60000", 2_000),
    (0, 161_365, "Tn_alignments_gene_map_p1_full", 10_000)    # EDIT AS REQUIRED
]

# ============================================================
# ====================== HELPER FUNCTIONS ====================
# ============================================================

def read_and_clean_summary_tables(files: dict) -> pd.DataFrame:
    """Read TnCentral 'alignments-table.csv' per contig, clean, and return merged summary df."""
    dfs = []
    for contig, fp in files.items():
        d = pd.read_csv(fp, sep=",")
        d["contig"] = contig
        dfs.append(d)

    merged_summary = pd.concat(dfs, ignore_index=True)
    merged_summary.columns = merged_summary.columns.str.lower().str.strip()

    # Split description into name and accession (if '-' exists)
    if "description" in merged_summary.columns:
        merged_summary[["name", "tn_accession"]] = merged_summary["description"].astype(str).str.split("-", n=1, expand=True)
    else:
        merged_summary["name"] = pd.NA
        merged_summary["tn_accession"] = pd.NA

    merged_summary["source"] = "tncentral"

    rename_map = {
        "identities": "identity",
        "score": "bitscore",
        "description": "subject",
    }
    merged_summary = merged_summary.rename(columns=rename_map)

    return merged_summary


def parse_fa_out(fa_out_path: str) -> pd.DataFrame:
    """Parse BLAST .fa.out into an HSP-level dataframe."""
    records = []

    with open(fa_out_path) as f:
        curr_query = None
        curr_subject = None
        curr_subject_len = None
        curr_score = None
        curr_evalue = None
        curr_identity = None
        q_positions = []
        s_positions = []

        def flush_current():
            if curr_subject and curr_score is not None and q_positions and s_positions:
                qstart = min(q_positions)
                qend   = max(q_positions)
                sstart = min(s_positions)
                send   = max(s_positions)

                records.append({
                    "query": curr_query,
                    "subject": curr_subject,
                    "subject_length": curr_subject_len,
                    "bitscore": curr_score,
                    "evalue": curr_evalue,
                    "identity": curr_identity,
                    "qstart": qstart,
                    "qend": qend,
                    "sstart": sstart,
                    "send": send,
                })

            q_positions.clear()
            s_positions.clear()

        for raw in f:
            line = raw.rstrip("\n")

            if line.lstrip().startswith("Query="):
                flush_current()
                curr_query = line.split("=", 1)[1].strip()

            elif line.startswith(">"):
                flush_current()
                curr_subject = line[1:].strip()
                curr_subject_len = None

            elif line.strip().startswith("Length=") and curr_subject is not None:
                m = re.search(r"Length\s*=\s*(\d+)", line)
                if m:
                    curr_subject_len = int(m.group(1))

            elif line.strip().startswith("Score ="):
                flush_current()
                curr_identity = None

                m = re.search(
                    r"Score =\s*([\d\.]+)\s*bits.*Expect(?:\(\d+\))?\s*=\s*([\deE\-\+\.]+)",
                    line
                )
                if m:
                    curr_score = float(m.group(1))
                    curr_evalue = float(m.group(2))

            elif "Identities =" in line:
                m = re.search(r"Identities\s*=\s*(\d+)\s*/\s*(\d+)\s*\((\d+)%\)", line)
                if m:
                    curr_identity = float(m.group(3))
                else:
                    m2 = re.search(r"Identities\s*=\s*(\d+)\s*/\s*(\d+)", line)
                    if m2:
                        num = int(m2.group(1))
                        den = int(m2.group(2))
                        curr_identity = 100.0 * num / den if den > 0 else None

            elif line.strip().startswith("Query "):
                parts = line.split()
                if len(parts) >= 4:
                    q_positions.append(int(parts[1]))
                    q_positions.append(int(parts[3]))

            elif line.strip().startswith("Sbjct "):
                parts = line.split()
                if len(parts) >= 4:
                    s_positions.append(int(parts[1]))
                    s_positions.append(int(parts[3]))

        flush_current()

    coords = pd.DataFrame(records)
    return coords


def clean_coords(coords: pd.DataFrame) -> pd.DataFrame:
    """Clean parsed .fa.out coords and compute align_length + pct_coverage."""
    coords = coords.copy()

    coords["name"] = coords["subject"].astype(str).str.split("-", n=1).str[0]

    coords = coords.rename(columns={
        "query": "contig",
        "qstart": "start",
        "qend": "end",
        "sstart": "ref_start",
        "send": "ref_end",
    })

    coords["contig"] = coords["contig"].astype(str).str.split().str[0]

    # numeric coercion
    for c in ["start", "end", "ref_start", "ref_end", "bitscore", "identity", "evalue", "subject_length"]:
        if c in coords.columns:
            coords[c] = pd.to_numeric(coords[c], errors="coerce")

    coords = coords.sort_values(by=["contig", "subject", "start"], ascending=[True, True, True]).reset_index(drop=True)

    coords["align_length"] = coords["end"] - coords["start"] + 1
    coords["pct_coverage"] = (coords["align_length"] / coords["subject_length"]) * 100
    coords["pct_coverage"] = coords["pct_coverage"].round(2)

    return coords


def dedup_exact_positions(coords: pd.DataFrame) -> pd.DataFrame:
    """Drop exact positional duplicates per (contig, subject, start, end), keeping best bitscore/identity."""
    coords_sorted = coords.sort_values(
        by=["contig", "subject", "start", "end", "ref_start", "ref_end", "bitscore", "identity", "evalue"],
        ascending=[True, True, True, True, True, True, False, False, True],
    )

    coords_nodup = coords_sorted.drop_duplicates(
        subset=["contig", "subject", "start", "end"],
        keep="first"
    ).reset_index(drop=True)

    return coords_nodup


def collapse_overlaps_keep_best(df: pd.DataFrame, count_touching_as_overlap: bool = True) -> pd.DataFrame:
    """
    Collapse overlaps within each (contig, subject) into overlap-clusters.
    Keep the row with best per-hit bitscore in each cluster.
    """
    df = df.sort_values(["contig", "subject", "start", "end"]).copy()

    out = []
    for (_, _), grp in df.groupby(["contig", "subject"], sort=False):
        g = grp.sort_values(["start", "end"]).copy()

        running_end = g["end"].cummax().shift(fill_value=-float("inf"))
        if count_touching_as_overlap:
            new_cluster = g["start"] > running_end
        else:
            new_cluster = g["start"] >= running_end  # touching splits clusters

        g["_cluster"] = new_cluster.cumsum()

        best = (g.sort_values("bitscore", ascending=False)
                  .groupby("_cluster", as_index=False)
                  .head(1))

        out.append(best.drop(columns=["_cluster"]))

    return pd.concat(out, ignore_index=True).sort_values(["contig", "subject", "start", "end"])


def build_gene_features_from_annotation(
    annotation_csv: str,
    contig_start: int,
    contig_end: int,
    flags_to_keep: set,
    flag_color_map: dict,
    default_color: str,
) -> list:
    """Build dna_features_viewer GraphicFeature list from annotation CSV."""
    ann = pd.read_csv(annotation_csv, sep=",")

    ann["start"] = pd.to_numeric(ann["start"], errors="coerce")
    ann["end"] = pd.to_numeric(ann["end"], errors="coerce")

    gene_features = []

    for _, row in ann.iterrows():
        flag = row.get("flag")
        if pd.isna(flag):
            continue

        flag_lower = str(flag).strip().lower()
        if flag_lower not in flags_to_keep:
            continue

        # Clip to contig region and convert to local coords (0-based)
        f_start_abs = max(row["start"], contig_start)
        f_end_abs   = min(row["end"], contig_end)
        if pd.isna(f_start_abs) or pd.isna(f_end_abs) or f_end_abs < f_start_abs:
            continue

        f_start = int(f_start_abs - contig_start)
        f_end   = int(f_end_abs - contig_start)

        # Strand
        strand_val = str(row.get("strand", "+")).strip()
        strand = 1 if strand_val == "+" else -1 if strand_val == "-" else 0

        # Label                                 ###################################### Flip label from gene/name here
        gene_label = None
        gene = row.get("gene")
        if isinstance(gene, str) and gene.strip():
            gene_label = gene.strip()

        feature_color = flag_color_map.get(flag_lower, default_color)

        gene_features.append(
            GraphicFeature(
                start=f_start,
                end=f_end,
                strand=strand,
                color=feature_color,
                label=gene_label,
            )
        )

    return gene_features


def make_zoom_plot(
    gene_features: list,
    hits_df: pd.DataFrame,
    zoom_start: int,
    zoom_end: int,
    outfile_prefix: str,
    tick_interval: int = 5_000,
    figsize=(22, 10),
    cmap_name: str = "plasma",
    title: Optional[str] = None,
    show: bool = True,
    # NEW: keep y-axis order fixed across ALL plots, based on full-plasmid ordering
    subject_order_global: Optional[List[str]] = None,
):
    """
    2-panel plot (zoomed):
      - top: dna_features_viewer map (ONLY features overlapping zoom window)
      - bottom: BLAST hits (ONLY hits overlapping zoom window),
               colored by per-hit 'bitscore'

    Key behavior:
      - If subject_order_global is provided, y-axis order is consistent across plots:
        subjects shown in the zoom window are displayed in the SAME order as in subject_order_global.
      - Features/hits are clipped to the zoom window and shifted so x-axis is 0..(zoom_end-zoom_start)
      - X tick labels show absolute coordinates by default.
    """
    if zoom_end <= zoom_start:
        raise ValueError("zoom_end must be > zoom_start")

    hits = hits_df.copy()

    # Ensure numeric types
    for c in ["bitscore", "start", "end"]:
        if c not in hits.columns:
            raise KeyError(f"Expected column '{c}' in hits_df but it was not found.")
        hits[c] = pd.to_numeric(hits[c], errors="coerce")

    # ---- build zoomed gene record (clip + shift) ----
    gf_zoom = []
    for feat in gene_features:
        if feat.end < zoom_start or feat.start > zoom_end:
            continue

        f_start = max(feat.start, zoom_start) - zoom_start
        f_end   = min(feat.end,   zoom_end)   - zoom_start

        gf_zoom.append(
            GraphicFeature(
                start=f_start,
                end=f_end,
                strand=feat.strand,
                color=feat.color,
                label=feat.label
            )
        )

    region_len = zoom_end - zoom_start
    gene_record_zoom = GraphicRecord(sequence_length=region_len, features=gf_zoom)

    # ---- subset hits to zoom window, clip + shift ----
    hits_zoom = hits[
        (hits["end"] >= zoom_start) &
        (hits["start"] <= zoom_end)
    ].copy()

    # ---- figure ----
    fig, (ax_map, ax_hits) = plt.subplots(
        2, 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [1, 2]}
    )

    # top panel
    gene_record_zoom.plot(ax=ax_map, with_ruler=False, strand_in_label_threshold=7)
    ax_map.set_xticks([])
    ax_map.set_title(title if title else f"Zoom {zoom_start:,}–{zoom_end:,} bp")

    if hits_zoom.empty:
        # Still save a gene-only plot for this window
        ax_hits.set_xlim(0, region_len)
        ax_hits.set_yticks([])
        ax_hits.set_ylabel("transposon")
        ax_hits.set_xlabel("bp")

        fig.tight_layout()
        fig.savefig(f"{outfile_prefix}.pdf", dpi=600, bbox_inches="tight", pad_inches=0.1)
        fig.savefig(f"{outfile_prefix}.png", dpi=600, bbox_inches="tight", pad_inches=0.1, facecolor="white")

        if show:
            plt.show()
        plt.close(fig)
        return

    # Clip and shift hit coordinates to zoom window
    hits_zoom["start_zoom"] = hits_zoom["start"].clip(lower=zoom_start) - zoom_start
    hits_zoom["end_zoom"]   = hits_zoom["end"].clip(upper=zoom_end) - zoom_start
    hits_zoom["width_zoom"] = hits_zoom["end_zoom"] - hits_zoom["start_zoom"] + 1

    # ---- y-axis ordering: KEEP GLOBAL ORDER ----
    present = set(hits_zoom["name"].dropna().unique())

    if subject_order_global is not None:
        subject_order = [s for s in subject_order_global if s in present]
    else:
        # fallback if user forgets to pass global order
        if "bitscore_alignment" in hits_zoom.columns:
            hits_zoom["bitscore_alignment"] = pd.to_numeric(hits_zoom["bitscore_alignment"], errors="coerce")
            subject_order = (
                hits_zoom.groupby("name")["bitscore_alignment"].max()
                .sort_values(ascending=False)
                .index.tolist()
            )
        else:
            subject_order = (
                hits_zoom.groupby("name")["bitscore"].max()
                .sort_values(ascending=False)
                .index.tolist()
            )

    if not subject_order:
        # No named subjects survived (e.g., name missing). Still plot gene track.
        ax_hits.set_xlim(0, region_len)
        ax_hits.set_yticks([])
        ax_hits.set_ylabel("transposon")
        ax_hits.set_xlabel("bp")

        fig.tight_layout()
        fig.savefig(f"{outfile_prefix}.pdf", dpi=600, bbox_inches="tight", pad_inches=0.1)
        fig.savefig(f"{outfile_prefix}.png", dpi=600, bbox_inches="tight", pad_inches=0.1, facecolor="white")

        if show:
            plt.show()
        plt.close(fig)
        return

    subject_to_y = {s: i for i, s in enumerate(subject_order)}

    # ---- colormap (per-hit bitscore) ----
    vmin = hits_zoom["bitscore"].min()
    vmax = hits_zoom["bitscore"].max()
    norm = Normalize(vmin=vmin, vmax=vmax)

    cmap = getattr(cm, cmap_name) if hasattr(cm, cmap_name) else cm.plasma

    # ---- draw hits ----
    for _, row in hits_zoom.iterrows():
        name = row.get("name")
        if pd.isna(name) or name not in subject_to_y:
            continue

        y = subject_to_y[name]
        color = cmap(norm(row["bitscore"]))

        ax_hits.broken_barh(
            [(row["start_zoom"], row["width_zoom"])],
            (y - 0.4, 0.8),
            facecolors=color
        )

    # ---- axes formatting ----
    ax_hits.set_xlim(0, region_len)

    ticks = np.arange(0, region_len + tick_interval, tick_interval)
    ax_hits.set_xticks(ticks)
    ax_hits.set_xticklabels([f"{(t + zoom_start):,}" for t in ticks])  # absolute coords

    ax_hits.set_yticks(range(len(subject_order)))
    ax_hits.set_yticklabels(subject_order)
    ax_hits.invert_yaxis()
    ax_hits.set_xlabel("bp (absolute coordinates)")
    ax_hits.set_ylabel("transposon")

    # ---- layout + colorbar ----
    fig.tight_layout(rect=[0, 0, 0.94, 1])

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.94, 0.12, 0.015, 0.5])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("bitscore (per hit)")

    # ---- save ----
    fig.savefig(f"{outfile_prefix}.pdf", dpi=600, bbox_inches="tight", pad_inches=0.1)
    fig.savefig(f"{outfile_prefix}.png", dpi=600, bbox_inches="tight", pad_inches=0.1, facecolor="white")

    if show:
        plt.show()
    plt.close(fig)


# ============================================================
# ========================= MAIN PIPELINE =====================
# ============================================================

# 1) Read/clean TnCentral summary tables
merged_summary = read_and_clean_summary_tables(FILES)
merged_summary.to_csv(MERGED_SUMMARY_OUTCSV, index=False)

# 2) Parse + clean .fa.out coords (HSPs)
coords_raw = parse_fa_out(FA_OUT_PATH)
coords = clean_coords(coords_raw)

# 3) De-dupe exact positional duplicates
coords_nodup = dedup_exact_positions(coords)

# 4) Merge coords with summary table columns (adds *_alignment where overlaps)
combined = pd.merge(
    coords_nodup,
    merged_summary,
    how="left",
    on=["contig", "name"],
    suffixes=("", "_alignment"),
)

# 5) Subset to target contig
hits_contig = combined[combined["contig"] == TARGET_CONTIG].copy()

# 6) Collapse overlaps (within contig+subject) by best per-hit bitscore
hits_collapsed = collapse_overlaps_keep_best(
    hits_contig,
    count_touching_as_overlap=COUNT_TOUCHING_AS_OVERLAP
).reset_index(drop=True)

hits_collapsed["bitscore_alignment"] = pd.to_numeric(hits_collapsed["bitscore_alignment"], errors="coerce")

GLOBAL_SUBJECT_ORDER = (
    hits_collapsed
    .groupby("name")["bitscore_alignment"]
    .max()
    .sort_values(ascending=False)
    .index
    .tolist()
)

# Optional: sort for readability
sort_cols = ["bitscore_alignment", "name"] if "bitscore_alignment" in hits_collapsed.columns else ["bitscore", "name"]
hits_collapsed = hits_collapsed.sort_values(by=sort_cols, ascending=[False, True]).reset_index(drop=True)

# Save collapsed hits
hits_collapsed.to_csv(HITS_COLLAPSED_OUTCSV, index=False)

# 7) Build gene features from annotations
gene_features = build_gene_features_from_annotation(
    annotation_csv=ANNOTATION_CSV,
    contig_start=CONTIG_START,
    contig_end=CONTIG_END,
    flags_to_keep=FLAGS_TO_KEEP,
    flag_color_map=FLAG_COLOR_MAP,
    default_color=DEFAULT_FEATURE_COLOR,
)

# 8) Plot zooms

for zoom in ZOOMS_TO_PLOT:
    if len(zoom) == 3:
        zs, ze, prefix = zoom
        tick_interval = TICK_INTERVAL_DEFAULT
    elif len(zoom) == 4:
        zs, ze, prefix, tick_interval = zoom
    else:
        raise ValueError("Each zoom must be (start, end, prefix[, tick_interval])")

    make_zoom_plot(
        gene_features=gene_features,
        hits_df=hits_collapsed,
        zoom_start=zs,
        zoom_end=ze,
        outfile_prefix=prefix,
        tick_interval=tick_interval,
        figsize=FIGSIZE,
        cmap_name=CMAP_NAME,
        title=f"{TARGET_CONTIG} | TnCentral hits (zoom {zs:,}–{ze:,} bp)",
        subject_order_global=GLOBAL_SUBJECT_ORDER,
        show=True
    )
