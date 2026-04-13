#!/usr/bin/env bash#
# 
# Wrapper for a fully-automated Autocycler assembly with output dir based on input FASTQ location/name.
# Usage:
# autocycler_bash.sh <read_fastq> <threads> <jobs> [read_type]
# Copyright 2025 Ryan Wick
# Licensed under GPL-3.0
#
# Modifications made for this project (2025): 
# Selected assemblers only (raven miniasm flye plassembler canu), dynamic output dir (next to reads), safer GNU parallel timeout, robustness tweaks,
# and automatic FASTA export from final GFAs (keeping default headers).
# Consensus weight of 2 given to Flye and Canu (better base-level accuracy than the other assemblers)
# Plassembler cluster weight 3 given to Plassembler (Plassembler is specialised to pick up small plasmids - often missed by other assemblers -
# so high trust given to structural topology from this assembler).
#
# usage:
# ./autocycler_bash.sh <read_fastq> <threads> <jobs> [read_type]

set -e

# Get arguments.
reads=$1                 # input reads FASTQ
threads=$2               # threads per job
jobs=$3                  # number of simultaneous jobs
read_type=${4:-ont_r10}  # read type (default = ont_r10)

# Input assembly jobs that exceed this time limit will be killed (GNU parallel expects seconds)
max_time="28800"         # 8 hours in seconds

# Validate input parameters
if [[ -z "$reads" || -z "$threads" || -z "$jobs" ]]; then
    echo "Usage: $0 <read_fastq> <threads> <jobs> [read_type]" 1>&2
    exit 1
fi
if [[ ! -f "$reads" ]]; then
    echo "Error: Input file '$reads' does not exist." 1>&2
    exit 1
fi
if (( threads > 128 )); then threads=128; fi  # Flye won't work with >128 threads
case $read_type in
    ont_r9|ont_r10|pacbio_clr|pacbio_hifi) ;;
    *) echo "Error: read_type must be ont_r9, ont_r10, pacbio_clr or pacbio_hifi" 1>&2; exit 1 ;;
esac

# Figure out output directory based on input FASTQ, in the SAME directory as the FASTQ
input_dir=$(dirname "$reads")
input_base=$(basename "$reads")

# Strip common FASTQ extensions to form a clean base
base_noext="$input_base"
case "$base_noext" in
    *.fastq.gz) base_noext="${base_noext%.fastq.gz}" ;;
    *.fq.gz)    base_noext="${base_noext%.fq.gz}" ;;
    *.fastq)    base_noext="${base_noext%.fastq}" ;;
    *.fq)       base_noext="${base_noext%.fq}" ;;
esac

out_dir="${input_dir}/${base_noext}_autocycler"
mkdir -p "$out_dir"
cd "$out_dir"

# Compute genome size
genome_size=$(autocycler helper genome_size --reads "$reads" --threads "$threads")

# Step 1: subsample the long-read set into multiple files
autocycler subsample \
  --reads "$reads" \
  --out_dir subsampled_reads \
  --genome_size "$genome_size" 2>> autocycler.stderr

# Step 2: assemble each subsampled file
mkdir -p assemblies
rm -f assemblies/jobs.txt

for assembler in raven miniasm flye plassembler canu; do
    for i in 01 02 03 04; do
        echo "autocycler helper $assembler --reads subsampled_reads/sample_$i.fastq --out_prefix assemblies/${assembler}_$i --threads $threads --genome_size $genome_size --read_type $read_type --min_depth_rel 0.1" >> assemblies/jobs.txt
    done
done

set +e
nice -n 19 parallel --jobs "$jobs" --joblog assemblies/joblog.tsv --results assemblies/logs --timeout "$max_time" < assemblies/jobs.txt
set -e

# Give circular contigs from Plassembler extra clustering weight (for better structural resolution of small plasmids)
shopt -s nullglob
for f in assemblies/plassembler*.fasta; do
    sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "$f"
done

# Give contigs from Canu and Flye extra consensus weight (highest sequence fidelity)
for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
done
shopt -u nullglob

# Remove the subsampled reads to save space (won't error if none)
rm -f subsampled_reads/*.fastq

# Step 3: compress the input assemblies into a unitig graph
autocycler compress -i assemblies -a autocycler_out 2>> autocycler.stderr

# Step 4: cluster the input contigs into putative genomic sequences
autocycler cluster -a autocycler_out 2>> autocycler.stderr

# Steps 5 and 6: trim and resolve each QC-pass cluster
for c in autocycler_out/clustering/qc_pass/cluster_*; do
    [[ -d "$c" ]] || continue
    autocycler trim -c "$c" 2>> autocycler.stderr
    autocycler resolve -c "$c" 2>> autocycler.stderr
done

# Step 7: combine resolved clusters into a final assembly
autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa 2>> autocycler.stderr

###############################################################################
# FASTA EXPORT (keeps default Autocycler/contig headers)
###############################################################################

# Function: convert GFA (segments) to FASTA using awk (no external deps)
gfa_to_fasta () {
    # Prints only segments with real sequences (S lines where field 3 != "*")
    awk 'BEGIN{FS="\t"} $1=="S" && $3!="*" {print ">"$2; print $3}' "$1"
}

echo "Exporting FASTA from GFA outputs..."

# 1) Per-cluster finals concatenated into a single file at the run root
consensus_fasta="final_consensus.fasta"
: > "$consensus_fasta"
found_any=0
shopt -s nullglob
for gfa in autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa; do
    if [[ -f "$gfa" ]]; then
        gfa_to_fasta "$gfa" >> "$consensus_fasta"
        found_any=1
    fi
done
shopt -u nullglob
if [[ "$found_any" -eq 1 ]]; then
    echo " - Wrote consolidated per-cluster FASTA: $consensus_fasta"
else
    echo " - Warning: no per-cluster 5_final.gfa files found; $consensus_fasta is empty."
fi

# 2) If a combined GFA exists, also emit a paired combined FASTA next to it
combined_gfa="autocycler_out/combined/final_assembly.gfa"
combined_fa="autocycler_out/combined/final_assembly.fasta"
if [[ -f "$combined_gfa" ]]; then
    gfa_to_fasta "$combined_gfa" > "$combined_fa"
    echo " - Wrote combined FASTA: $combined_fa"
else
    echo " - Combined GFA not found; skipping $combined_fa"
fi

echo "Run complete."
echo "Results are in: $out_dir"
echo "Key outputs:"
echo " - Autocycler project: $out_dir/autocycler_out/"
echo " - Final per-cluster GFA: $out_dir/autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa"
echo " - Consolidated FASTA of per-cluster finals: $out_dir/$consensus_fasta"
echo " - Combined GFA/FASTA (if produced): $out_dir/$combined_gfa and $out_dir/$combined_fa"
