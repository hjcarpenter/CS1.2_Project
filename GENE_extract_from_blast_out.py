# -*- coding: utf-8 -*-
"""
# 1: Takes the full genbank output of plasmid1 fasta blastn 
(megablastN, target 5000, default E)
Search the multigenbank file for accessions containing GENE
Extract accessions and metadata

# 2:
Make a filtered gbk with the GENE set
Make a fasta from the filtered gbk. 

** Run for CTX-M-15 and NDM-5
***Remember to EDIT out file names to match GENE***

Created on Sat Jan 17 10:23:32 2026

@author: Heath
"""

from Bio import SeqIO
import pandas as pd
import re

def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r"(19|20)\d{2}", str(date_str))
    return int(match.group()) if match else None

gb_file = r"D:\CRE_PROJECT_DATA_ANALYSIS\BLAST\sequence.gb" # multigenbank file

GENE = "NDM-5"  # <----- set GENE here!!!!

rows = []

for record in SeqIO.parse(gb_file, "genbank"):
    accession = record.id
    seq_len = len(record)
    organism = record.annotations.get("organism", None)

    # Default metadata
    collection_date = None
    country = None
    isolation_source = None
    host = None
    strain = None

    # Pull metadata from source feature
    for feature in record.features:
        if feature.type == "source":
            q = feature.qualifiers
            collection_date = q.get("collection_date", [None])[0]
            country = q.get("geo_loc_name", [None])[0]
            isolation_source = q.get("isolation_source", [None])[0]
            host = q.get("host", [None])[0]
            strain = q.get("strain", [None])[0]
            break

    # Look for GENE
    found = False
    for feature in record.features:
        if feature.type in ["CDS", "gene"]:
            text = " ".join(
                feature.qualifiers.get("gene", []) +
                feature.qualifiers.get("product", [])
            ).upper()

            if GENE in text:
                found = True
                break

    if found:
        rows.append({
            "accession": accession,
            "length": seq_len,
            "organism": organism,
            "strain": strain,
            "host": host,
            "collection_date": collection_date,
            "geo_loc_name": country,
            "isolation_source": isolation_source
        })

df = pd.DataFrame(rows)

print(df)
# get year
df["collection_year"] = df["collection_date"].apply(extract_year)
# save
df.to_csv(f"{GENE}_metadata_table.csv", index=False)

# ===================== 2 ==========================================

# ---- MAKE A FILTERED MULTIGENBANK -----------

# filter the gbk file for just the plasmids with GENE
gb_in = gb_file
gb_out = r"D:\CRE_PROJECT_DATA_ANALYSIS\BLAST\sequence_NDM5_only.gb"


# Set of accessions with GENE (versioned)
keep = set(df["accession"].dropna().astype(str))

count = 0
with open(gb_out, "w") as out_handle:
    for record in SeqIO.parse(gb_in, "genbank"):
        if record.id in keep:
            SeqIO.write(record, out_handle, "genbank")
            count += 1

print(f"Wrote {count} {GENE}-positive plasmids to:\n{gb_out}")

# check
# written = [r.id for r in SeqIO.parse(gb_out, "genbank")]
# set(written) == keep

# ---- MAKE A FASTA FROM THE GBK --------

# Get fasta files for a blast database
fasta_out = r"D:\CRE_PROJECT_DATA_ANALYSIS\BLAST\sequence_NDM5_only.fasta"
SeqIO.convert(gb_out, "genbank", fasta_out, "fasta")
print("Wrote:", fasta_out)

