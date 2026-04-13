# CLI Code to extract specific regions from genbank files (as specified by start/end coordinates)
# NOTE! Make sure this script is in the same directory as the genbank files to be processed!

# --> format for region extraction
#{
#        "input": "CS1_2_chromosome.gbff", 	# GenBank file name
#        "start": 258569, 			# region start coordinate
#        "end": 277391,				# region end coordinate
#        "output": "CS1_2_chr_region_NDM-5.gbk"  # Name for regional GenBank file 
#    }


from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path

# Base directory
base_dir = Path(".")

## sequences to extract --> EDIT AS REQUIRED ##

regions = [

    # =========================
    # CS1_2 chromosome regions
    # =========================
    {
        "input": "CS1_2_chromosome.gbff",
        "start": 258569,
        "end": 277931,
        "output": "CS1_2_chr_ctxmregion_258569_277931.gbk",
    },
    {
        "input": "CS1_2_chromosome.gbff",
        "start": 1885079,
        "end": 1925430,
        "output": "CS1_2_chr_HPI_1885079_1925430.gbk",
    },
    {
        "input": "CS1_2_chromosome.gbff",
        "start": 1858000,
        "end": 1883000,
        "output": "CS1_2_chr_cag_1858000_1883000.gbk",
    },

    # =========================
    # CS1_2 plasmid1 regions
    # =========================
    {
        "input": "CS1_2_plasmid1.gbff",
        "start": 31317,
        "end": 58619,
        "output": "CS1_2_plasmid1_NDMregion_31317_58619.gbk",
    },
    {
        "input": "CS1_2_plasmid1.gbff",
        "start": 143000,
        "end": 159000,
        "output": "CS1_2_plasmid1_iucsit_143000_159000.gbk",
    },
    {
        "input": "CS1_2_plasmid1.gbff",
        "start": 13000,
        "end": 17000,
        "output": "CS1_2_plasmid1_oxacat_13000_17000.gbk",
    },

    # =========================
    # Zurfluh strain 675SK2
    # =========================
    # chromosome
    {
        "input": "chr_bakta_675SK2.gbff",
        "start": 4567919,
        "end": 4587284,
        "output": "675SK2_chr_ctxmregion_4567919_4587284.gbk",
    },

    # plasmid
    {
        "input": "plasmid_bakta_675SK2.gbff",
        "start": 54940,
        "end": 78314,
        "output": "675SK2_plasmid_ndmregion_54940_78314.gbk",
    },
]

# Implementation #
for r in regions:
    record = SeqIO.read(base_dir / r["input"], "genbank")

    # Convert to 0-based indexing
    start = r["start"] - 1
    end = r["end"]

    sub_seq = record.seq[start:end]

    new_features = []
    for f in record.features:
        if f.location.start >= start and f.location.end <= end:
            new_loc = FeatureLocation(
                f.location.start - start,
                f.location.end - start,
                strand=f.location.strand
            )
            new_feat = SeqFeature(
                location=new_loc,
                type=f.type,
                qualifiers=f.qualifiers
            )
            new_features.append(new_feat)

    record.seq = sub_seq
    record.features = new_features
    record.annotations["molecule_type"] = "DNA"
    record.id = record.id + "_region"
    record.name = record.name + "_region"
    record.description = f"{r['input']}:{r['start']}..{r['end']}"

    SeqIO.write(record, base_dir / r["output"], "genbank")

print("Region extraction complete.")
