# CS1.2_Project
Scripts used in the processing and analysis of outputs from various tools used in the characterisation of E. coli strain CS1.2

This is a brief summary of the files and their use.

# 1. Command Line Interface commands for utilised bioinformatics tools and supplementary information
- **Set up information**
- **Utilised tools and versions**
- **Data information**
- **Commands/ example commands covering:**
  - QC and filter data
  - Assembly
  - Pre-polish QC
  - Polish
  - Post-polish QC
  - Reorientation
  - Fasta manipulation
  - Create a new BAM file for reoriented
  - Final QC + CheckM
  - Plasmid detection/ characterisation
  - Strain characterisation
  - Bakta annotation and hit extraction
  - AMR and unification
  - Integrons
  - Local BlastN
  - MUMmer alignment
  - Clinker gene cluster alignment
  - pyCirclize MUMmer plot
  - Pling
  - MOBsuite cluster

# 2. Supplementary Command Line Scripts
- environment.yml : this is part of the Autocycler conda installation pipeline
- autocycler_bash.sh : file to run Autocycler
- extract_regions_from_genbank.py : script to extract specified loci by coordinates from genbank files
- convert_gbk_to_fasta_CLI.txt : Call to convert specified gbk files to fasta files
- mummer_plot_tune.py : mummer circos plot (synteny) for a single query/ref pair
- mummer_plot_batch.py : batch option for mummer circos plots 

# 3. Python Scripts: Cleaning and Processing Data
## A: VIRULENCE
- unify_virulence.py: clean and combine outputs from virulence gene detection tools 
## B: AMR
- clean_hAMR_output.py: clean up and fix issues with the hAMRonisation output
- AMR_exploratory.py: Exploratory analysis of AMR genes and distribution
- ast_heatmap.py: Takes AST data, resfinder phenotype prediction and genotype inferred resistance and creates a heatmap linking the observed resistance to detected resistance genes/determinants.  
## C: MGEs
- integron_finder_results.py: clean IntegronFinder output
- ISEScan_results.py: clean and filter the ISEScan output
- ISFinder_results.py: clean and filter ISFinder BLASTn results
- MobileElementFinder_results.py : clean and filter the MobileElementFinder results
- MGE_consolidation.py: combines the clean/filtered results from Integronfinder, ISEScan, ISFinder and MobileElementFinder
- TnCentral_Results_Clean_and_Plot.py : processes TnCentral Blast results, makes a plot of cleaned hits aligned to a gene map of the contig with hits coloured by bitscore.
## D: Toxin-antitoxin-Systems
- filter_clean_tadb_results.py: cleans and filters Toxin Antitoxin Database results.
## E: Processing of BlastN (plasmid1) results  
- GENE_extract_from_blast_out.py: Takes the genbank output of a blastN search, searches for accessions containing a specific “GENE”, extracts these accessions with their metadata and creates filtered genbank and fasta files for just these accessions.
- Process_local_blast_CS12ctxm15region_in_CTXM_plasmids.py: takes output of a local blastN of the CTXM15 region from the chromosome of CS1.2 (query) against a database of CTXM15-containing plasmid homologue hits. Adds in metadata. Flags those containing NDM5, filters out chromosomal hits and filters for high identity hits with available dates. 

# 4. Python Plotting Scripts
- LINEAR_dnafeatures_viewer.py: Creates linear gene plots for specified regions from the flagged annotation csv.
- plasmid_pycirclize_map.py: Pycirclize circos plot of plasmid based on flagged annotation CSV
- BRIG_multi_arrows.py: Create BRIG-like plots for plasmids
- Chromosome_GC_plot_final.py: Create a pyCirclize circos plot of the chromosome with GC content and AMR/ virulence feature flags
- pling_plasmid_distance_matrix.py: create a cluster heatmap from the Pling distance matrix

# Supplementary
- split_data_by_contig.py: Splits the Bakta output, AMR and virulence results by contig to facilitate annotation csv generation in Excel.
- Flagged annotation CSVs used for plotting in this project
- Fasta files for the final assembly (full), the chromosome and the IncF plasmid

