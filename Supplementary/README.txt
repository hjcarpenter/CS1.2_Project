Assembly and Annotation of E. coli strain CS1.2 SUPPLEMENTARY FILES
-------------------------------------------------------------------

ANNOTATION_FOR_PLOTS: 
Flagged annotation .csv are the basis for pyCirclize, DNAFeaturesViewer and Clinker visualisations generated in this project. They consist of the processed and combined outputs of the following tools:

Bakta 1.11.0 db 6.0 full
RGI 6.0.5, CARD 4.0.1: from https://card.mcmaster.ca/analyze/rgi (hits: Perfect and Strict) 
starAMR 0.11.1+galaxy0; starAMR databases with ResFinder: 2.4.0_d1e607b_2024-08-06, PointFinder: 4.1.1_694919f_2024-08-08, PlasmidFinder: 3e77502_2024-03-07: from https://galaxy-main.usegalaxy.org/
AMRFinderPlus 4.0.23, database 2025-07-16.1
hAMRonization 1.1.9
VirulenceFinder 2.0.5, database 2022-12-02: from https://cge.food.dtu.dk/services/VirulenceFinder/ (species: Escherichia coli)
TAfinder 2.0, database 3.0: from https://bioinfo-mml.sjtu.edu.cn/TADB3/
MobileElementFinder 1.0.3, database 1.0.2: from https://cge.food.dtu.dk/services/MobileElementFinder/
IntegronFinder 2.0.6
ISEScan Galaxy Version 1.7.3+galaxy0: from https://galaxy-main.usegalaxy.org/
TnCentral BLAST from: https://tncentral.ncc.unesp.br/blast/

Manual edits were largely limited to antimicrobial resistance-associated regions and were based on results from these databases:
ISFinder database 2025-09-30: from isfinder.biotoul.fr
TnCentral database: from https://tncentral.ncc.unesp.br/
NCBI GenBank: from https://www-ncbi-nlm-nih-gov/genbank/
UniProt: from uniport.org

### contents ###
chromosome_flagged_all.csv : Combined Flagged annotations for the entire CS1.2 chromosome (contig01)

chromosome_CTXM15_region.csv: Flagged annotation with manual edit for the chromosomal (contig01) blaCTX-M-15 insertion region.

pCS1.2IncF-NDM_flagged_fullmapUTF8.csv : Combined Flagged annotations with manual edits for plasmid pCS1.2IncF-NDM (contig02)

-----------------------

FASTA FILES:
final_assembly.fasta: Complete Autocycler assembly post-polish and reorientation (contigs01-05)
chromosome.fasta: individual fasta file for chromosomal contig01
plasmid1.fasta: individual fasta file for pCS1.2IncF-NDM (contig02)

