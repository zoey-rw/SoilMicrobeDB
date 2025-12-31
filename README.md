# Soil Microbe Database (SMDB) - a whole-genome reference database with over 50,000 genomes, focused on soil biodiversity 

Includes over 20,000 metagenome-assembled genomes (MAGs) from soils, and 1400 fungi from Joint Genome Institute's Mycocosm

For full details and genome sources, check out our preprint: https://www.biorxiv.org/content/10.1101/2025.03.21.644662v1

## Usage
Explore/download abundances of a subset of these genomes across NEON soil samples using the interactive portal, developed in partnership with the Boston University Software & Application Innovation Lab (SAIL): https://zoeywerbin.shinyapps.io/soil_microbe_db/

For local classification against metagenome samples, the database index files can be downloaded from an AWS bucket, along with Snakemake rules for each step of classification (Kraken2 for read classification, Architeuthis for kmer-filtering, and Bracken for abundance estimation). This is hosted through the AWS Open Data program - details will be added here as soon as available.

To add genomes to the database, you will need to download the database along with the Struo2 software, which contains a "db-update" workflow. The input file must have the filepaths of local genomes, along with their NCBI taxonomy, for compatibility with SoilGenomeDB.

## Running Analysis Scripts

### Dependencies

**R Packages:**
- `tidyverse`, `data.table` (core data manipulation)
- `pavian` (reading Kraken/Bracken files)
- Additional packages as needed by individual scripts

**Python Packages:**
- `biopython`

**External Tools:**
- Kraken2, Architeuthis, Bracken (classification pipeline)
- NCBI BLAST (optional, for validation scripts - can use web API)

### Script Organization

Scripts are organized in `scripts/summarize_outputs/` and numbered sequentially (01-19) to reflect execution order and dependencies.

**Foundation Scripts (01-13):** Data summarization and integration
- Process classification outputs
- Merge with external data (qPCR, PLFA, ITS)
- Calculate taxonomic summaries

**Database Comparison (14-15):** Compare classification results across databases
- `14_compare_db_classification_differences.r`: Compare GTDB vs SoilMicrobeDB classification metrics
- `15_compare_taxonomic_assignments.r`: Broad read-level comparison (identifies both PlusPF "Homo sapiens" and SoilMicrobeDB fungal → GTDB bacteria misclassifications, saves read IDs)

**BLAST Verification (16-18):** Validate taxonomic assignments 
**Note:** Scripts 16-18 require raw sample data, rather than summarized files
- `16_extract_reads_to_fasta.py`: Automatically extracts reads from both misclassification types to FASTA format
- `17_blast_reads_ncbi.py`: Automatically BLASTs all FASTA files against NCBI nt database (web API)
- `18_analyze_blast_results.r`: Automatically analyzes all BLAST results

### Running Scripts

Most scripts can be run directly:
```bash
Rscript scripts/summarize_outputs/14_compare_db_classification_differences.r
Rscript scripts/summarize_outputs/15_compare_taxonomic_assignments.r [sampleID]
python3 scripts/summarize_outputs/17_blast_reads_ncbi.py [sampleID] [sample_size]
```

## Credits & Acknowledgements

The National Ecological Observatory Network is a project solely funded by the National Science Foundation and managed under cooperative agreement by Battelle. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

## Disclaimer
Information and documents contained within this repository are available as-is. Codes or documents, or their use, may not be supported or maintained under any program or service and may not be compatible with data currently available from the NEON Data Portal.