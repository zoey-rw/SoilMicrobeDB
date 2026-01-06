# Document organization: AWS S3 bucket file descriptions 

This repository includes development and analysis code for the SoilMicrobeDB, as well as a subdirectory with information on how to use the database for your own analysis via an AWS bucket.

The tutorial for using the SoilMicrobeDB is available here: https://github.com/zoey-rw/SoilMicrobeDB/blob/master/AWS_conda_environment_example/tutorial.md

## S3 Bucket Structure

The S3 bucket `s3://kraken2-soil-microbe-database/` contains the following files and directories:

### Database Files (in `kraken/` subdirectory)

All Kraken2 and Bracken database files are located in the `kraken/` subdirectory. When downloading, these should be placed in a local `kraken2/` directory (note: the local directory name is `kraken2`, while the S3 subdirectory is `kraken`).

**Essential files for Kraken2 classification (Step 2):**
These files are required for basic read classification. Total size: ~289 GB.
- `kraken/taxo.k2d` - Taxonomy index file
- `kraken/seqid2taxid.map` - Sequence ID to taxonomy ID mapping file
- `kraken/opts.k2d` - Database options file
- `kraken/hash.k2d` - Hash table file (largest file, ~280+ GB)
- `kraken/taxonomy/` - Directory containing NCBI taxonomy files:
  - `taxonomy/names.dmp` - Taxonomy names
  - `taxonomy/nodes.dmp` - Taxonomy tree structure
  - Additional taxonomy files as needed

**Additional files for Bracken abundance estimation (Step 4):**
These files are only needed if you plan to estimate relative abundances using Bracken. They are included in the full database download (~400 MB additional).
- `kraken/database150mers.kraken` - Bracken database for 150bp reads
- `kraken/database150mers.kmer_distrib` - K-mer distribution file for 150bp reads
- `kraken/database.kraken` - Bracken database file

**Note:** If you only need classification (Step 2) without abundance estimation, you can skip downloading the Bracken files to save ~400 MB. However, the full download includes all files.

### Reference Data Files (at bucket root)

**Genome list:**
- `soil_microbe_db_genome_table.csv` - Complete genome list used as input for the Struo2 database creation pipeline. This TSV file contains columns for genome name, taxonomy ID, and filepath.

**Test sample files:**
- `test_sample_R1.fastq.gz` - Forward reads for test sample
- `test_sample_R2.fastq.gz` - Reverse reads for test sample

**Analysis results (optional reference):**
- `soil_microbe_db_filtered_genus_merged_lineage.csv` - Example results file created from mapping the NEON soil sample collection against the SoilMicrobeDB, with all output files merged at the genus taxonomic rank and lineage information added.

## Download Instructions

### Full Database Download (Recommended)

To download the complete database including all Kraken2 and Bracken files (~289 GB):

```bash
# Create local directory
mkdir -p /databases/SoilMicrobeDB/kraken2

# Download all database files from the kraken/ subdirectory
aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/ /databases/SoilMicrobeDB/kraken2/
```

### Minimal Download (Classification Only)

If you only need Kraken2 classification (Step 2) and will skip Bracken abundance estimation (Step 4), you can download only the essential files:

```bash
# Create local directory
mkdir -p /databases/SoilMicrobeDB/kraken2

# Download only essential Kraken2 files
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/kraken/taxo.k2d /databases/SoilMicrobeDB/kraken2/
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/kraken/seqid2taxid.map /databases/SoilMicrobeDB/kraken2/
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/kraken/opts.k2d /databases/SoilMicrobeDB/kraken2/
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/kraken/hash.k2d /databases/SoilMicrobeDB/kraken2/

# Download taxonomy directory
aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/taxonomy/ /databases/SoilMicrobeDB/kraken2/taxonomy/
```

**Storage Requirements:**
- Full database (Kraken2 + Bracken): ~289 GB
- Kraken2 only (classification without abundance estimation): ~289 GB (Bracken files are small relative to hash.k2d)
- RAM requirements: ~300 GB for optimal performance (can use `--memory-mapping` flag for lower RAM usage)

### Test Sample Files

To download test sample files:

```bash
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R1.fastq.gz ./
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R2.fastq.gz ./
```

For complete usage instructions, see the [tutorial](https://github.com/zoey-rw/SoilMicrobeDB/blob/master/AWS_conda_environment_example/tutorial.md).
