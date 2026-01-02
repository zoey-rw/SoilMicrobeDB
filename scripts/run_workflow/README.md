# Workflow Scripts - Classification Pipeline

## Script Overview

1. **01_run_kraken.sh** - Runs Kraken2 classification
2. **02_run_architeuthis.sh** - Runs Architeuthis filtering and Bracken abundance estimation
3. **03_add_lineage.sh** - Merges Bracken outputs and adds lineage information
4. **04_reshape_score_reads.r** - Extracts scoring information from `_scores.output` files
5. **06_calculate_classification_pct_by_rank.r** - Calculates classification percentages by taxonomic rank
6. **05_cleanup_intermediate_files.sh** - Cleanup that checks dependencies before deleting files

## Workflow Execution Order

### Step 1: Kraken2 Classification
**Script**: `01_run_kraken.sh`
- **Input**: FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)
- **Output**: 
  - `*_kraken.kreport` files → `data/classification/01_kraken_output/` (preserved)
  - `*_kraken.output` files → `data/classification/01_kraken_output/` (preserved)
- **Usage**: Edit the `DBNAME` variable at the top of the script to process different databases. Run once per database.

### Step 2: Architeuthis Filtering and Bracken
**Script**: `02_run_architeuthis.sh`
- **Input**: `*_kraken.output` files from Step 1
- **Output**: 
  - `*_scores.output` files → `data/classification/02_bracken_output/` (preserved until extraction)
  - `*_filtered_kraken.kreport` files → `data/classification/02_bracken_output/` (preserved)
  - `*.b2` files (species/genus/domain/phylum) → `data/classification/02_bracken_output/` (preserved until merged)
- **Process**: 
  - Runs Architeuthis mapping score
  - Converts filtered output to kreport format
  - Runs Bracken at species, genus, domain, and phylum levels
- **Usage**: Edit the `DBNAME` variable at the top of the script (lines 23-26) to process different databases. Run once per database.

### Step 3: Merge and Add Lineage
**Script**: `03_add_lineage.sh`
- **Input**: 
  - `*.b2` files (species and phylum levels) from Step 2
- **Output**: 
  - `*_species_merged_lineage.csv` files → `data/classification/taxonomic_rank_summaries/species/` (preserved - final output)
  - `*_filtered_phylum_merged_lineage.csv` files → `data/classification/taxonomic_rank_summaries/phylum/` (preserved - final output)
- **Process**: 
  - Merges all samples using `architeuthis merge`
  - Adds lineage using `architeuthis lineage`
- **Databases**: Processes all databases automatically (soil_microbe_db, gtdb_207, gtdb_207_unfiltered, pluspf)
- **Usage**: No changes needed - processes all databases automatically

### Step 4: Extract Scoring Information
**Script**: `04_reshape_score_reads.r`
- **Input**: 
  - `*_scores.output` files from Step 2 
  - `data/classification/analysis_files/seq_depth_df.rds` (from `01_calculate_sequencing_depth.r`)
- **Output**: 
  - `data/classification/analysis_files/filter_results_summary.csv` (preserved, used by visualization scripts)
  - `data/classification/analysis_files/filter_results_processed_files.txt` (log of processed files)
- **Usage**: `Rscript scripts/run_workflow/04_reshape_score_reads.r`
- **Note**: Run this before cleanup to extract all scoring information

### Step 5: Calculate Classification Percentages by Rank
**Script**: `06_calculate_classification_pct_by_rank.r`
- **Input**: 
  - `*_kraken.kreport` files from Step 1
  - `*_filtered_kraken.kreport` files from Step 2
  - `data/classification/analysis_files/seq_depth_df.rds` (from `01_calculate_sequencing_depth.r`)
- **Output**: 
  - `data/classification/analysis_files/classification_pct_by_rank.csv` (preserved, used by visualization scripts)
  - `data/classification/analysis_files/classification_pct_by_rank_per_sample.csv` (preserved)
- **Process**: 
  - Calculates percentage of reads classified at each taxonomic rank (phylum → class → order → family → genus → species → strain)
  - Calculates for both before and after Architeuthis filtering
  - Calculates for all domains and fungi-specific
- **Usage**: `Rscript scripts/run_workflow/06_calculate_classification_pct_by_rank.r`
- **Note**: Can run after Step 2 (only needs kreport files). Run before cleanup if you want to enable kreport file deletion.

### Step 6: Cleanup (Optional)
**Script**: `05_cleanup_intermediate_files.sh`
- **Purpose**: Smart cleanup that checks dependencies before deleting files
- **Deletes**:
  - `*_scores.output` files (only if scores are extracted to CSV)
  - `*_filtered.output` files (only if `_filtered_kraken.kreport` exists)
  - `*_summary.output` files (always safe, not used)
  - `*.b2` files (only if corresponding merged CSV exists)
  - `*_kraken.output` files (only if `*_kraken.kreport` exists AND script 04 has processed the sample)
- **Preserves**:
  - All `*_kraken.kreport` files
  - All `*_filtered_kraken.kreport` files
  - All `*_merged_lineage.csv` files (final outputs)
  - `filter_results_summary.csv`
- **Usage**: `bash scripts/run_workflow/05_cleanup_intermediate_files.sh`
- **Note**: Run after Steps 3, 4, and 5 to free storage space safely

## File Preservation Strategy

### Never Deleted (unless script 6 has been run):
- `*_kraken.kreport` files (used by summarize scripts, script 6, and for re-runs)
- `*_filtered_kraken.kreport` files (needed for Bracken regeneration and script 6)
- `*_merged_lineage.csv` files (final outputs used by all downstream analysis)
- `filter_results_summary.csv` (extracted scoring information)

### Files Deleted by Cleanup Script (After Verification)
- `*_scores.output` files (only if scores extracted to CSV)
- `*_filtered.output` files (only if kreport exists)
- `*_summary.output` files (not used by any script)
- `*.b2` files (only if merged CSV exists)
- `*_kraken.output` files (only if `*_kraken.kreport` exists AND script 04 has processed the sample)

## Usage Examples

### Running for Different Databases

**01_run_kraken.sh**: Edit the `DBNAME` variable at the top of the script to process different databases. Run once per database.

**02_run_architeuthis.sh**: Edit the `DBNAME` variable at the top of the script (lines 23-26) to process different databases. Run once per database.

**03_add_lineage.sh**: Processes all databases automatically. No changes needed.

**04_reshape_score_reads.r**: Processes all scoring files automatically.

**06_calculate_classification_pct_by_rank.r**: Processes all kreport files automatically.

**05_cleanup_intermediate_files.sh**: Processes all files automatically. Checks dependencies before deletion.

## Dependencies

### Required Before Running Workflow
- Kraken2 databases must be built and available
- Conda environment `struo2` must be activated
- Required R packages: `tidyverse`, `data.table`, `future.apply` (for `04_reshape_score_reads.r`)

### Required Before Cleanup
- `03_add_lineage.sh` must be run to create merged CSV files
- `04_reshape_score_reads.r` must be run to extract scoring information
- `06_calculate_classification_pct_by_rank.r` should be run if you want to enable kreport file deletion (optional)
- `01_calculate_sequencing_depth.r` must be run (required by `04_reshape_score_reads.r` and `06_calculate_classification_pct_by_rank.r`)