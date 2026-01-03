# Genome download and processing for database

Scripts for creating and maintaining the Soil Microbe Database for short-read metagenomic classification using Kraken2. Scripts are organized by genome source.

## Directory Structure

```
create_database/
├── README.md 
├── config.R                    # Path configuration and quality settings
├── helper_functions.r          # Shared utilities
└── process_genomes/
    |── gtdb_mapping # Run first to generate mapping keys
    ├── GTDB_207/
    ├── SMAG/                
    ├── SPIRE/                
    ├── Mycocosm/
    ├── JGI_Gold/
    └── GEM/
└── create_soil_genome_db_input.R # Merges all sources into one input file for Struo2 pipeline
```

### Required Execution Order

Scripts must be run in a specific order due to dependencies on shared files and data:

#### Step 1: Create GTDB-NCBI Mapping (REQUIRED FIRST)
**Script**: `process_genomes/gtdb_mapping/create_all_gtdb_mappings.r`

**Purpose**: Creates the unified `GTDB_NCBI_key.rds` mapping file used by all downstream scripts.

**Dependencies**: 
- NCBI taxdump files (2021 and 2023 versions)
- Internet access to download GTDB Excel mapping files

**Output**: `data/genome_database/gtdb_mapping/GTDB_NCBI_key.rds`

**Rationale**: This mapping file is required by all genome source scripts for taxonomy conversion. It contains mappings for GTDB versions 214, 95, and 207, which different sources filter as needed.

#### Step 2: Process Genome Sources (Can run in parallel after Step 1)

The following scripts can be run in any order **after** the GTDB mapping file is created:

- **SMAG**: `process_genomes/SMAG/download_SMAGs.r`
  - Uses: `GTDB_NCBI_key.rds` (all GTDB versions)
  - Optionally uses: GTDB 214 metadata (if available in global environment)
  
- **SPIRE**: `process_genomes/SPIRE/download_spire_MAGs.r`
  - Uses: `GTDB_NCBI_key.rds` (filters to GTDB 207)
  - Uses: GTDB 207 metadata files (`bac120_metadata_r207.tsv`, `ar53_metadata_r207.tsv`)
  
- **GEM/Nayfach**: `process_genomes/GEM/download_nayfach_MAGs.r`
  - Uses: `GTDB_NCBI_key.rds` (all GTDB versions)
  - Uses: GTDB 95 and 214 metadata files (optional, for Strategy 1 mapping)
  
- **GTDB 207**: `process_genomes/GTDB_207/create_GTDB_download_list.r`
  - Uses: GTDB 207 metadata files
  - Note: This script creates a download list but doesn't download genomes
  
- **JGI GOLD**: `process_genomes/JGI_Gold/download_jgi_gold_genomes.r`
  - Uses: `GTDB_NCBI_key.rds` (filters to GTDB 207)
  - Uses: GTDB 207 metadata files
  - Uses: JGI GOLD data files (`goldData.xlsx`, `GOLDs5levelEcosystemClassificationPaths.xlsx`)
  
- **RefSoil**: `process_genomes/RefSoil/download_refsoil_genomes.r`
  - Uses: `GTDB_NCBI_key.rds` (filters to GTDB 207)
  - Uses: GTDB 207 metadata files
  - Uses: JGI GOLD data files (to find matching genomes)
  - **Note**: Should be run after JGI GOLD script if you want to reuse JGI GOLD data files, but this is not strictly required
  
- **Mycocosm**: Two scripts in sequence:
  1. `process_genomes/Mycocosm/download_actual_myco_genomes.r` (downloads genomes)
  2. `process_genomes/Mycocosm/create_mycocosm_database_input.R` (processes downloaded genomes)
  - **Note**: The download script must run before the processing script
  - Uses: NCBI taxonomy files (for processing script)

#### Step 3: Merge All Sources
**Script**: `create_soil_genome_db_input.R`

**Purpose**: Combines output from all genome source scripts into a single input file for the Struo2 pipeline.

**Dependencies**: All genome source scripts must complete successfully first.

### Shared Dependencies

#### GTDB-NCBI Mapping File (`GTDB_NCBI_key.rds`)
- **Created by**: `gtdb_mapping/create_all_gtdb_mappings.r`
- **Used by**: SMAG, SPIRE, GEM, JGI_Gold, RefSoil
- **Rationale**: Provides taxonomy mapping across all GTDB versions (214, 95, 207) and NCBI versions (2021, 2023)

#### GTDB 207 Metadata Files
- **Files**: `bac120_metadata_r207.tsv`, `ar53_metadata_r207.tsv`
- **Used by**: SPIRE, GTDB_207, JGI_Gold, RefSoil
- **Rationale**: Provides direct taxonomy matching for genomes using GTDB 207 taxonomy
- **Note**: These files are large and should be placed in `data/genome_database/` or accessed via mounted remote server

#### GTDB 95/214 Metadata Files (Optional)
- **Files**: `bac120_metadata_r95.tsv`, `ar53_metadata_r95.tsv`, `bac120_metadata_r214.tsv`, `ar53_metadata_r214.tsv`
- **Used by**: GEM (optional, for Strategy 1 mapping), SMAG (optional, checks global environment)
- **Rationale**: Provides additional precision for taxonomy matching, but mapping key is sufficient for most cases

#### JGI GOLD Data Files
- **Files**: `goldData.xlsx`, `GOLDs5levelEcosystemClassificationPaths.xlsx`
- **Used by**: JGI_Gold, RefSoil
- **Rationale**: RefSoil uses JGI GOLD to identify genomes matching RefSoil tax IDs

#### NCBI Taxonomy Files
- **Directory**: `data/genome_database/ncbi_taxonomy/`
- **Used by**: All scripts that add NCBI taxonomy information (SPIRE, JGI_Gold, RefSoil, Mycocosm, GTDB_207)
- **Rationale**: Required for building full NCBI taxonomy strings from tax IDs

### Summary of Execution Order

```
1. gtdb_mapping/create_all_gtdb_mappings.r          [REQUIRED FIRST]
   ↓
2. Genome source scripts (can run in parallel):
   - SMAG/download_SMAGs.r
   - SPIRE/download_spire_MAGs.r
   - GEM/download_nayfach_MAGs.r
   - GTDB_207/create_GTDB_download_list.r
   - JGI_Gold/download_jgi_gold_genomes.r
   - RefSoil/download_refsoil_genomes.r (after JGI_Gold if reusing data files)
   - Mycocosm/download_actual_myco_genomes.r
   - Mycocosm/create_mycocosm_database_input.R (after download script)
   ↓
3. create_soil_genome_db_input.R                  [REQUIRED LAST]
```

## Usage

### Test Mode
All scripts support test mode to limit downloads:
```bash
cd scripts/create_database/process_genomes/SPIRE
TEST_MODE=TRUE MAX_GENOMES=10 Rscript download_spire_MAGs.r
```

### Configuration
All scripts use `config.R`. Source-specific settings are in named lists (SMAG_CONFIG, SPIRE_CONFIG, etc.).

### Notes
- **GTDB mapping**: Single `GTDB_NCBI_key.rds` file contains all GTDB versions (214, 95, 207); scripts filter to relevant version as needed
- **Test mode**: Use `MAX_GENOMES=10` to limit local downloads
- **Mycocosm downloads**: Require RSelenium and browser setup (see `process_genomes/Mycocosm/README.md` for details)
- **GTDB metadata files**: Large files should be placed in `data/genome_database/` or accessed via mounted remote server
