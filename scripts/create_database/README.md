# Database Creation Scripts

Scripts for creating and maintaining the Soil Microbe Database. Scripts are organized by genome source.

## Directory Structure

```
create_database/
├── README.md (this file)
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
- **GTDB metadata required**. Will be downloaded automatically within process_genomes/gtdb_mapping/create_all_gtdb_mappings.r
- **GTDB mapping**: Single `GTDB_NCBI_key.rds` file contains all GTDB versions (214, 95, 207); scripts filter to relevant version as needed
- **Test mode**: Use `MAX_GENOMES=10` to limit local downloads
