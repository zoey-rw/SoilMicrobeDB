# Reference Data Files

This directory contains reference data files used for building the Soil Microbe Database.

## Files Available on GitHub

The following small reference files are automatically downloaded from GitHub if not found locally:

- **GTDB Mapping**: `gtdb_mapping/GTDB_NCBI_key.rds` (1.0M)
- **SPIRE Taxonomy**: `SPIRE/spire.per_cluster.taxonomy.tsv` (15M)
- **SPIRE Microntology**: `SPIRE/spire_v1_microntology.tsv` (39M)
- **GEM Supplementary**: `GEM/41587_2020_718_MOESM3_ESM.xlsx` (33M)
- **RefSoil Data**: `RefSoil/ref_data/RefSoil_s1.xlsx` (144K)

## Large Files (Not on GitHub)

The following files exceed GitHub's 100MB file size limit and must be obtained separately:

### JGI GOLD Files
- `JGI_Gold/ref_data/goldData.xlsx` (224M)
  - Download from: https://gold.jgi.doe.gov/download?mode=site_excel
  - Place in: `data/genome_database/JGI_Gold/ref_data/`

- `JGI_Gold/ref_data/GOLDs5levelEcosystemClassificationPaths.xlsx` (208M)
  - Download from: https://gold.jgi.doe.gov/download?mode=site_excel
  - Place in: `data/genome_database/JGI_Gold/ref_data/`

### SPIRE Metadata
- `SPIRE/spire_v1_genome_metadata.tsv.gz` (96M)
  - Download from: http://spire.embl.de/
  - Place in: `data/genome_database/SPIRE/`

### GTDB Metadata Files
- `bac120_metadata_r207.tsv` (large)
- `ar53_metadata_r207.tsv` (large)
- `bac120_metadata_r214.tsv` (large)
- `bac120_metadata_r95.tsv` (large)
  - Download from: https://gtdb.ecogenomic.org/
  - Place in: `data/genome_database/`

## Setup Instructions

1. **Small files**: Will be automatically downloaded from GitHub when running scripts (if configured)
2. **Large files**: Download manually from the sources listed above and place in the appropriate directories
3. **Environment variable**: Set `SOIL_MICROBE_DB_GITHUB` to override the default GitHub repository URL

## Directory Structure

```
data/genome_database/
├── gtdb_mapping/          # GTDB-NCBI mapping files
├── SPIRE/                 # SPIRE MAG metadata and genomes
├── GEM/                   # GEM/Nayfach supplementary data
├── JGI_Gold/              # JGI GOLD reference data
│   └── ref_data/
├── RefSoil/               # RefSoil reference data
│   └── ref_data/
└── ncbi_genomes/          # NCBI genome downloads
```


