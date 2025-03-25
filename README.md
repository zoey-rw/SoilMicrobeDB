# Soil Microbe Database (SMDB) - a whole-genome reference database with over 50,000 genomes, focused on soil biodiversity 

Includes over 20,000 metagenome-assembled genomes (MAGs) from soils, and 1400 fungi from Joint Genome Institute's Mycocosm

For full details and genome sources, check out our preprint: https://www.biorxiv.org/content/10.1101/2025.03.21.644662v1

## Usage
Explore/download abundances of a subset of these genomes across NEON soil samples using the interactive portal, developed in partnership with the Boston University Software & Application Innovation Lab (SAIL): https://zoeywerbin.shinyapps.io/soil_microbe_db/

For local classification against metagenome samples, the database index files can be downloaded from an AWS bucket, along with Snakemake rules for each step of classification (Kraken2 for read classification, Architeuthis for kmer-filtering, and Bracken for abundance estimation). This is hosted through the AWS Open Data program - details will be added here as soon as available.

To add genomes to the database, you will need to download the database along with the Struo2 software, which contains a "db-update" workflow. The input file must have the filepaths of local genomes, along with their NCBI taxonomy, for compatibility with SoilGenomeDB.

## Credits & Acknowledgements

The National Ecological Observatory Network is a project solely funded by the National Science Foundation and managed under cooperative agreement by Battelle. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

## Disclaimer
Information and documents contained within this repository are available as-is. Codes or documents, or their use, may not be supported or maintained under any program or service and may not be compatible with data currently available from the NEON Data Portal.
