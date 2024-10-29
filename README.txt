# SOIL MICROBE DATABASE (SMDB) - a whole-genome reference database with over 50,000 genomes, focused on soil biodiversity 

Includes over 20,000 metagenome-assembled genomes (MAGs) from soils, and 1400 fungi from Joint Genome Institute's Mycocosm

For full details and genome sources, check out our preprint!

# Usage
Explore/download abundances of each genome across soil samples using our interactive portal!

For local classification against metagenome samples, the database files can be downloaded, along with Snakemake rules for each step of classification (Kraken2 for read classification, Architeuthis for kmer-filtering, and Bracken for abundance estimation)

To add genomes to the database, you will need to download the Struo2 software and have the filepaths of local genomes, along with their NCBI taxonomy 