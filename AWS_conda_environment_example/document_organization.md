# Document organization: AWS S3 bucket file descriptions 

This repository includes development and analysis code for the SoilMicrobeDB, as well as a subdirectory with information on how to use the database for your own analysis via an AWS bucket.

The tutorial for using the SoilMicrobeDB is available here: https://github.com/zoey-rw/SoilMicrobeDB/blob/master/AWS_conda_environment_example/tutorial.md

The S3 Bucket contains the following files. The first four are necessary for running Kraken2, and the last three are for running Bracken for abundance estimation for reads of length 150bp.

Kraken2 files:
    kraken/taxo.k2d
    kraken/seqid2taxid.map
    kraken/opts.k2d
    kraken/hash.k2d
    
Bracken files:
    kraken/database150mers.kraken
    kraken/database150mers.kmer_distrib
    kraken/database.kraken
    
SoilMicrobeDB genome list, used as input for the Struo2 database creation pipeline:
        soil_microbe_db_genome_table.csv
    
Test files:
    test_sample_R1.fastq.gz
    test_sample_R2.fastq.gz
    
NEON results file, created from mapping the NEON soil sample collection against the SoilMicrobeDB, and merging all output files at the level of genus abundance:
    soil_microbe_db_filtered_genus_merged_lineage.csv
    
