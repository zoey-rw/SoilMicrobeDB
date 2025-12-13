#!/bin/bash -l

# Set SCC project
#$ -P dietzelab

# 512 GB node
#$ -pe omp 28
#$ -l mem_per_core=18G

#$ -j y # Merge stderr into the stdout file, to reduce clutter
#$ -o /projectnb/frpmars/soil_microbe_db/log_files/run_kraken_mars.log

# load kraken2

module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/
architeuthis_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/02_bracken_output
kraken_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/01_kraken_output

# Set database to process - change this for different databases
# Options: soil_microbe_db, pluspf, gtdb_207, gtdb_207_unfiltered
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxo.k2d

# Uncomment and modify for other databases:
# DBNAME=pluspf
# DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
# DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxo.k2d
# DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
#
# DBNAME=gtdb_207
# DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
# DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
# DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxo.k2d
#
# DBNAME=gtdb_207_unfiltered
# DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
# DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy
# DB_taxo=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2/taxo.k2d


time_with_seconds=$(date +%T)
echo "Beginning Architeuthis/Bracken loop at: $time_with_seconds"

for samp_file in /projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/01_kraken_output/*.kreport; do

    samp_file_basename="$(basename -- $samp_file)"
    
    samp_name=${samp_file_basename::-15}
    echo "$samp_name architeuthis loop"
    
    # Only process files for the current database
    if [[ $samp_name != *$DBNAME* ]]; then
        continue
    fi
    
    KRAKEN_OUTPUT=${kraken_dir_path}/${samp_name}_kraken.output
    
    ARCHITEUTHIS_FILTERED=${architeuthis_dir_path}/${samp_name}_filtered.output
    ARCHITEUTHIS_SUMMARY=${architeuthis_dir_path}/${samp_name}_summary.output
    ARCHITEUTHIS_SCORES=${architeuthis_dir_path}/${samp_name}_scores.output
    ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport
    
    BRACKEN_OUTPUT=${architeuthis_dir_path}/${samp_name}_filtered.b2
    BRACKEN_OUTPUT_GENUS=${architeuthis_dir_path}/${samp_name}_genus_filtered.b2
    BRACKEN_OUTPUT_DOMAIN=${architeuthis_dir_path}/${samp_name}_domain_filtered.b2
    BRACKEN_OUTPUT_PHYLUM=${architeuthis_dir_path}/${samp_name}_phylum_filtered.b2
    
    # Run architeuthis mapping score
    if test -e $ARCHITEUTHIS_SCORES; then
        echo "$ARCHITEUTHIS_SCORES exists; skipping this run."
    else
        # Run filter
        #architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED 
        #architeuthis mapping summary $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SUMMARY
        architeuthis mapping score $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SCORES
    fi
    
    # Convert filtered output to Kraken report format
    if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
        echo "$ARCHITEUTHIS_FILTERED_REPORT exists; skipping this run."
    else
        # Convert back to Kraken output
        cd /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/misc_scripts/mock_community/kraken2_report/daydream-boost-kraken2-e533c50/src
        g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
        ./kraken2-report $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT
        
        # Clean up ARCHITEUTHIS_FILTERED after kreport is created (not needed after kreport is created)
        # ARCHITEUTHIS_SCORES will be processed by 04_reshape_score_reads.r to extract scoring information
        if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
            rm -f $ARCHITEUTHIS_FILTERED
        fi
    fi
    
    # Run Bracken at species, genus, domain, and phylum levels
    if test -e $BRACKEN_OUTPUT; then
        echo "$BRACKEN_OUTPUT exists; skipping this run."
    else
        /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
        /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_GENUS -l G
        /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_DOMAIN -l D
        /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_PHYLUM -l P
    fi
done

time_with_seconds=$(date +%T)
echo "Finished Architeuthis/Bracken loop at: $time_with_seconds"

# Note: Cleanup of intermediate files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh
# This script intelligently checks dependencies before deleting files:
# - Deletes _scores.output only if scores are extracted to CSV
# - Deletes _filtered.output only if _filtered_kraken.kreport exists
# - Deletes .b2 files only if merged CSV files exist
# Run: bash scripts/run_workflow/05_cleanup_intermediate_files.sh
