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
DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxo.k2d
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
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
    
    # Run architeuthis mapping filter (creates filtered.output)
    if test -e $ARCHITEUTHIS_FILTERED; then
        echo "$ARCHITEUTHIS_FILTERED exists; skipping this run."
    else
        architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED
    fi
    
    # Run architeuthis mapping score
    if test -e $ARCHITEUTHIS_SCORES; then
        echo "$ARCHITEUTHIS_SCORES exists; skipping this run."
    else
        architeuthis mapping score $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SCORES
    fi
    
    # Convert filtered output to Kraken report format
    if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
        echo "$ARCHITEUTHIS_FILTERED_REPORT exists; checking validity..."
        # Remove blank lines from existing file (in case it has them)
        sed -i '/^$/d' $ARCHITEUTHIS_FILTERED_REPORT
        # Check if file is still valid after removing blank lines
        if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT is empty after removing blank lines. Will attempt to recreate."
            rm -f $ARCHITEUTHIS_FILTERED_REPORT
        fi
    fi
    
    if [ ! -e "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
        # Check if filtered.output exists and is not empty
        if [ ! -s "$ARCHITEUTHIS_FILTERED" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED is empty or does not exist. Skipping kreport conversion."
        else
            # Convert back to Kraken output
            cd /projectnb/frpmars/soil_microbe_db/scripts/kraken2-master/src
            g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
            ./kraken2-report $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT
            
            # Validate kreport file was created and is not empty
            if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
                echo "  ERROR: $ARCHITEUTHIS_FILTERED_REPORT is empty or was not created properly."
            else
                # Remove leading blank lines and ensure file starts with valid content
                # Bracken expects first line to start with C or U
                sed -i '/^$/d' $ARCHITEUTHIS_FILTERED_REPORT
                # Check if file is still valid after removing blank lines
                if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
                    echo "  ERROR: $ARCHITEUTHIS_FILTERED_REPORT is empty after removing blank lines."
                else
                    # Clean up ARCHITEUTHIS_FILTERED after kreport is created (not needed after kreport is created)
                    # ARCHITEUTHIS_SCORES will be processed by 04_reshape_score_reads.r to extract scoring information
                    rm -f $ARCHITEUTHIS_FILTERED
                fi
            fi
        fi
    fi
    
    # Run Bracken at species, genus, domain, and phylum levels
    # Only run if kreport file exists and is valid
    if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
        echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT is empty or missing. Skipping Bracken for $samp_name"
    else
        # Check if first line of kreport is valid (should NOT start with C or U - that would be .output format)
        # Valid kreport format starts with percentage (number) or space, not C/U
        first_char=$(head -n 1 "$ARCHITEUTHIS_FILTERED_REPORT" | cut -c1)
        if [ -z "$first_char" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT has empty first line. Skipping Bracken for $samp_name"
        elif [ "$first_char" = "C" ] || [ "$first_char" = "U" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT appears to be in .output format (starts with '$first_char'), not .kreport format. Skipping Bracken for $samp_name"
        else
            # Valid kreport format - proceed with Bracken
            if test -e $BRACKEN_OUTPUT; then
                echo "$BRACKEN_OUTPUT exists; skipping this run."
            else
                /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
                /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_GENUS -l G
                /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_DOMAIN -l D
                /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_PHYLUM -l P
            fi
        fi
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
