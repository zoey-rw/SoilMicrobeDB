#!/bin/bash -l
# Run Architeuthis filtering and Bracken abundance estimation
# Processes Kraken2 output files to filter reads and estimate abundances
#
# Usage: Edit DBNAME variable below, then: bash scripts/run_workflow/02_run_architeuthis.sh
#        Or submit as job: qsub scripts/run_workflow/02_run_architeuthis.sh
#
# Input:  *_kraken.output files from Step 1
# Output: *_scores.output, *_filtered_kraken.kreport, *.b2 files

# Set SCC project
#$ -P dietzelab

# 512 GB node
#$ -pe omp 28
#$ -l mem_per_core=18G

#$ -j y # Merge stderr into the stdout file, to reduce clutter
#$ -o /projectnb/frpmars/soil_microbe_db/log_files/run_kraken_mars.log

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
#DBNAME=pluspf
#DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
#DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxo.k2d
#DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
#
#DBNAME=gtdb_207
#DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
#DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
#DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxo.k2d
#
DBNAME=gtdb_207_unfiltered
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy
DB_taxo=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2/taxo.k2d

# Force reprocessing option - set to "true" to regenerate _scores.output files even if they exist
# Useful when files need to be regenerated or are incomplete
# Note: _scores.output files are required for script 4 (reshape_score_reads.r)
FORCE_REPROCESS_SCORES=false


time_with_seconds=$(date +%T)
echo "Beginning Architeuthis/Bracken loop at: $time_with_seconds"

# Compile kraken2-report once before the loop (not on every iteration)
KRAKEN2_REPORT_BIN=/projectnb/frpmars/soil_microbe_db/scripts/kraken2-master/src/kraken2-report
if [ ! -f "$KRAKEN2_REPORT_BIN" ]; then
    echo "Compiling kraken2-report..."
    cd /projectnb/frpmars/soil_microbe_db/scripts/kraken2-master/src
    g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
    if [ ! -f "$KRAKEN2_REPORT_BIN" ]; then
        echo "ERROR: Failed to compile kraken2-report"
        exit 1
    fi
    echo "kraken2-report compiled successfully"
fi

# Set process limits to prevent fork failures
# Increase process limit if possible (check current limit first)
CURRENT_ULIMIT=$(ulimit -u)
echo "Current process limit (ulimit -u): $CURRENT_ULIMIT"
# Try to increase if possible (some systems may not allow this)
ulimit -u 32768 2>/dev/null || ulimit -u 16384 2>/dev/null || ulimit -u 8192 2>/dev/null || true
NEW_ULIMIT=$(ulimit -u)
echo "Process limit set to: $NEW_ULIMIT"

# Process files in batches to avoid hitting process limits
# Batch size: process N files, then wait for all to complete
# Smaller batch size for systems with strict process limits
BATCH_SIZE=5
BATCH_COUNT=0

# Collect all files first to avoid glob expansion issues with many files
# This prevents shell from trying to expand a huge glob pattern
echo "Collecting kraken.output files..."
mapfile -t kraken_files < <(find ${kraken_dir_path} -name "*_kraken.output" -type f | sort)

echo "Found ${#kraken_files[@]} total kraken.output files"
echo "Processing in batches of $BATCH_SIZE files"

# Function to clean up zombie processes
cleanup_processes() {
    # Wait for any background jobs
    wait 2>/dev/null || true
    # Kill any zombie processes (jobs that have finished but not reaped)
    jobs -p 2>/dev/null | xargs -r kill -0 2>/dev/null || true
    # Small delay to allow system to recover
    sleep 0.5
}

# Loop through kraken.output files directly (more reliable than kreport files)
# This ensures we process all samples that have kraken output, even if kreport files are missing
for kraken_file in "${kraken_files[@]}"; do
    
    # Extract samp_name from kraken.output filename
    samp_file_basename="$(basename -- $kraken_file)"
    samp_name=${samp_file_basename::-13}  # Remove "_kraken.output" (13 chars)
    
    # Only process files for the current database
    if [[ $samp_name != *$DBNAME* ]]; then
        continue
    fi
    
    echo "$samp_name architeuthis loop"
    
    KRAKEN_OUTPUT=$kraken_file
    
    ARCHITEUTHIS_FILTERED=${architeuthis_dir_path}/${samp_name}_filtered.output
    ARCHITEUTHIS_SUMMARY=${architeuthis_dir_path}/${samp_name}_summary.output
    ARCHITEUTHIS_SCORES=${architeuthis_dir_path}/${samp_name}_scores.output
    ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport
    
    BRACKEN_OUTPUT=${architeuthis_dir_path}/${samp_name}_filtered.b2
    BRACKEN_OUTPUT_GENUS=${architeuthis_dir_path}/${samp_name}_genus_filtered.b2
    BRACKEN_OUTPUT_DOMAIN=${architeuthis_dir_path}/${samp_name}_domain_filtered.b2
    BRACKEN_OUTPUT_PHYLUM=${architeuthis_dir_path}/${samp_name}_phylum_filtered.b2
    BRACKEN_OUTPUT_CLASS=${architeuthis_dir_path}/${samp_name}_class_filtered.b2
    BRACKEN_OUTPUT_ORDER=${architeuthis_dir_path}/${samp_name}_order_filtered.b2
    BRACKEN_OUTPUT_FAMILY=${architeuthis_dir_path}/${samp_name}_family_filtered.b2
    
    # Run architeuthis mapping filter (creates filtered.output)
    if test -e $ARCHITEUTHIS_FILTERED; then
        echo "$ARCHITEUTHIS_FILTERED exists; skipping this run."
    else
        architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED
    fi
    
    # Run architeuthis mapping score
    # Note: _scores.output files are required for script 4 (reshape_score_reads.r)
    if test -e $ARCHITEUTHIS_SCORES && [ "$FORCE_REPROCESS_SCORES" != "true" ]; then
        echo "$ARCHITEUTHIS_SCORES exists; skipping this run."
    else
        if [ "$FORCE_REPROCESS_SCORES" = "true" ] && test -e $ARCHITEUTHIS_SCORES; then
            echo "FORCE_REPROCESS_SCORES=true: Regenerating $ARCHITEUTHIS_SCORES"
            rm -f $ARCHITEUTHIS_SCORES
        fi
        architeuthis mapping score $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SCORES
    fi
    
    # Convert filtered output to Kraken report format
    NEED_KREPORT_RECONVERSION=false
    
    if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
        echo "$ARCHITEUTHIS_FILTERED_REPORT exists; checking validity..."
        # Remove blank lines from existing file (in case it has them)
        sed -i '/^$/d' $ARCHITEUTHIS_FILTERED_REPORT
        # Check if file is still valid after removing blank lines
        if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT is empty after removing blank lines. Will attempt to recreate."
            rm -f $ARCHITEUTHIS_FILTERED_REPORT
            NEED_KREPORT_RECONVERSION=true
        else
            # Check if kreport has suspiciously few reads compared to filtered.output
            # This indicates the conversion may have failed
            # Note: filtered.output may have been deleted, so check if it still exists
            if [ -s "$ARCHITEUTHIS_FILTERED" ]; then
                filtered_reads=$(grep -c "^C" "$ARCHITEUTHIS_FILTERED" 2>/dev/null || echo "0")
                kreport_reads=$(awk '{sum += $2} END {print sum}' "$ARCHITEUTHIS_FILTERED_REPORT" 2>/dev/null || echo "0")
                
                # If kreport has < 1% of filtered reads, something went wrong with conversion
                if [ "$filtered_reads" -gt 1000 ] && [ "$kreport_reads" -gt 0 ]; then
                    percent=$(awk "BEGIN {printf \"%.2f\", ($kreport_reads / $filtered_reads) * 100}")
                    threshold=$(awk "BEGIN {printf \"%.0f\", $filtered_reads * 0.01}")
                    if [ "$kreport_reads" -lt "$threshold" ]; then
                        echo "  WARNING: Kreport has only $kreport_reads reads vs $filtered_reads in filtered.output ($percent%)"
                        echo "  This suggests kreport conversion failed. Will attempt to recreate."
                        rm -f $ARCHITEUTHIS_FILTERED_REPORT
                        NEED_KREPORT_RECONVERSION=true
                    fi
                elif [ "$filtered_reads" -gt 1000 ] && [ "$kreport_reads" -eq 0 ]; then
                    echo "  WARNING: Kreport has 0 reads but filtered.output has $filtered_reads reads"
                    echo "  This suggests kreport conversion failed. Will attempt to recreate."
                    rm -f $ARCHITEUTHIS_FILTERED_REPORT
                    NEED_KREPORT_RECONVERSION=true
                fi
            else
                # filtered.output was deleted, but kreport exists - check if kreport looks valid
                kreport_reads=$(awk '{sum += $2} END {print sum}' "$ARCHITEUTHIS_FILTERED_REPORT" 2>/dev/null || echo "0")
                if [ "$kreport_reads" -eq 0 ]; then
                    echo "  WARNING: Kreport has 0 reads and filtered.output is missing. Cannot re-convert."
                fi
            fi
        fi
    else
        NEED_KREPORT_RECONVERSION=true
    fi
    
    if [ "$NEED_KREPORT_RECONVERSION" = "true" ]; then
        # Check if filtered.output exists and is not empty
        if [ ! -s "$ARCHITEUTHIS_FILTERED" ]; then
            echo "  WARNING: $ARCHITEUTHIS_FILTERED is empty or does not exist. Cannot create kreport."
        else
            echo "  Converting filtered.output to kreport format..."
            # Convert back to Kraken output using pre-compiled binary
            $KRAKEN2_REPORT_BIN $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT
            
            # Validate kreport file was created and is not empty
            if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
                echo "  ERROR: $ARCHITEUTHIS_FILTERED_REPORT is empty or was not created properly."
            else
                # Remove leading blank lines and ensure file starts with valid content
                sed -i '/^$/d' $ARCHITEUTHIS_FILTERED_REPORT
                # Check if file is still valid after removing blank lines
                if [ ! -s "$ARCHITEUTHIS_FILTERED_REPORT" ]; then
                    echo "  ERROR: $ARCHITEUTHIS_FILTERED_REPORT is empty after removing blank lines."
                else
                    # Verify conversion was successful by checking read counts
                    filtered_reads=$(grep -c "^C" "$ARCHITEUTHIS_FILTERED" 2>/dev/null || echo "0")
                    kreport_reads=$(awk '{sum += $2} END {print sum}' "$ARCHITEUTHIS_FILTERED_REPORT" 2>/dev/null || echo "0")
                    echo "  Conversion complete: filtered.output has $filtered_reads reads, kreport has $kreport_reads reads"
                    
                    # Check if conversion looks reasonable
                    # Kreport read count should be similar to filtered read count (within reasonable range)
                    # If kreport has < 1% of filtered reads and filtered has > 1000 reads, conversion likely failed
                    if [ "$filtered_reads" -gt 1000 ] && [ "$kreport_reads" -gt 0 ]; then
                        percent=$(awk "BEGIN {printf \"%.2f\", ($kreport_reads / $filtered_reads) * 100}")
                        if [ "$(echo "$percent < 1.0" | bc 2>/dev/null || echo "0")" = "1" ]; then
                            echo "  ERROR: Kreport conversion appears to have failed (only $percent% of reads converted)"
                            echo "  Keeping filtered.output for potential re-conversion"
                            # Don't delete filtered.output if conversion failed
                        else
                            # Clean up ARCHITEUTHIS_FILTERED after successful kreport creation
                            # ARCHITEUTHIS_SCORES will be processed by 04_reshape_score_reads.r to extract scoring information
                            rm -f $ARCHITEUTHIS_FILTERED
                        fi
                    else
                        # For small files, just clean up
                        rm -f $ARCHITEUTHIS_FILTERED
                    fi
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
            # Check if kreport has any reads (sum of column 2 should be > 0)
            # kreport format: percentage, cladeReads, taxonReads, ...
            total_reads=$(awk '{sum += $2} END {print sum}' "$ARCHITEUTHIS_FILTERED_REPORT" 2>/dev/null)
            if [ -z "$total_reads" ] || [ "$total_reads" = "0" ] || [ "$total_reads" = "" ]; then
                echo "  WARNING: $ARCHITEUTHIS_FILTERED_REPORT has 0 reads. Skipping Bracken for $samp_name"
                echo "  This may indicate that architeuthis filter removed all reads, or the kreport conversion failed."
            else
                # Valid kreport format with reads - proceed with Bracken
                # Check and run each rank independently (ensures all ranks are generated even if some already exist)
                # Note: Only check for .b2 files (not kreport files) - .b2 is the required format for architeuthis merge
                if test -e $BRACKEN_OUTPUT; then
                    echo "$BRACKEN_OUTPUT exists; skipping species-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
                fi
                
                if test -e $BRACKEN_OUTPUT_GENUS; then
                    echo "$BRACKEN_OUTPUT_GENUS exists; skipping genus-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_GENUS -l G
                fi
                
                if test -e $BRACKEN_OUTPUT_DOMAIN; then
                    echo "$BRACKEN_OUTPUT_DOMAIN exists; skipping domain-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_DOMAIN -l D
                fi
                
                if test -e $BRACKEN_OUTPUT_PHYLUM; then
                    echo "$BRACKEN_OUTPUT_PHYLUM exists; skipping phylum-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_PHYLUM -l P
                fi
                
                if test -e $BRACKEN_OUTPUT_CLASS; then
                    echo "$BRACKEN_OUTPUT_CLASS exists; skipping class-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_CLASS -l C
                fi
                
                if test -e $BRACKEN_OUTPUT_ORDER; then
                    echo "$BRACKEN_OUTPUT_ORDER exists; skipping order-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_ORDER -l O
                fi
                
                if test -e $BRACKEN_OUTPUT_FAMILY; then
                    echo "$BRACKEN_OUTPUT_FAMILY exists; skipping family-level Bracken."
                else
                    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_FAMILY -l F
                fi
            fi
        fi
    fi
    
    # Batch processing: wait every BATCH_SIZE files to prevent process accumulation
    BATCH_COUNT=$((BATCH_COUNT + 1))
    if [ $((BATCH_COUNT % BATCH_SIZE)) -eq 0 ]; then
        echo "Processed $BATCH_COUNT files, cleaning up processes..."
        cleanup_processes
        echo "Continuing with next batch..."
    fi
done

# Final cleanup
echo "Final cleanup of processes..."
cleanup_processes

time_with_seconds=$(date +%T)
echo "Finished Architeuthis/Bracken loop at: $time_with_seconds"

# Note: Cleanup of intermediate files is handled by scripts/run_workflow/05_cleanup_intermediate_files.sh
# - Deletes _scores.output only if scores are extracted to CSV
# - Deletes _filtered.output only if _filtered_kraken.kreport exists
# - Deletes .b2 files only if merged CSV files exist
# Run: bash scripts/run_workflow/05_cleanup_intermediate_files.sh
