#!/bin/bash
# Cleanup script that tracks inputs and outputs for each sample/database/rank
# Only deletes files if the subsequent processed file already exists
# Deletes scoring files only if scores are already in the scoring CSV
#
# Usage: bash scripts/run_workflow/05_cleanup_intermediate_files.sh
#        Or run from project root: cd /projectnb/frpmars/soil_microbe_db/ && bash scripts/run_workflow/05_cleanup_intermediate_files.sh

set -e  # Exit on error

# Set base directory - adjust if running from different location
if [ -d "/projectnb/frpmars/soil_microbe_db" ]; then
    # Running on server
    BASE_DIR="/projectnb/frpmars/soil_microbe_db"
    architeuthis_dir_path="$BASE_DIR/data/NEON_metagenome_classification/02_bracken_output"
    summary_dir="$BASE_DIR/data/NEON_metagenome_classification/summary_files"
    scoring_csv="$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv"
    processed_log="$BASE_DIR/data/classification/analysis_files/filter_results_processed_files.txt"
else
    # Running locally - use relative paths
    architeuthis_dir_path="data/classification/02_bracken_output"
    summary_dir="data/classification/summary_files"
    scoring_csv="data/classification/analysis_files/filter_results_summary.csv"
    processed_log="data/classification/analysis_files/filter_results_processed_files.txt"
fi

# Also check HARDDRIVE location for scoring files
architeuthis_dir_harddrive="/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"

echo "=========================================="
echo "Smart Cleanup of Intermediate Files"
echo "=========================================="
echo ""

# Function to check if a file exists
file_exists() {
    [ -f "$1" ]
}

# Function to check if scoring file is in CSV
score_in_csv() {
    local score_file="$1"
    local samp_name=$(basename "$score_file" | sed 's/_scores\.output$//')
    
    if [ ! -f "$scoring_csv" ]; then
        return 1  # CSV doesn't exist, so score is not in it
    fi
    
    # Check if samp_name appears in the CSV (simple grep check)
    grep -q "$samp_name" "$scoring_csv" 2>/dev/null
}

# Function to check if scoring file is in processed log
score_in_log() {
    local score_file="$1"
    
    if [ ! -f "$processed_log" ]; then
        return 1  # Log doesn't exist
    fi
    
    grep -q "$score_file" "$processed_log" 2>/dev/null
}

# Track statistics
deleted_count=0
skipped_count=0
total_size_freed=0

# ============================================================================
# 1. Clean up _scores.output files (only if scores are in CSV)
# ============================================================================
echo "1. Checking _scores.output files..."
score_files=""

# Check local directory
if [ -d "$architeuthis_dir_path" ]; then
    score_files="$score_files $(find "$architeuthis_dir_path" -name "*_scores.output" -type f 2>/dev/null)"
fi

# Check HARDDRIVE directory if it exists
if [ -d "$architeuthis_dir_harddrive" ]; then
    score_files="$score_files $(find "$architeuthis_dir_harddrive" -name "*_scores.output" -type f 2>/dev/null)"
fi

if [ -n "$score_files" ]; then
    for score_file in $score_files; do
        # Skip empty entries
        [ -z "$score_file" ] && continue
        
        samp_name=$(basename "$score_file" | sed 's/_scores\.output$//')
        
        # Check if score is in CSV or processed log
        if score_in_csv "$score_file" || score_in_log "$score_file"; then
            file_size=$(stat -f%z "$score_file" 2>/dev/null || stat -c%s "$score_file" 2>/dev/null || echo 0)
            rm -f "$score_file"
            deleted_count=$((deleted_count + 1))
            total_size_freed=$((total_size_freed + file_size))
            echo "  ✓ Deleted: $(basename "$score_file") (score extracted)"
        else
            skipped_count=$((skipped_count + 1))
            echo "  ⊘ Skipped: $(basename "$score_file") (score not yet extracted)"
        fi
    done
else
    echo "  No _scores.output files found"
fi
echo ""

# ============================================================================
# 2. Clean up _filtered.output files (only if _filtered_kraken.kreport exists)
# ============================================================================
echo "2. Checking _filtered.output files..."
if [ -d "$architeuthis_dir_path" ]; then
    filtered_files=$(find "$architeuthis_dir_path" -name "*_filtered.output" -type f 2>/dev/null)
    
    if [ -n "$filtered_files" ]; then
        for filtered_file in $filtered_files; do
            samp_name=$(basename "$filtered_file" | sed 's/_filtered\.output$//')
            kreport_file="$architeuthis_dir_path/${samp_name}_filtered_kraken.kreport"
            
            if file_exists "$kreport_file"; then
                file_size=$(stat -f%z "$filtered_file" 2>/dev/null || stat -c%s "$filtered_file" 2>/dev/null || echo 0)
                rm -f "$filtered_file"
                deleted_count=$((deleted_count + 1))
                total_size_freed=$((total_size_freed + file_size))
                echo "  ✓ Deleted: $(basename "$filtered_file") (kreport exists)"
            else
                skipped_count=$((skipped_count + 1))
                echo "  ⊘ Skipped: $(basename "$filtered_file") (kreport not found)"
            fi
        done
    else
        echo "  No _filtered.output files found"
    fi
else
    echo "  Directory not found: $architeuthis_dir_path"
fi
echo ""

# ============================================================================
# 3. Clean up _summary.output files (always safe to delete, not used)
# ============================================================================
echo "3. Checking _summary.output files..."
if [ -d "$architeuthis_dir_path" ]; then
    summary_files=$(find "$architeuthis_dir_path" -name "*_summary.output" -type f 2>/dev/null)
    
    if [ -n "$summary_files" ]; then
        for summary_file in $summary_files; do
            file_size=$(stat -f%z "$summary_file" 2>/dev/null || stat -c%s "$summary_file" 2>/dev/null || echo 0)
            rm -f "$summary_file"
            deleted_count=$((deleted_count + 1))
            total_size_freed=$((total_size_freed + file_size))
            echo "  ✓ Deleted: $(basename "$summary_file") (not used)"
        done
    else
        echo "  No _summary.output files found"
    fi
else
    echo "  Directory not found: $architeuthis_dir_path"
fi
echo ""

# ============================================================================
# 4. Clean up species-level .b2 files (only if merged CSV exists)
# ============================================================================
echo "4. Checking species-level .b2 files..."
databases=("soil_microbe_db" "gtdb_207" "gtdb_207_unfiltered" "pluspf")

for dbname in "${databases[@]}"; do
    # Check for merged lineage file
    merged_file="$summary_dir/${dbname}_species_merged_lineage.csv"
    
    if file_exists "$merged_file"; then
        # Pattern depends on database name
        if [ "$dbname" = "soil_microbe_db" ]; then
            pattern="*${dbname}_filtered.b2"
        elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
            pattern="*gtdb_207_unfiltered.b2"
        elif [ "$dbname" = "gtdb_207" ]; then
            pattern="*gtdb_207_filtered.b2"
        elif [ "$dbname" = "pluspf" ]; then
            pattern="*pluspf_filtered.b2"
        fi
        
        b2_files=$(find "$architeuthis_dir_path" -name "$pattern" -type f 2>/dev/null)
        
        if [ -n "$b2_files" ]; then
            for b2_file in $b2_files; do
                file_size=$(stat -f%z "$b2_file" 2>/dev/null || stat -c%s "$b2_file" 2>/dev/null || echo 0)
                rm -f "$b2_file"
                deleted_count=$((deleted_count + 1))
                total_size_freed=$((total_size_freed + file_size))
            done
            echo "  ✓ Deleted species .b2 files for $dbname (merged CSV exists)"
        fi
    else
        echo "  ⊘ Skipped species .b2 files for $dbname (merged CSV not found)"
    fi
done
echo ""

# ============================================================================
# 5. Clean up phylum-level .b2 files (only if merged CSV exists)
# ============================================================================
echo "5. Checking phylum-level .b2 files..."
for dbname in "${databases[@]}"; do
    # Check for merged lineage file (phylum uses "filtered_phylum" naming)
    merged_file="$summary_dir/${dbname}_filtered_phylum_merged_lineage.csv"
    
    if file_exists "$merged_file"; then
        # Pattern depends on database name
        if [ "$dbname" = "soil_microbe_db" ]; then
            pattern="*${dbname}_phylum_filtered.b2"
        elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
            pattern="*gtdb_207_unfiltered_phylum_filtered.b2"
        elif [ "$dbname" = "gtdb_207" ]; then
            pattern="*gtdb_207_phylum_filtered.b2"
        elif [ "$dbname" = "pluspf" ]; then
            pattern="*pluspf_phylum_filtered.b2"
        fi
        
        b2_files=$(find "$architeuthis_dir_path" -name "$pattern" -type f 2>/dev/null)
        
        if [ -n "$b2_files" ]; then
            for b2_file in $b2_files; do
                file_size=$(stat -f%z "$b2_file" 2>/dev/null || stat -c%s "$b2_file" 2>/dev/null || echo 0)
                rm -f "$b2_file"
                deleted_count=$((deleted_count + 1))
                total_size_freed=$((total_size_freed + file_size))
            done
            echo "  ✓ Deleted phylum .b2 files for $dbname (merged CSV exists)"
        fi
    else
        echo "  ⊘ Skipped phylum .b2 files for $dbname (merged CSV not found)"
    fi
done
echo ""

# ============================================================================
# Summary
# ============================================================================
echo "=========================================="
echo "Cleanup Summary"
echo "=========================================="
echo "Files deleted: $deleted_count"
echo "Files skipped: $skipped_count"
if [ $total_size_freed -gt 0 ]; then
    # Convert bytes to human-readable format
    if command -v numfmt >/dev/null 2>&1; then
        size_freed=$(numfmt --to=iec-i --suffix=B $total_size_freed)
    else
        # Fallback: convert to MB
        size_freed_mb=$((total_size_freed / 1024 / 1024))
        size_freed="${size_freed_mb}MB"
    fi
    echo "Space freed: $size_freed"
else
    echo "Space freed: 0 bytes"
fi
echo ""
