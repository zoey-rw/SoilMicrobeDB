#!/bin/bash
# Cleanup intermediate files from classification pipeline
# Intelligently checks dependencies before deleting files
# Only deletes files if the subsequent processed file already exists
#
# Usage: bash scripts/run_workflow/05_cleanup_intermediate_files.sh
#        Or run from project root: cd /projectnb/frpmars/soil_microbe_db/ && bash scripts/run_workflow/05_cleanup_intermediate_files.sh
#
# Deletes:
#   - *_scores.output files (only if scores extracted to CSV AND CSV is valid)
#   - *_filtered.output files (only if _filtered_kraken.kreport exists AND is valid)
#   - *.b2 files (only if merged CSV files exist AND are valid)
#
# Safety checks:
#   - Files must be non-empty
#   - Files must not have been modified in the last 5 minutes (may still be writing)
#   - Files must be readable

# Don't use set -e here - we want to continue even if some checks fail
# Individual operations use error handling to continue processing

# Set base directory - adjust if running from different location
if [ -d "/projectnb/frpmars/soil_microbe_db" ]; then
    # Running on server
    BASE_DIR="/projectnb/frpmars/soil_microbe_db"
    architeuthis_dir_path="$BASE_DIR/data/NEON_metagenome_classification/02_bracken_output"
    summary_dir="$BASE_DIR/data/NEON_metagenome_classification/summary_files"
    # Check both possible locations for scoring CSV (script 4 uses summary_files)
    if [ -f "$BASE_DIR/data/summary_files/filter_results_summary.csv" ]; then
        scoring_csv="$BASE_DIR/data/summary_files/filter_results_summary.csv"
        processed_log="$BASE_DIR/data/summary_files/filter_results_processed_files.txt"
    elif [ -f "$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv" ]; then
        scoring_csv="$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv"
        processed_log="$BASE_DIR/data/classification/analysis_files/filter_results_processed_files.txt"
    else
        scoring_csv="$BASE_DIR/data/summary_files/filter_results_summary.csv"
        processed_log="$BASE_DIR/data/summary_files/filter_results_processed_files.txt"
    fi
else
    # Running locally - use relative paths, check both locations
    architeuthis_dir_path="data/NEON_metagenome_classification/02_bracken_output"
    summary_dir="data/NEON_metagenome_classification/summary_files"
    # Check both possible locations for scoring CSV
    if [ -f "data/summary_files/filter_results_summary.csv" ]; then
        scoring_csv="data/summary_files/filter_results_summary.csv"
        processed_log="data/summary_files/filter_results_processed_files.txt"
    elif [ -f "data/classification/analysis_files/filter_results_summary.csv" ]; then
        scoring_csv="data/classification/analysis_files/filter_results_summary.csv"
        processed_log="data/classification/analysis_files/filter_results_processed_files.txt"
    else
        scoring_csv="data/summary_files/filter_results_summary.csv"
        processed_log="data/summary_files/filter_results_processed_files.txt"
    fi
fi

# Also check HARDDRIVE location for scoring files
architeuthis_dir_harddrive="/Volumes/HARDDRIVE/SoilMicrobeDB/data/classification/02_bracken_output"

echo "=========================================="
echo "Smart Cleanup of Intermediate Files"
echo "=========================================="
echo ""

# Verify critical directories exist
if [ ! -d "$architeuthis_dir_path" ] && [ ! -d "$architeuthis_dir_harddrive" ]; then
    echo "WARNING: Neither architeuthis directory found:"
    echo "  Checked: $architeuthis_dir_path"
    echo "  Checked: $architeuthis_dir_harddrive"
    echo "  Continuing anyway (may be running before Step 2)..."
fi

if [ ! -d "$summary_dir" ]; then
    echo "WARNING: Summary directory not found: $summary_dir"
    echo "  This is expected if Step 3 hasn't been run yet."
    echo "  Continuing anyway..."
fi

echo "Using paths:"
echo "  Architeuthis dir: $architeuthis_dir_path"
[ -d "$architeuthis_dir_harddrive" ] && echo "  HARDDRIVE dir: $architeuthis_dir_harddrive"
echo "  Summary dir: $summary_dir"
echo "  Scoring CSV: $scoring_csv"
echo "  Processed log: $processed_log"
echo ""

# Function to normalize samp_name (remove trailing underscores)
normalize_samp_name() {
    local samp_name="$1"
    # Remove trailing underscores
    samp_name=$(echo "$samp_name" | sed 's/_\+$//')
    # Remove trailing underscores after database names
    samp_name=$(echo "$samp_name" | sed 's/\(_gtdb_207_unfiltered\|_gtdb_207\|_soil_microbe_db\|_pluspf\)_\+$/\1/')
    echo "$samp_name"
}

# Function to check if a file exists
file_exists() {
    [ -f "$1" ]
}

# Check if file is valid and complete (not empty, not recently modified, readable)
file_is_valid() {
    local file_path="$1"
    local grace_period_minutes="${2:-5}"  # Default 5 minute grace period
    
    # Check if file exists
    if [ ! -f "$file_path" ]; then
        return 1
    fi
    
    # Check if file is readable
    if [ ! -r "$file_path" ]; then
        return 1
    fi
    
    # Check if file is non-empty
    file_size=$(stat -f%z "$file_path" 2>/dev/null || stat -c%s "$file_path" 2>/dev/null || echo 0)
    if [ "$file_size" -eq 0 ]; then
        return 1  # File is empty
    fi
    
    # Check if file was recently modified (might still be writing)
    if command -v stat >/dev/null 2>&1; then
        # Get file modification time
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            mod_time=$(stat -f %m "$file_path" 2>/dev/null)
        else
            # Linux
            mod_time=$(stat -c %Y "$file_path" 2>/dev/null)
        fi
        
        if [ -n "$mod_time" ]; then
            current_time=$(date +%s)
            age_seconds=$((current_time - mod_time))
            age_minutes=$((age_seconds / 60))
            
            # If file was modified less than grace_period_minutes ago, consider it still being written
            if [ "$age_minutes" -lt "$grace_period_minutes" ]; then
                return 1  # File too recently modified
            fi
        fi
    fi
    
    # For CSV files, check if they have at least one data line (not just header)
    if [[ "$file_path" == *.csv ]]; then
        # Count non-empty lines (skip header)
        line_count=$(tail -n +2 "$file_path" 2>/dev/null | grep -v '^$' | wc -l | tr -d ' ')
        if [ "${line_count:-0}" -eq 0 ]; then
            return 1  # CSV has no data rows
        fi
    fi
    
    return 0  # File is valid
}

# Function to check if scoring file is in CSV
score_in_csv() {
    local score_file="$1"
    local samp_name=$(basename "$score_file" | sed 's/_scores\.output$//')
    samp_name=$(normalize_samp_name "$samp_name")
    
    # Check both possible CSV locations (relative and absolute paths)
    local csv_locations=("$scoring_csv")
    
    # Check relative paths
    if [ -f "data/summary_files/filter_results_summary.csv" ] && [ "$scoring_csv" != "data/summary_files/filter_results_summary.csv" ]; then
        csv_locations+=("data/summary_files/filter_results_summary.csv")
    fi
    if [ -f "data/classification/analysis_files/filter_results_summary.csv" ] && [ "$scoring_csv" != "data/classification/analysis_files/filter_results_summary.csv" ]; then
        csv_locations+=("data/classification/analysis_files/filter_results_summary.csv")
    fi
    
    # Check absolute server paths if on server
    if [ -d "/projectnb/frpmars/soil_microbe_db" ]; then
        BASE_DIR="/projectnb/frpmars/soil_microbe_db"
        if [ -f "$BASE_DIR/data/summary_files/filter_results_summary.csv" ] && [ "$scoring_csv" != "$BASE_DIR/data/summary_files/filter_results_summary.csv" ]; then
            csv_locations+=("$BASE_DIR/data/summary_files/filter_results_summary.csv")
        fi
        if [ -f "$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv" ] && [ "$scoring_csv" != "$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv" ]; then
            csv_locations+=("$BASE_DIR/data/classification/analysis_files/filter_results_summary.csv")
        fi
    fi
    
    for csv_file in "${csv_locations[@]}"; do
        if [ -f "$csv_file" ]; then
            # Check if samp_name appears in the CSV (normalized version)
            # samp_name values are unique identifiers, so a simple grep is safe
            # Also check for version with trailing underscore in case CSV has it
            if grep -q "$samp_name" "$csv_file" 2>/dev/null || grep -q "${samp_name}_" "$csv_file" 2>/dev/null; then
                return 0  # Found in CSV
            fi
        fi
    done
    
    return 1  # Not found in any CSV
}

# Function to check if scoring file is in processed log
score_in_log() {
    local score_file="$1"
    
    # Check both possible log locations (relative and absolute paths)
    local log_locations=("$processed_log")
    
    # Check relative paths
    if [ -f "data/summary_files/filter_results_processed_files.txt" ] && [ "$processed_log" != "data/summary_files/filter_results_processed_files.txt" ]; then
        log_locations+=("data/summary_files/filter_results_processed_files.txt")
    fi
    if [ -f "data/classification/analysis_files/filter_results_processed_files.txt" ] && [ "$processed_log" != "data/classification/analysis_files/filter_results_processed_files.txt" ]; then
        log_locations+=("data/classification/analysis_files/filter_results_processed_files.txt")
    fi
    
    # Check absolute server paths if on server
    if [ -d "/projectnb/frpmars/soil_microbe_db" ]; then
        BASE_DIR="/projectnb/frpmars/soil_microbe_db"
        if [ -f "$BASE_DIR/data/summary_files/filter_results_processed_files.txt" ] && [ "$processed_log" != "$BASE_DIR/data/summary_files/filter_results_processed_files.txt" ]; then
            log_locations+=("$BASE_DIR/data/summary_files/filter_results_processed_files.txt")
        fi
        if [ -f "$BASE_DIR/data/classification/analysis_files/filter_results_processed_files.txt" ] && [ "$processed_log" != "$BASE_DIR/data/classification/analysis_files/filter_results_processed_files.txt" ]; then
            log_locations+=("$BASE_DIR/data/classification/analysis_files/filter_results_processed_files.txt")
        fi
    fi
    
    for log_file in "${log_locations[@]}"; do
        if [ -f "$log_file" ]; then
            if grep -q "$score_file" "$log_file" 2>/dev/null; then
                return 0  # Found in log
            fi
        fi
    done
    
    return 1  # Not found in any log
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
        samp_name=$(normalize_samp_name "$samp_name")
        
        # Only delete if score is confirmed in CSV (not just in processed log)
        # The CSV is the source of truth - processed log is just for tracking
        if score_in_csv "$score_file"; then
            # Verify CSV file is valid before deleting score file
            if file_is_valid "$scoring_csv" 5; then
                file_size=$(stat -f%z "$score_file" 2>/dev/null || stat -c%s "$score_file" 2>/dev/null || echo 0)
                rm -f "$score_file"
                deleted_count=$((deleted_count + 1))
                total_size_freed=$((total_size_freed + file_size))
                echo "  ✓ Deleted: $(basename "$score_file") (score confirmed in CSV)"
            else
                skipped_count=$((skipped_count + 1))
                echo "  ⊘ Skipped: $(basename "$score_file") (CSV file not valid or recently modified)"
            fi
        else
            skipped_count=$((skipped_count + 1))
            echo "  ⊘ Skipped: $(basename "$score_file") (score not yet in CSV)"
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
filtered_files=""

# Check both local and HARDDRIVE locations
if [ -d "$architeuthis_dir_path" ]; then
    filtered_files="$filtered_files $(find "$architeuthis_dir_path" -name "*_filtered.output" -type f 2>/dev/null)"
fi

if [ -d "$architeuthis_dir_harddrive" ]; then
    filtered_files="$filtered_files $(find "$architeuthis_dir_harddrive" -name "*_filtered.output" -type f 2>/dev/null)"
fi

if [ -n "$filtered_files" ]; then
    for filtered_file in $filtered_files; do
        # Skip empty entries
        [ -z "$filtered_file" ] && continue
        
        raw_samp_name=$(basename "$filtered_file" | sed 's/_filtered\.output$//')
        samp_name=$(normalize_samp_name "$raw_samp_name")
        
        # Check for kreport in same directory as filtered file
        # Try normalized version first, then version with trailing underscore
        filtered_dir=$(dirname "$filtered_file")
        kreport_file="$filtered_dir/${samp_name}_filtered_kraken.kreport"
        
        # If normalized version doesn't exist, try with trailing underscore
        if [ ! -f "$kreport_file" ] && [[ "$raw_samp_name" != "$samp_name" ]]; then
            kreport_file_alt="$filtered_dir/${raw_samp_name}_filtered_kraken.kreport"
            if [ -f "$kreport_file_alt" ]; then
                kreport_file="$kreport_file_alt"
            fi
        fi
        
        # Verify kreport file is valid (non-empty, not recently modified) before deleting filtered.output
        if file_is_valid "$kreport_file" 5; then
            file_size=$(stat -f%z "$filtered_file" 2>/dev/null || stat -c%s "$filtered_file" 2>/dev/null || echo 0)
            rm -f "$filtered_file"
            deleted_count=$((deleted_count + 1))
            total_size_freed=$((total_size_freed + file_size))
            echo "  ✓ Deleted: $(basename "$filtered_file") (kreport exists and is valid)"
        elif file_exists "$kreport_file"; then
            skipped_count=$((skipped_count + 1))
            echo "  ⊘ Skipped: $(basename "$filtered_file") (kreport exists but may still be writing)"
        else
            skipped_count=$((skipped_count + 1))
            echo "  ⊘ Skipped: $(basename "$filtered_file") (kreport not found: $(basename "$kreport_file"))"
        fi
    done
else
    echo "  No _filtered.output files found"
fi
echo ""

# ============================================================================
# 3. Clean up _summary.output files (always safe to delete, not used)
# ============================================================================
echo "3. Checking _summary.output files..."
summary_files=""

# Check both local and HARDDRIVE locations
if [ -d "$architeuthis_dir_path" ]; then
    summary_files="$summary_files $(find "$architeuthis_dir_path" -name "*_summary.output" -type f 2>/dev/null)"
fi

if [ -d "$architeuthis_dir_harddrive" ]; then
    summary_files="$summary_files $(find "$architeuthis_dir_harddrive" -name "*_summary.output" -type f 2>/dev/null)"
fi

if [ -n "$summary_files" ]; then
    for summary_file in $summary_files; do
        # Skip empty entries
        [ -z "$summary_file" ] && continue
        
        file_size=$(stat -f%z "$summary_file" 2>/dev/null || stat -c%s "$summary_file" 2>/dev/null || echo 0)
        rm -f "$summary_file"
        deleted_count=$((deleted_count + 1))
        total_size_freed=$((total_size_freed + file_size))
        echo "  ✓ Deleted: $(basename "$summary_file") (not used)"
    done
else
    echo "  No _summary.output files found"
fi
echo ""

# ============================================================================
# 4. Clean up .b2 files (only if merged CSV exists) - All ranks
# ============================================================================
echo "4. Checking .b2 files for all ranks..."
databases=("soil_microbe_db" "gtdb_207" "gtdb_207_unfiltered" "pluspf")
ranks=("species" "genus" "domain" "phylum")

for dbname in "${databases[@]}"; do
    for rank in "${ranks[@]}"; do
        # All ranks use consistent naming: ${dbname}_${rank}_merged_lineage.csv
        merged_file="$summary_dir/${dbname}_${rank}_merged_lineage.csv"
        
        # Only proceed if merged file exists and is valid (not empty, not recently modified)
        if file_is_valid "$merged_file" 5; then
            # Determine pattern based on database and rank (must match script 3 patterns)
            # Also create alternate pattern for double underscores
            case "$rank" in
                "species")
                    if [ "$dbname" = "soil_microbe_db" ]; then
                        pattern="*${dbname}_filtered.b2"
                        pattern_alt="*${dbname}__filtered.b2"
                    elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
                        pattern="*gtdb_207_unfiltered*filtered.b2"
                        pattern_alt="*gtdb_207_unfiltered__filtered.b2"
                    elif [ "$dbname" = "gtdb_207" ]; then
                        pattern="*gtdb_207__filtered.b2"
                        pattern_alt="*gtdb_207_filtered.b2"
                    elif [ "$dbname" = "pluspf" ]; then
                        pattern="*pluspf__filtered.b2"
                        pattern_alt="*pluspf_filtered.b2"
                    else
                        pattern="*${dbname}_filtered.b2"
                        pattern_alt="*${dbname}__filtered.b2"
                    fi
                    ;;
                "genus")
                    if [ "$dbname" = "soil_microbe_db" ]; then
                        pattern="*${dbname}_genus_filtered.b2"
                        pattern_alt="*${dbname}__genus_filtered.b2"
                    elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
                        pattern="*gtdb_207_unfiltered_genus_filtered.b2"
                        pattern_alt="*gtdb_207_unfiltered__genus_filtered.b2"
                    elif [ "$dbname" = "gtdb_207" ]; then
                        pattern="*gtdb_207__genus_filtered.b2"
                        pattern_alt="*gtdb_207_genus_filtered.b2"
                    elif [ "$dbname" = "pluspf" ]; then
                        pattern="*pluspf__genus_filtered.b2"
                        pattern_alt="*pluspf_genus_filtered.b2"
                    else
                        pattern="*${dbname}_genus_filtered.b2"
                        pattern_alt="*${dbname}__genus_filtered.b2"
                    fi
                    ;;
                "domain")
                    if [ "$dbname" = "soil_microbe_db" ]; then
                        pattern="*${dbname}_domain_filtered.b2"
                        pattern_alt="*${dbname}__domain_filtered.b2"
                    elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
                        pattern="*gtdb_207_unfiltered_domain_filtered.b2"
                        pattern_alt="*gtdb_207_unfiltered__domain_filtered.b2"
                    elif [ "$dbname" = "gtdb_207" ]; then
                        pattern="*gtdb_207__domain_filtered.b2"
                        pattern_alt="*gtdb_207_domain_filtered.b2"
                    elif [ "$dbname" = "pluspf" ]; then
                        pattern="*pluspf__domain_filtered.b2"
                        pattern_alt="*pluspf_domain_filtered.b2"
                    else
                        pattern="*${dbname}_domain_filtered.b2"
                        pattern_alt="*${dbname}__domain_filtered.b2"
                    fi
                    ;;
                "phylum")
                    if [ "$dbname" = "soil_microbe_db" ]; then
                        pattern="*${dbname}_phylum_filtered.b2"
                        pattern_alt="*${dbname}__phylum_filtered.b2"
                    elif [ "$dbname" = "gtdb_207_unfiltered" ]; then
                        pattern="*gtdb_207_unfiltered_phylum_filtered.b2"
                        pattern_alt="*gtdb_207_unfiltered__phylum_filtered.b2"
                    elif [ "$dbname" = "gtdb_207" ]; then
                        pattern="*gtdb_207__phylum_filtered.b2"
                        pattern_alt="*gtdb_207_phylum_filtered.b2"
                    elif [ "$dbname" = "pluspf" ]; then
                        pattern="*pluspf__phylum_filtered.b2"
                        pattern_alt="*pluspf_phylum_filtered.b2"
                    else
                        pattern="*${dbname}_phylum_filtered.b2"
                        pattern_alt="*${dbname}__phylum_filtered.b2"
                    fi
                    ;;
            esac
            
            # Find and delete matching .b2 files (check both single and double underscore patterns)
            # Check both local and HARDDRIVE
            b2_files=""
            if [ -d "$architeuthis_dir_path" ]; then
                if [ -n "$pattern_alt" ] && [ "$pattern_alt" != "$pattern" ]; then
                    b2_files="$b2_files $(find "$architeuthis_dir_path" \( -name "$pattern" -o -name "$pattern_alt" \) -type f 2>/dev/null)"
                else
                    b2_files="$b2_files $(find "$architeuthis_dir_path" -name "$pattern" -type f 2>/dev/null)"
                fi
            fi
            if [ -d "$architeuthis_dir_harddrive" ]; then
                if [ -n "$pattern_alt" ] && [ "$pattern_alt" != "$pattern" ]; then
                    b2_files="$b2_files $(find "$architeuthis_dir_harddrive" \( -name "$pattern" -o -name "$pattern_alt" \) -type f 2>/dev/null)"
                else
                    b2_files="$b2_files $(find "$architeuthis_dir_harddrive" -name "$pattern" -type f 2>/dev/null)"
                fi
            fi
            
            if [ -n "$b2_files" ]; then
                # Get merged file modification time
                merged_file_time=$(stat -f %m "$merged_file" 2>/dev/null || stat -c %Y "$merged_file" 2>/dev/null || echo 0)
                
                # First pass: check if ANY .b2 file is newer than merged file
                has_newer=false
                for b2_file in $b2_files; do
                    [ -z "$b2_file" ] && continue
                    b2_file_time=$(stat -f %m "$b2_file" 2>/dev/null || stat -c %Y "$b2_file" 2>/dev/null || echo 0)
                    if [ "$b2_file_time" -gt "$merged_file_time" ]; then
                        has_newer=true
                        break
                    fi
                done
                
                # Only delete .b2 files if ALL are older than merged file
                # If any are newer, keep ALL files so script 03 can merge everything together
                if [ "$has_newer" = false ]; then
                    # All .b2 files are older than merged file, safe to delete all
                    file_count=0
                    for b2_file in $b2_files; do
                        [ -z "$b2_file" ] && continue
                        file_size=$(stat -f%z "$b2_file" 2>/dev/null || stat -c%s "$b2_file" 2>/dev/null || echo 0)
                        rm -f "$b2_file"
                        deleted_count=$((deleted_count + 1))
                        total_size_freed=$((total_size_freed + file_size))
                        file_count=$((file_count + 1))
                    done
                    if [ $file_count -gt 0 ]; then
                        echo "  ✓ Deleted $file_count $rank .b2 file(s) for $dbname (all files older than merged CSV)"
                    fi
                else
                    # Some .b2 files are newer, keep ALL files for re-merge
                    newer_count=0
                    older_count=0
                    for b2_file in $b2_files; do
                        [ -z "$b2_file" ] && continue
                        b2_file_time=$(stat -f %m "$b2_file" 2>/dev/null || stat -c %Y "$b2_file" 2>/dev/null || echo 0)
                        if [ "$b2_file_time" -gt "$merged_file_time" ]; then
                            newer_count=$((newer_count + 1))
                        else
                            older_count=$((older_count + 1))
                        fi
                    done
                    echo "  ⊘ Kept all $((newer_count + older_count)) .b2 file(s) for $dbname ($newer_count newer, $older_count older) - will re-merge all together"
                fi
            fi
        else
            echo "  ⊘ Skipped $rank .b2 files for $dbname (merged CSV not found)"
        fi
    done
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
