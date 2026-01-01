#!/usr/bin/env python3
"""
16: Extract specific reads from FASTQ files and convert to FASTA for BLAST
Automatically processes both misclassification types from script 15

Usage: python3 scripts/summarize_outputs/16_extract_reads_to_fasta.py <sampleID> <R1.fastq.gz> <R2.fastq.gz> [output_dir]

Input:  
  - sampleID: Sample identifier (e.g., ORNL_046-O-20170621-COMP)
  - R1.fastq.gz, R2.fastq.gz: Original FASTQ files
  - output_dir: Optional output directory (default: data/classification/analysis_files)

Output: 
  - {sampleID}_homo_sapiens_R1.fasta, {sampleID}_homo_sapiens_R2.fasta (if PlusPF misclassifications found)
  - {sampleID}_smdb_fungal_gtdb_bacteria_R1.fasta, {sampleID}_smdb_fungal_gtdb_bacteria_R2.fasta (if SMD fungal mismatches found)

Workflow:
  1. Run script 15 → generates both read ID lists
  2. Run this script (16) → automatically extracts reads from both lists to FASTA format
  3. Run script 17 → automatically BLASTs all FASTA files
  4. Run script 18 → automatically analyzes all BLAST results
"""

import sys
import gzip
import os

def read_read_ids(read_ids_file):
    """Read read IDs from file"""
    with open(read_ids_file, 'r') as f:
        read_ids = set(line.strip() for line in f if line.strip())
    
    # Also add variants with /1 and /2 suffixes
    read_id_set = set(read_ids)
    for rid in read_ids:
        read_id_set.add(f"{rid}/1")
        read_id_set.add(f"{rid}/2")
    
    return read_id_set, read_ids

def extract_reads(fastq_file, read_id_set, output_fasta):
    """Extract reads from FASTQ and write to FASTA"""
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'
    
    extracted_count = 0
    with open_func(fastq_file, mode) as f_in, open(output_fasta, 'w') as f_out:
        line_num = 0
        current_read_id = None
        current_seq = None
        
        for line in f_in:
            line_num += 1
            line = line.rstrip('\n\r')
            
            if line_num % 4 == 1:
                current_read_id = line[1:].split()[0]
                base_read_id = current_read_id.split('/')[0]
            elif line_num % 4 == 2:
                current_seq = line
            elif line_num % 4 == 0:
                if current_read_id in read_id_set or base_read_id in read_id_set:
                    f_out.write(f">{current_read_id}\n")
                    f_out.write(f"{current_seq}\n")
                    extracted_count += 1
                
                current_read_id = None
                current_seq = None
    
    return extracted_count

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 scripts/summarize_outputs/16_extract_reads_to_fasta.py <sampleID> <R1.fastq.gz> <R2.fastq.gz> [output_dir]")
        sys.exit(1)
    
    sampleID = sys.argv[1]
    r1_file = sys.argv[2]
    r2_file = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "data/classification/analysis_files"
    
    # Find all read ID files from script 15
    read_id_files = {
        'homo_sapiens': os.path.join(output_dir, f"{sampleID}_homo_sapiens_read_ids.txt"),
        'fungal_bacteria': os.path.join(output_dir, f"{sampleID}_smdb_fungal_gtdb_bacteria_read_ids.txt"),
        'smdb_confident_fungi': os.path.join(output_dir, f"{sampleID}_smdb_confident_fungi_read_ids.txt")
    }
    
    found_files = []
    for key, read_ids_file in read_id_files.items():
        if os.path.exists(read_ids_file):
            found_files.append((key, read_ids_file))
        else:
            print(f"Note: {read_ids_file} not found, skipping {key} misclassifications")
    
    if not found_files:
        print(f"Error: No read ID files found for sample {sampleID}")
        print(f"Expected files:")
        for key, path in read_id_files.items():
            print(f"  - {path}")
        sys.exit(1)
    
    # Process each read ID file
    for key, read_ids_file in found_files:
        print(f"\nProcessing {key} misclassifications from {read_ids_file}")
        
        read_id_set, original_read_ids = read_read_ids(read_ids_file)
        print(f"  Found {len(original_read_ids)} read IDs")
        
        # Determine output prefix based on key
        if key == 'homo_sapiens':
            output_prefix = os.path.join(output_dir, f"{sampleID}_homo_sapiens")
        elif key == 'fungal_bacteria':
            output_prefix = os.path.join(output_dir, f"{sampleID}_smdb_fungal_gtdb_bacteria")
        elif key == 'smdb_confident_fungi':
            output_prefix = os.path.join(output_dir, f"{sampleID}_smdb_confident_fungi")
        else:
            output_prefix = os.path.join(output_dir, f"{sampleID}_{key}")
        
        r1_fasta = f"{output_prefix}_R1.fasta"
        n_r1 = extract_reads(r1_file, read_id_set, r1_fasta)
        print(f"  Extracted {n_r1} reads to {r1_fasta}")
        
        r2_fasta = f"{output_prefix}_R2.fasta"
        n_r2 = extract_reads(r2_file, read_id_set, r2_fasta)
        print(f"  Extracted {n_r2} reads to {r2_fasta}")
    
    print(f"\nCompleted extraction for {len(found_files)} misclassification type(s)")
