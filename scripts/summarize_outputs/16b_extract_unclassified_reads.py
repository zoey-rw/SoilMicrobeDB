#!/usr/bin/env python3
"""
16b: Extract unclassified reads (not in any scores file, pre-scoring)
Samples reads that were not classified by any database for BLAST verification

Usage: python3 scripts/summarize_outputs/16b_extract_unclassified_reads.py <sampleID> <R1.fastq.gz> <R2.fastq.gz> [sample_size] [output_dir]

Input:  
  - sampleID: Sample identifier (e.g., ORNL_046-O-20170621-COMP)
  - R1.fastq.gz, R2.fastq.gz: Original FASTQ files
  - sample_size: Number of unclassified reads to sample (default: 500)
  - output_dir: Optional output directory (default: data/classification/analysis_files)

Output: 
  - {sampleID}_unclassified_read_ids.txt (sampled read IDs)
  - {sampleID}_unclassified_R1.fasta (extracted sequences for BLAST)
"""

import sys
import gzip
import os
import random

def load_classified_read_ids(classified_ids_file):
    """Load classified read IDs from file"""
    classified = set()
    if os.path.exists(classified_ids_file):
        with open(classified_ids_file, 'r') as f:
            for line in f:
                read_id = line.strip()
                if read_id:
                    classified.add(read_id)
                    # Also add variants with /1 and /2 suffixes
                    classified.add(f"{read_id}/1")
                    classified.add(f"{read_id}/2")
    return classified

def extract_unclassified_reads(fastq_file, classified_set, sample_size, output_fasta):
    """Extract unclassified reads from FASTQ and write to FASTA"""
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'
    
    unclassified_reads = []
    extracted_count = 0
    line_num = 0
    
    print(f"  Scanning {fastq_file} for unclassified reads...")
    
    with open_func(fastq_file, mode) as f_in:
        current_read_id = None
        current_seq = None
        
        for line in f_in:
            line_num += 1
            line = line.rstrip('\n\r')
            
            if line_num % 4 == 1:
                # Read ID line
                current_read_id = line[1:].split()[0]
                base_read_id = current_read_id.split('/')[0]
            elif line_num % 4 == 2:
                # Sequence line
                current_seq = line
            elif line_num % 4 == 0:
                # Quality line - end of read
                # Check if this read is unclassified
                if current_read_id not in classified_set and base_read_id not in classified_set:
                    unclassified_reads.append((current_read_id, current_seq))
                
                current_read_id = None
                current_seq = None
                
                # Progress update
                if line_num % 1000000 == 0:
                    print(f"    Processed {line_num // 4:,} reads, found {len(unclassified_reads):,} unclassified")
    
    print(f"  Found {len(unclassified_reads):,} unclassified reads")
    
    # Sample if needed
    if len(unclassified_reads) > sample_size:
        print(f"  Sampling {sample_size} reads from {len(unclassified_reads):,} unclassified reads")
        unclassified_reads = random.sample(unclassified_reads, sample_size)
    
    # Write to FASTA
    print(f"  Writing {len(unclassified_reads)} reads to {output_fasta}")
    with open(output_fasta, 'w') as f_out:
        for read_id, seq in unclassified_reads:
            f_out.write(f">{read_id}\n")
            f_out.write(f"{seq}\n")
            extracted_count += 1
    
    return extracted_count, len(unclassified_reads)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 scripts/summarize_outputs/16b_extract_unclassified_reads.py <sampleID> <R1.fastq.gz> <R2.fastq.gz> [sample_size] [output_dir]")
        print("  sampleID: Sample identifier (e.g., ORNL_046-O-20170621-COMP)")
        print("  R1.fastq.gz, R2.fastq.gz: Original FASTQ files")
        print("  sample_size: Number of unclassified reads to sample (default: 500)")
        print("  output_dir: Optional output directory (default: data/classification/analysis_files)")
        sys.exit(1)
    
    sampleID = sys.argv[1]
    r1_file = sys.argv[2]
    r2_file = sys.argv[3] if len(sys.argv) > 3 else None
    sample_size = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    output_dir = sys.argv[5] if len(sys.argv) > 5 else "data/classification/analysis_files"
    
    print(f"=== Extracting Unclassified Reads ===")
    print(f"Sample: {sampleID}")
    print(f"Sample size: {sample_size}")
    print(f"")
    
    # Load classified read IDs
    classified_ids_file = os.path.join(output_dir, f"{sampleID}_classified_read_ids.txt")
    print(f"Loading classified read IDs from: {classified_ids_file}")
    classified_set = load_classified_read_ids(classified_ids_file)
    print(f"  Loaded {len(classified_set):,} classified read IDs (including /1 and /2 variants)")
    print(f"")
    
    # Extract unclassified reads from R1
    output_prefix = os.path.join(output_dir, f"{sampleID}_unclassified")
    r1_fasta = f"{output_prefix}_R1.fasta"
    
    print(f"Processing R1 file: {r1_file}")
    n_r1, total_unclassified_r1 = extract_unclassified_reads(r1_file, classified_set, sample_size, r1_fasta)
    print(f"  Extracted {n_r1} reads to {r1_fasta}")
    print(f"")
    
    # Save sampled read IDs
    read_ids_file = f"{output_prefix}_read_ids.txt"
    with open(read_ids_file, 'w') as f:
        # Read from FASTA to get read IDs
        with open(r1_fasta, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    read_id = line[1:].strip().split()[0]
                    base_id = read_id.split('/')[0]
                    f.write(f"{base_id}\n")
    
    print(f"Saved sampled read IDs to: {read_ids_file}")
    print(f"")
    print(f"Next steps:")
    print(f"  1. BLAST the FASTA file:")
    print(f"     python3 scripts/summarize_outputs/17_blast_reads_ncbi.py {sampleID} {sample_size} {output_dir} 4")
    print(f"     (This will automatically detect the unclassified FASTA file)")
    print(f"  2. Analyze results with script 18")





