#!/usr/bin/env python3
"""
17: BLAST reads against NCBI nt database using web API
Automatically processes all FASTA files from script 16

Usage: python3 scripts/summarize_outputs/17_blast_reads_ncbi.py <sampleID> [sample_size] [output_dir] [num_threads]

Input:  
  - sampleID: Sample identifier (e.g., ORNL_046-O-20170621-COMP)
  - sample_size: Optional number of sequences to sample per FASTA file (default: 100)
  - output_dir: Optional output directory (default: data/classification/analysis_files)
  - num_threads: Optional number of parallel threads (default: 4)

Output: 
  - {sampleID}_homo_sapiens_R1_blast.txt (if PlusPF misclassifications found)
  - {sampleID}_smdb_fungal_gtdb_bacteria_R1_blast.txt (if SMD fungal mismatches found)

Workflow:
  1. Run script 15 → generates both read ID lists
  2. Run script 16 → automatically extracts reads from both lists to FASTA format
  3. Run this script (17) → automatically BLASTs all FASTA files
  4. Run script 18 → automatically analyzes all BLAST results

Note: Uses NCBI web API with rate limiting (10 seconds between requests per thread)
      Estimated time: ~10 seconds per sequence (divided by number of threads)
      Expect threshold: 1e-4 (relaxed from 1e-5 to capture more matches)
      Appends to existing BLAST files (does not overwrite) - skips already processed sequences
      Supports multi-threading for faster processing
"""

import sys
import time
import os
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import random
import threading
from queue import Queue

NCBIWWW.email = "zoeywerbin@example.com"

def blast_sequence(sequence, sequence_id, max_attempts=3):
    """BLAST a single sequence against NCBI nt database"""
    for attempt in range(max_attempts):
        try:
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=sequence,
                hitlist_size=10,
                expect=1e-2,  # Very lenient (1e-2) to capture more matches, including those missed by 1e-4
                format_type="XML"
            )
            return result_handle
        except Exception as e:
            if attempt < max_attempts - 1:
                wait_time = (attempt + 1) * 30
                time.sleep(wait_time)
            else:
                return None

def parse_blast_xml(result_handle):
    """Parse BLAST XML results and return top hit information"""
    try:
        blast_record = NCBIXML.read(result_handle)
        
        if not blast_record.alignments:
            return None
        
        top_hit = blast_record.alignments[0]
        top_hsp = top_hit.hsps[0]
        
        # Extract taxonomy info if available
        taxid = None
        sciname = None
        if hasattr(top_hit, 'title'):
            # Try to extract from title
            title = top_hit.title
        else:
            title = str(top_hit)
        
        return {
            'accession': top_hit.accession if hasattr(top_hit, 'accession') else 'unknown',
            'title': title,
            'evalue': top_hsp.expect,
            'bitscore': top_hsp.bits,
            'identity': top_hsp.identities / top_hsp.align_length * 100 if top_hsp.align_length > 0 else 0,
            'align_length': top_hsp.align_length,
            'query_start': top_hsp.query_start,
            'query_end': top_hsp.query_end,
            'subject_start': top_hsp.sbjct_start,
            'subject_end': top_hsp.sbjct_end
        }
    except Exception:
        return None

def worker_thread(sequence_queue, output_file, file_lock, progress_counter, progress_lock, total_count):
    """Worker thread that processes sequences from the queue"""
    while True:
        item = sequence_queue.get()
        if item is None:  # Poison pill to stop thread
            sequence_queue.task_done()
            break
        
        record, seq_num = item
        result_handle = blast_sequence(str(record.seq), record.id)
        
        if result_handle:
            hit_info = parse_blast_xml(result_handle)
            result_handle.close()
            
            if hit_info:
                result_line = f"{record.id}\t{hit_info['accession']}\t{hit_info['title']}\t" \
                            f"{hit_info['evalue']}\t{hit_info['bitscore']}\t{hit_info['identity']:.2f}\t" \
                            f"{hit_info['align_length']}\t{hit_info['query_start']}\t{hit_info['query_end']}\t" \
                            f"{hit_info['subject_start']}\t{hit_info['subject_end']}\n"
            else:
                result_line = f"{record.id}\tno_hits\t\t\t\t\t\t\t\n"
        else:
            result_line = f"{record.id}\tblast_failed\t\t\t\t\t\t\t\n"
        
        # Thread-safe file writing
        with file_lock:
            with open(output_file, 'a') as out:
                out.write(result_line)
                out.flush()
        
        # Thread-safe progress reporting
        with progress_lock:
            progress_counter[0] += 1
            count = progress_counter[0]
            if count % 10 == 0 or count == total_count:
                print(f"  Progress: {count}/{total_count} sequences processed")
        
        sequence_queue.task_done()
        time.sleep(10)  # Rate limiting per thread

def blast_fasta_file(fasta_file, output_file, sample_size, num_threads=4):
    """BLAST a single FASTA file, appending to existing results if file exists"""
    print(f"\nBLASTing {fasta_file}...")
    print(f"  Using {num_threads} threads")
    
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Total sequences in FASTA: {len(sequences)}")
    
    # Check for existing results - only skip sequences with successful hits
    # Retry sequences that had "no_hits" or "blast_failed" with new threshold
    successfully_processed = set()
    file_exists = os.path.exists(output_file)
    
    if file_exists:
        print(f"  Existing BLAST file found: {output_file}")
        # Check what sequences have results (including "no_hits" which are valid results)
        with open(output_file, 'r') as f:
            f.readline()  # Skip header
            all_processed = set()
            failed_only = set()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    read_id = parts[0]
                    result = parts[1] if len(parts) > 1 else ""
                    if result == "blast_failed":
                        # Only retry actual failures, not "no_hits" (which are valid results)
                        failed_only.add(read_id)
                    else:
                        # Any other result (successful hit or "no_hits") means sequence was processed
                        all_processed.add(read_id)
        
        # If file has sufficient results (>=500), treat as complete and skip entirely
        # "no_hits" is a valid result and shouldn't be retried
        # This handles cases where FASTA has more sequences than were originally BLASTed
        if len(all_processed) >= 500:
            print(f"  File has {len(all_processed)} results (sufficient for analysis)")
            if len(failed_only) == 0:
                print(f"  All sequences processed (including 'no_hits' which are valid). Skipping.")
                return
            else:
                # Only retry the specific "blast_failed" sequences, not all missing ones
                print(f"  Will retry {len(failed_only)} sequences with 'blast_failed' results only")
                sequences_to_retry = [s for s in sequences if s.id in failed_only]
                if len(sequences_to_retry) == 0:
                    print(f"  No sequences to retry. Skipping.")
                    return
                # Skip to processing only the failed sequences
                sequences_to_process = sequences_to_retry
                successfully_processed = all_processed  # Mark all processed as done
        elif (len(all_processed) + len(failed_only)) >= len(sequences) * 0.9:
            # File has results for 90%+ of sequences (including failed) - treat as complete
            # Skip entirely, even if there are some "blast_failed" results
            # (User indicated these files are already complete)
            total_with_results = len(all_processed) + len(failed_only)
            print(f"  File has results for {total_with_results}/{len(sequences)} sequences (90%+ complete)")
            print(f"  Treating as complete. Skipping (including {len(failed_only)} 'blast_failed' results).")
            return
        else:
            # File incomplete - process missing sequences
            successfully_processed = all_processed
            print(f"  Already processed: {len(successfully_processed)} sequences")
            print(f"  Will retry: {len(failed_only)} 'blast_failed' sequences")
            # Filter out successfully processed sequences (retry failed ones and missing ones)
            sequences_to_process = [s for s in sequences if s.id not in successfully_processed]
    
    # If sequences_to_process wasn't set above (for incomplete files), set it now
    if 'sequences_to_process' not in locals():
        sequences_to_process = [s for s in sequences if s.id not in successfully_processed]
    print(f"  Sequences to process: {len(sequences_to_process)}")
    
    if len(sequences_to_process) == 0:
        print(f"  All sequences already processed. Skipping.")
        return
    
    # Sample if needed
    num_remaining = len(sequences_to_process)
    if num_remaining > sample_size:
        sequences_to_process = random.sample(sequences_to_process, sample_size)
        print(f"  Sampling {sample_size} sequences from {num_remaining} remaining")
    
    # Write header if file is new
    if not file_exists:
        with open(output_file, 'w') as out:
            out.write("read_id\taccession\ttitle\tevalue\tbitscore\tidentity\talign_length\tquery_start\tquery_end\tsubject_start\tsubject_end\n")
    
    # Create queue and threads
    sequence_queue = Queue()
    file_lock = threading.Lock()
    progress_counter = [0]  # Use list for thread-safe counter
    progress_lock = threading.Lock()
    
    # Add sequences to queue
    for i, record in enumerate(sequences_to_process, 1):
        sequence_queue.put((record, i))
    
    # Start worker threads
    threads = []
    for _ in range(num_threads):
        t = threading.Thread(target=worker_thread, 
                           args=(sequence_queue, output_file, file_lock, progress_counter, progress_lock, len(sequences_to_process)))
        t.start()
        threads.append(t)
    
    # Wait for all sequences to be processed
    sequence_queue.join()
    
    # Stop all threads
    for _ in range(num_threads):
        sequence_queue.put(None)
    for t in threads:
        t.join()
    
    print(f"  Completed: {output_file}")
    print(f"  Total sequences in file: {len(successfully_processed) + len(sequences_to_process)}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/summarize_outputs/17_blast_reads_ncbi.py <sampleID> [sample_size] [output_dir] [num_threads]")
        print("  sampleID: Sample identifier (e.g., ORNL_046-O-20170621-COMP)")
        print("  sample_size: Optional number of sequences to sample per FASTA file (default: 100)")
        print("  output_dir: Optional output directory (default: data/classification/analysis_files)")
        print("  num_threads: Optional number of parallel threads (default: 4)")
        sys.exit(1)
    
    sampleID = sys.argv[1]
    sample_size = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "data/classification/analysis_files"
    num_threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    
    # Find all R1 FASTA files from script 16 and 16b
    fasta_patterns = [
        os.path.join(output_dir, f"{sampleID}_homo_sapiens_R1.fasta"),
        os.path.join(output_dir, f"{sampleID}_smdb_fungal_gtdb_bacteria_R1.fasta"),
        os.path.join(output_dir, f"{sampleID}_smdb_confident_fungi_R1.fasta"),
        os.path.join(output_dir, f"{sampleID}_unclassified_R1.fasta")
    ]
    
    found_fastas = []
    for fasta_file in fasta_patterns:
        if os.path.exists(fasta_file):
            found_fastas.append(fasta_file)
        else:
            print(f"Note: {fasta_file} not found, skipping")
    
    if not found_fastas:
        print(f"Error: No FASTA files found for sample {sampleID}")
        print(f"Expected files:")
        for pattern in fasta_patterns:
            print(f"  - {pattern}")
        sys.exit(1)
    
    # BLAST each FASTA file
    for fasta_file in found_fastas:
        # Determine output file name
        base_name = os.path.basename(fasta_file).replace("_R1.fasta", "")
        output_file = os.path.join(output_dir, f"{base_name}_R1_blast.txt")
        
        # For unclassified reads, use a larger sample size if not specified
        file_sample_size = sample_size
        if "unclassified" in base_name and sample_size == 100:
            file_sample_size = 500  # Default to 500 for unclassified
        
        blast_fasta_file(fasta_file, output_file, file_sample_size, num_threads)
    
    print(f"\nCompleted BLAST for {len(found_fastas)} FASTA file(s)")

if __name__ == "__main__":
    main()
