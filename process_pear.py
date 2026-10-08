#!/usr/bin/env python3
import os
import subprocess
import sys
import regex
from Bio import SeqIO
from collections import Counter

def process_file(filename):
    """Process a single FASTQ file"""
    
    # Anchor sequences
    ANCHOR_1 = "GATGTGTAATCGCCTGGAACTCC"
    ANCHOR_2 = "GAGCTTTTGAATGCCCCG"
    MAX_MISMATCHES = 1

    # Directories
    input_dir = "pear_output"
    output_dir = "frequency_tables_barcode"
    os.makedirs(output_dir, exist_ok=True)

    input_path = os.path.join(input_dir, filename)
    base_name = filename.split(".assembled.fastq")[0]
    output_path = os.path.join(output_dir, f"{base_name}_frequency.txt")

    sequence_counts = Counter()
    total_reads = 0
    matched_reads = 0

    # Pre-compile fuzzy patterns
    pat1 = regex.compile(r"(?:" + ANCHOR_1 + r"){e<=" + str(MAX_MISMATCHES) + r"}")
    pat2 = regex.compile(r"(?:" + ANCHOR_2 + r"){e<=" + str(MAX_MISMATCHES) + r"}")

    print(f"Processing {filename}...")

    with open(input_path, "r") as f:
        for record in SeqIO.parse(f, "fastq"):
            total_reads += 1
            full_seq = str(record.seq).upper()

            m1 = pat1.search(full_seq)
            m2 = pat2.search(full_seq)

            if m1 and m2 and m2.start() > m1.end():
                trimmed_seq = full_seq[m1.end():m2.start()]
                sequence_counts[trimmed_seq] += 1
                matched_reads += 1

    match_rate = (matched_reads / total_reads * 100) if total_reads > 0 else 0
    unique_seqs = len(sequence_counts)
    total_frequency = sum(sequence_counts.values())

    # Write frequency table, sorted by most common
    with open(output_path, "w") as out:
        out.write(f"Sequence\t{base_name}_freq\t{base_name}_ratio\n")
        for seq, freq in sequence_counts.most_common():
            ratio = freq / total_frequency
            out.write(f"{seq}\t{freq}\t{ratio:.6f}\n")

    # Write individual summary
    summary_path = os.path.join(output_dir, f"{base_name}_summary.txt")
    with open(summary_path, "w") as s:
        s.write(f"filename\ttotal_reads\tmatched_reads\tmatch_rate\tunique_sequences\n")
        s.write(f"{base_name}\t{total_reads}\t{matched_reads}\t{match_rate:.1f}%\t{unique_seqs}\n")

    print(f"Completed {filename}: {matched_reads:,}/{total_reads:,} reads ({match_rate:.1f}%)")
    return base_name, total_reads, matched_reads, match_rate, unique_seqs

if __name__ == "__main__":
    # Install regex if needed
    subprocess.check_call([sys.executable, "-m", "pip", "install", "regex", "-q"])
    
    # Get filename from command line argument
    if len(sys.argv) != 2:
        print("Usage: python process_single_file.py <filename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    process_file(filename)