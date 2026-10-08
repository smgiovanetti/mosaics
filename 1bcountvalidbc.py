import os
import pandas as pd
from Bio.Seq import Seq
import re

def extract_barcode(seq):
    match = re.search(r'CAGCCCTACAGGGT(.{15})', seq)
    if match:
        # Get the 15bp sequence after motif
        extracted = match.group(1)
        # Take reverse complement
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rc_extracted = ''.join(complement[base] for base in extracted[::-1])
        # Take first 13bp of reverse complement as the designed barcode
        return rc_extracted[:13] if len(rc_extracted) >= 13 else None
    return None

def check_barcode_pattern(barcode):
    """Check if barcode matches NGANAGNACNACN pattern"""
    if not barcode or len(barcode) != 13:
        return False
    
    # Pattern: NGANAGNACNACN
    # Fixed positions and their expected bases:
    expected_fixed = {
        1: 'G',   # position 1 = G
        2: 'A',   # position 2 = A  
        4: 'A',   # position 4 = A
        5: 'G',   # position 5 = G
        7: 'A',   # position 7 = A
        8: 'C',   # position 8 = C
        10: 'A',  # position 10 = A
        11: 'C'   # position 11 = C
    }
    
    for pos, expected_base in expected_fixed.items():
        if barcode[pos] != expected_base:
            return False
    
    return True

def classify_sequence(nt_seq):
    # Trim to multiple of 3 to avoid partial codon warning
    trimmed = nt_seq[:len(nt_seq) - (len(nt_seq) % 3)] if len(nt_seq) % 3 != 0 else nt_seq
    
    try:
        aa_seq = str(Seq(trimmed).translate())
    except Exception:
        return {pos: 'invalid_nt' for pos in positions}

    calls = {}
    for pos, aa_defs in positions.items():
        idx = pos - codon_start
        if idx < len(aa_seq):
            aa = aa_seq[idx]
            if aa == aa_defs['R']:
                calls[pos] = 'R'
            elif aa == aa_defs['S']:
                calls[pos] = 'S'
            else:
                calls[pos] = 'other'
        else:
            calls[pos] = 'truncated'
    return calls

# Define RaTG13 (R) and SARS-CoV-2 (S) amino acids at each position
positions = {
    439: {'R': 'K', 'S': 'N'},
    440: {'R': 'H', 'S': 'N'},
    441: {'R': 'I', 'S': 'L'},
    443: {'R': 'A', 'S': 'S'},
    445: {'R': 'E', 'S': 'V'},
    449: {'R': 'F', 'S': 'Y'},
    459: {'R': 'A', 'S': 'S'},
    478: {'R': 'K', 'S': 'T'},
    483: {'R': 'Q', 'S': 'V'},
    484: {'R': 'T', 'S': 'E'},
    486: {'R': 'L', 'S': 'F'},
    490: {'R': 'Y', 'S': 'F'},
    493: {'R': 'Y', 'S': 'Q'},
    494: {'R': 'R', 'S': 'S'},
    498: {'R': 'Y', 'S': 'Q'},
    501: {'R': 'D', 'S': 'N'},
    505: {'R': 'H', 'S': 'Y'}
}

codon_start = 439

# Directories
input_dir = 'frequency_tables_barcode'
output_dir = 'valid_barcode_tables'
os.makedirs(output_dir, exist_ok=True)

# Process only specific files
target_files = [
    'm1-0_S35_R1_001.fastq.gz_frequency.txt',
    'w1-0_S33_R1_001.fastq.gz_frequency.txt'
]

summary_stats = []

# Process only the target files
for filename in target_files:
    if not os.path.exists(os.path.join(input_dir, filename)):
        print(f"File not found: {filename}")
        continue
    
    input_path = os.path.join(input_dir, filename)
    df = pd.read_csv(input_path, sep="\t")

    if 'Sequence' not in df.columns:
        print(f"Skipping {filename}: no 'Sequence' column found.")
        continue

    print(f"\nProcessing {filename}...")
    
    # Extract barcodes
    df['barcode'] = df['Sequence'].apply(extract_barcode)
    
    # Filter to sequences with valid barcodes matching the pattern
    df['valid_barcode'] = df['barcode'].apply(check_barcode_pattern)
    
    # Stats before filtering
    total_sequences = len(df)
    total_with_barcode = df['barcode'].notna().sum()
    valid_barcode_sequences = df['valid_barcode'].sum()
    
    # Get frequency column
    freq_col = [col for col in df.columns if col.endswith('_freq')][0]
    total_reads = df[freq_col].sum()
    reads_with_barcode = df[df['barcode'].notna()][freq_col].sum()
    reads_valid_barcode = df[df['valid_barcode']][freq_col].sum()
    
    print(f"  Total sequences: {total_sequences:,}")
    print(f"  Sequences with barcodes: {total_with_barcode:,} ({100*total_with_barcode/total_sequences:.1f}%)")
    print(f"  Sequences with valid barcodes: {valid_barcode_sequences:,} ({100*valid_barcode_sequences/total_sequences:.1f}%)")
    print(f"  Reads in valid barcodes: {reads_valid_barcode:,} / {total_reads:,} ({100*reads_valid_barcode/total_reads:.1f}%)")
    
    # Filter to only valid barcode sequences
    filtered_df = df[df['valid_barcode']].copy()
    
    if len(filtered_df) == 0:
        print(f"  WARNING: No valid barcodes found in {filename}")
        continue
    
    # Drop the helper columns
    filtered_df = filtered_df.drop(columns=['valid_barcode'])
    
    # Add position classifications
    annotations = filtered_df['Sequence'].apply(classify_sequence).apply(pd.Series)
    annotated_df = pd.concat([filtered_df, annotations], axis=1)
    
    # Get frequency columns
    freq_columns = [col for col in annotated_df.columns if col.endswith('_freq')]
    
    # Group by barcode and aggregate
    agg_dict = {
        'Sequence': 'first',  # Take the most frequent sequence for each barcode
        **{col: 'sum' for col in freq_columns},
        **{pos: 'first' for pos in positions.keys()}
    }
    
    # Sort by frequency and group by barcode
    annotated_df = annotated_df.sort_values(freq_columns[0], ascending=False)
    final_df = annotated_df.groupby('barcode').agg(agg_dict).reset_index()
    
    # Reorder columns
    final_columns = ['Sequence', 'barcode'] + freq_columns + list(positions.keys())
    final_df = final_df[final_columns]
    
    # Save filtered results
    output_path = os.path.join(output_dir, filename)
    final_df.to_csv(output_path, sep="\t", index=False)
    
    print(f"  Saved {len(final_df):,} valid barcodes to {output_path}")
    
    # Store summary stats
    summary_stats.append({
        'filename': filename,
        'total_sequences': total_sequences,
        'valid_barcode_sequences': valid_barcode_sequences,
        'valid_barcode_pct': 100*valid_barcode_sequences/total_sequences,
        'total_reads': total_reads,
        'valid_reads': reads_valid_barcode,
        'valid_reads_pct': 100*reads_valid_barcode/total_reads,
        'unique_valid_barcodes': len(final_df)
    })

# Save summary
summary_df = pd.DataFrame(summary_stats)
summary_path = os.path.join(output_dir, '_filtering_summary.txt')
summary_df.to_csv(summary_path, sep='\t', index=False)

print(f"\n=== FILTERING SUMMARY ===")
print(f"Processed {len(summary_stats)} files")
print(f"Average valid barcode rate: {summary_df['valid_barcode_pct'].mean():.1f}%")
print(f"Average valid read rate: {summary_df['valid_reads_pct'].mean():.1f}%")
print(f"Total unique valid barcodes across all files: {summary_df['unique_valid_barcodes'].sum():,}")
print(f"Summary saved to: {summary_path}")