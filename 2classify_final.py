import os
import pandas as pd
from Bio.Seq import Seq
import regex

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


def classify_sequence(nt_seq):
    try:
        aa_seq = str(Seq(nt_seq).translate())
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


def is_valid_barcode(bc):
    pattern = "NGTNGTNCTNTCNTT"

    if len(bc) != len(pattern):
        return False

    for b, p in zip(bc, pattern):
        if p == "N":
            if b not in "ACGT":
                return False
        else:
            if b != p:
                return False
    return True


def extract_barcode(seq):
    seq = str(seq).upper()

    match = regex.search(r'(?:CAGCCCTACAGGGT){e<=2}(.{15})', seq)
    if match:
        bc = match.group(1)
        if is_valid_barcode(bc):
            return bc

    match = regex.search(r'(?:ACCCTGTAGGGCTG){e<=2}(.{15})', seq)
    if match:
        bc = match.group(1)
        if is_valid_barcode(bc):
            return bc

    return 'NOT_FOUND'


input_dir = 'frequency_tables_barcode'
output_dir = 'annotated_tables_barcode_new'
os.makedirs(output_dir, exist_ok=True)

for filename in sorted(os.listdir(input_dir)):
    if not filename.endswith(".txt"):
        continue

    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, filename)

    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
        print(f"Skipping {filename}: output already exists")
        continue

    df = pd.read_csv(input_path, sep="\t", encoding="latin1")

    if 'Sequence' not in df.columns:
        print(f"Skipping {filename}: no 'Sequence' column found.")
        continue

    df['Sequence'] = df['Sequence'].astype(str).str.upper()
    df = df.loc[:, ~df.columns.str.endswith('_ratio')].copy()

    df['barcode'] = df['Sequence'].apply(extract_barcode)

    annotations = pd.DataFrame(df['Sequence'].apply(classify_sequence).tolist())
    annotated_df = pd.concat([df, annotations], axis=1)

    spike_columns = [col for col in df.columns if col.endswith('freq')]

    collapsed_df = annotated_df.groupby('Sequence').agg(
        {**{col: 'sum' for col in spike_columns},
         **{pos: 'first' for pos in positions.keys()},
         'barcode': 'first'}
    ).reset_index()

    final_columns = ['Sequence', 'barcode'] + spike_columns + list(positions.keys())
    collapsed_df = collapsed_df[final_columns]

    collapsed_df.to_csv(output_path, sep="\t", index=False)

    n_not_found = (df['barcode'] == 'NOT_FOUND').sum()
    print(f"{filename}: NOT_FOUND = {n_not_found:,} / {len(df):,}")
    print(f"Processed and collapsed {filename}")