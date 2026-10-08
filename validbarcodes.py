import re

def check_barcode_pattern(barcode):
    pattern = r'^.GA.AG.AC.AC.$'  # N=any nucleotide (.), fixed letters
    return bool(re.match(pattern, barcode))

# Check how many of your barcodes match the expected pattern
valid_count = 0
invalid_examples = []

with open("common_barcodes.txt", "r") as f:
    for line in f:
        barcode = line.strip()
        if check_barcode_pattern(barcode):
            valid_count += 1
        else:
            invalid_examples.append(barcode)

print(f"Valid pattern barcodes: {valid_count}")
print(f"Invalid pattern barcodes: {len(invalid_examples)}")
print("Examples of invalid barcodes:")
for ex in invalid_examples[:10]:
    print(f"  {ex}")