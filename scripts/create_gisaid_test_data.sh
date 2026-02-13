#!/bin/bash
set -e

SOURCE_DIR="/home4/lm305z/IAV_DB/flu_vgtk_integrations/tmp/gisaid-data"
DEST_DIR="/home3/oml4h/RABV-gTK/test_data/gisaid_data"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

echo "Sampling data from $SOURCE_DIR to $DEST_DIR"

# 1. Metadata: Header + 100 first rows
echo "Extracting metadata..."
head -n 101 "$SOURCE_DIR/metadata.tsv" > "$DEST_DIR/metadata.tsv"

# 2. Extract IDs (Column 44: Segment_Id)
# Skip the header line for ID extraction
tail -n +2 "$DEST_DIR/metadata.tsv" | awk -F'\t' '{print $44}' > "$DEST_DIR/ids.txt"

# 3. Extract Sequences
echo "Extracting sequences matching IDs..."
if command -v seqkit &> /dev/null; then
    seqkit grep -f "$DEST_DIR/ids.txt" "$SOURCE_DIR/all_nuc.fas" > "$DEST_DIR/sequences.fasta"
else
    echo "seqkit not found, falling back to python..."
    python -c "
from Bio import SeqIO
import sys

ids = set(line.strip() for line in open('$DEST_DIR/ids.txt'))
records = []
for record in SeqIO.parse('$SOURCE_DIR/all_nuc.fas', 'fasta'):
    if record.id in ids:
        records.append(record)
SeqIO.write(records, '$DEST_DIR/sequences.fasta', 'fasta')
"
fi

# 4. Clean up
rm "$DEST_DIR/ids.txt"

echo "Done."
echo "Metadata entries: $(wc -l < "$DEST_DIR/metadata.tsv")"
echo "Sequence entries: $(grep -c "^>" "$DEST_DIR/sequences.fasta")"
