#!/usr/bin/env python3
"""
Collect all sequences that were filtered/failed during nextalign alignment.
Reads .errors.csv files from nextalign output and writes a summary of filtered sequences.
"""

import os
import csv
import argparse
from pathlib import Path


def collect_filtered_sequences(nextalign_dir: str, output_file: str) -> dict:
    """
    Scan nextalign directory for .errors.csv files and collect filtered sequence IDs.
    
    Returns dict with structure:
        {seq_name: {"reference": ref_id, "error": error_message, "warnings": warnings}}
    """
    filtered = {}
    nextalign_path = Path(nextalign_dir)
    
    # Search in both query_aln and reference_aln subdirectories
    for errors_file in nextalign_path.rglob("*.errors.csv"):
        # Get reference ID from parent directory name
        ref_id = errors_file.parent.name
        
        try:
            with open(errors_file, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    seq_name = row.get("seqName", "").strip()
                    errors = row.get("errors", "").strip()
                    warnings = row.get("warnings", "").strip()
                    
                    # Only include if there's an actual error (not just warnings)
                    if seq_name and errors:
                        filtered[seq_name] = {
                            "reference": ref_id,
                            "error": errors,
                            "warnings": warnings,
                        }
        except Exception as e:
            print(f"Warning: Could not read {errors_file}: {e}")
    
    return filtered


def write_filtered_list(filtered: dict, output_file: str):
    """Write filtered sequences to TSV file."""
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else ".", exist_ok=True)
    
    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["seq_name", "reference", "error", "warnings"])
        for seq_name, info in sorted(filtered.items()):
            writer.writerow([
                seq_name,
                info["reference"],
                info["error"],
                info["warnings"],
            ])


def write_filtered_ids_only(filtered: dict, output_file: str):
    """Write just the filtered sequence IDs, one per line (for easy exclusion)."""
    ids_file = output_file.replace(".tsv", "_ids.txt")
    if ids_file == output_file:
        ids_file = output_file + ".ids.txt"
    
    with open(ids_file, "w", encoding="utf-8") as f:
        for seq_name in sorted(filtered.keys()):
            f.write(seq_name + "\n")
    
    return ids_file


def main():
    parser = argparse.ArgumentParser(
        description="Collect sequences filtered by nextalign due to errors"
    )
    parser.add_argument(
        "-n", "--nextalign_dir",
        required=True,
        help="Path to Nextalign output directory containing query_aln/ and reference_aln/"
    )
    parser.add_argument(
        "-o", "--output",
        default="filtered_sequences.tsv",
        help="Output TSV file with filtered sequences and reasons"
    )
    parser.add_argument(
        "-b", "--base_dir",
        default=".",
        help="Base directory for output"
    )
    
    args = parser.parse_args()
    
    output_path = os.path.join(args.base_dir, args.output)
    
    print(f"Scanning {args.nextalign_dir} for filtered sequences...")
    filtered = collect_filtered_sequences(args.nextalign_dir, output_path)
    
    print(f"Found {len(filtered)} filtered sequences")
    
    if filtered:
        write_filtered_list(filtered, output_path)
        ids_file = write_filtered_ids_only(filtered, output_path)
        print(f"Wrote filtered sequences to: {output_path}")
        print(f"Wrote filtered IDs to: {ids_file}")
        
        # Print summary
        print("\nFiltered sequences:")
        for seq_name, info in sorted(filtered.items()):
            # Truncate long error messages for display
            error_short = info["error"][:80] + "..." if len(info["error"]) > 80 else info["error"]
            print(f"  {seq_name}: {error_short}")
    else:
        # Write empty files so downstream processes have something to read
        write_filtered_list({}, output_path)
        write_filtered_ids_only({}, output_path)
        print("No filtered sequences found")


if __name__ == "__main__":
    main()
