import os
import csv
from os.path import join, dirname, abspath
from argparse import ArgumentParser
from datetime import datetime
from collections import OrderedDict
import pandas as pd

"""
Transform an external TSV + FASTA into the gB_matrix-compatible format and append rows to the base gB_matrix.
- Key in input TSV: configurable with --key, default 'Segment_Id'
- Copies the key into: primary_accession, accession_version, locus columns
- Adds sequences from FASTA (header must match the key exactly)
- If a TSV row has no matching FASTA, sets exclusion = "missing fasta sequence"
- Duplicates in input by key: keeps the first, writes the others to merge_matrix_duplicates.tsv and log count
- Column mapping file: two-column TSV (source_col<TAB>target_col); mapped columns are copied/renamed
- Merge strategy: append-only (no existence check in base); writes a new merged file
"""

class NormalizeAndMerge:
    def __init__(self, gb_matrix, input_tsv, input_fasta, mapping_file, output, key="Segment_Id", dataset_source="external", prefer="input", drop_unmapped=False, log_file=None):
        self.gb_matrix = gb_matrix
        self.input_tsv = input_tsv
        self.input_fasta = input_fasta
        self.mapping_file = mapping_file
        self.output = output
        self.key = key
        self.dataset_source = dataset_source
        self.prefer = prefer                 # unused in append-only
        self.drop_unmapped = drop_unmapped

        out_dir = dirname(abspath(self.output)) or "."
        os.makedirs(out_dir, exist_ok=True)

        self.log_file = log_file or abspath("gisaid_processing.log")

    def timestamp(self):
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(self, message):
        with open(self.log_file, "a") as logf:
            logf.write(f"[{self.timestamp()}] {message}\n")
        print(f"[{self.timestamp()}] {message}")

    def read_mapping(self):
        """Read source_col<TAB>target_col pairs."""
        mapping = OrderedDict()
        with open(self.mapping_file, mode='r') as mf:
            reader = csv.reader(mf, delimiter='\t')
            for row in reader:
                if not row or len(row) < 2:
                    continue
                src, dst = row[0].strip(), row[1].strip()
                if src and dst:
                    mapping[src] = dst
        return mapping

    def parse_fasta(self):
        """Parse FASTA into sequence dict"""
        seqs = {}
        if not self.input_fasta:
            return seqs

        current_header = None
        chunks = []

        with open(self.input_fasta, 'r') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header is not None:
                        seqs[current_header] = ''.join(chunks)
                    current_header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if current_header is not None:
                seqs[current_header] = ''.join(chunks)

        return seqs

    def count_real_length(self, seq):
        if not isinstance(seq, str):
            return 0
        seq_lower = seq.lower()
        return sum(seq_lower.count(n) for n in ("a", "t", "g", "c"))
    
    def read_base_matrix_df(self) -> pd.DataFrame:
        base_df = pd.read_csv(self.gb_matrix, sep="\t", dtype=str)
        return base_df

    def read_input_tsv_df(self) -> pd.DataFrame:
        in_df = pd.read_csv(self.input_tsv, sep="\t", dtype=str)
        return in_df

    def write_duplicates_report_df(self, dup_df: pd.DataFrame):
        if dup_df is None or dup_df.empty:
            return None
        dup_file = join(dirname(abspath(self.output)) or ".", 'merge_matrix_duplicates.tsv')
        dup_df.to_csv(dup_file, sep="\t", index=False)
        return dup_file



    def process(self):
        base_df = self.read_base_matrix_df()
        self.log(f"[INFO] Loaded base matrix: {self.gb_matrix} (rows={len(base_df)})")
        base_df["dataset_source"] = "genbank" # define source name
        mapping = self.read_mapping()
        in_df = self.read_input_tsv_df()
        self.log(f"[INFO] Loaded input TSV: {self.input_tsv} (rows={len(in_df)})")
        self.log(f"[INFO] Loaded mapping file: {self.mapping_file} (pairs={len(mapping)})")
        fasta_map = self.parse_fasta()
        self.log(f"[INFO] Loaded FASTA: {self.input_fasta} (entries={len(fasta_map)})")

        if self.key not in in_df.columns:
            self.log(f"[WARN] Key column '{self.key}' not found in input TSV. No rows will be appended.")
            base_df.to_csv(self.output, sep="\t", index=False)
            self.log(f"[INFO] Merged matrix written: {self.output}")
            return

        # Duplicates handling in input tsv (keep first)
        key_series = in_df[self.key].fillna("").astype(str)
        missing_key_mask = key_series.eq("")
        dup_mask = in_df.duplicated(subset=[self.key], keep="first")

        dup_report = pd.DataFrame()
        if missing_key_mask.any():
            r1 = in_df[missing_key_mask].copy()
            r1.insert(0, "row_key", "")
            r1.insert(0, "reason", "missing-key")
            dup_report = pd.concat([dup_report, r1], ignore_index=True)
        if dup_mask.any():
            r2 = in_df[dup_mask].copy()
            r2.insert(0, "row_key", r2[self.key])
            r2.insert(0, "reason", "duplicate")
            dup_report = pd.concat([dup_report, r2], ignore_index=True)

        if not dup_report.empty:
            dup_file = self.write_duplicates_report_df(dup_report)
            self.log(f"[WARN] Duplicates/missing-key in input: {len(dup_report)} (report: {dup_file})")
        else:
            self.log("[INFO] No duplicates in input.")

        keep_mask = (~missing_key_mask) & (~dup_mask)
        norm_df = in_df.loc[keep_mask].copy()

        if mapping:
            norm_df = norm_df.rename(columns=mapping)

        # copy key to gB fields
        norm_df["primary_accession"] = norm_df[self.key]
        norm_df["accession_version"] = norm_df[self.key]
        norm_df["locus"] = norm_df[self.key]
        norm_df["gi_number"] = norm_df[self.key]

        # import sequences from FASTA to column 'sequence'
        norm_df["sequence"] = norm_df[self.key].map(fasta_map).fillna("")
        norm_df.loc[norm_df["sequence"] == "", "exclusion"] = "missing fasta sequence"
        norm_df["real_length"] = norm_df["sequence"].apply(self.count_real_length)
        
        norm_df["dataset_source"] = self.dataset_source

        base_cols = list(base_df.columns)

        if self.drop_unmapped:
            mapping_targets = set(mapping.values())
            required_cols = {'primary_accession', 'accession_version', 'locus', 'sequence', 'dataset_source', 'exclusion'}
            keep_cols = list(set(base_cols) | mapping_targets | required_cols | {self.key})
            ordered_keep = [c for c in base_cols if c in keep_cols] + [c for c in norm_df.columns if c in keep_cols and c not in base_cols]
            norm_df = norm_df[[c for c in ordered_keep if c in norm_df.columns]]
            extras = [c for c in norm_df.columns if c not in base_cols]
            final_cols = base_cols + extras
        else:
            extras = [c for c in norm_df.columns if c not in base_cols]
            final_cols = base_cols + extras

        norm_df = norm_df.reindex(columns=final_cols, fill_value="")
        merged_df = pd.concat([base_df.reindex(columns=final_cols, fill_value=""), norm_df], ignore_index=True)

        merged_df.to_csv(self.output, sep="\t", index=False)
        self.log(f"[INFO] Merged matrix written: {self.output}")


if __name__ == "__main__":
    parser = ArgumentParser(description='Normalize external TSV+FASTA to gB_matrix format and append to base matrix.')
    parser.add_argument('-g', '--gb_matrix', help='Base gB_matrix file', default='gB_matrix.tsv')
    parser.add_argument('-t', '--input_tsv', help='Input TSV to transform', default='./tmp/gisaid-data/metadata.tsv')
    parser.add_argument('-f', '--input_fasta', help='Input FASTA with sequences', default='./tmp/gisaid-data/all_nuc.fas')
    parser.add_argument('-m', '--mapping_file', help='Column mapping file (source_col<TAB>target_col)',
                        default='./generic-influenza/column_mapping.tsv')
    parser.add_argument('-k', '--key', help='Key column name in input TSV', default='Segment_Id')
    parser.add_argument('-s', '--dataset_source', 
                        help='Source name value to annotate in the "dataset_source" column in "gB_matrix_merged.tsv"', default='external')
    parser.add_argument('-o', '--output', help='Full output path including filename',
                        default='./tmp/gisaid-data/gB_matrix_merged.tsv')
    parser.add_argument('--log_file', help='Path to main process log (default: gisaid_processing.log in CWD)', default=None)
    parser.add_argument('--prefer', choices=['input', 'base'], default='input',
                        help='(accepted for compatibility) preferred side when upsert mode is enabled')
    parser.add_argument('--drop_unmapped', action='store_true',
                        help='If set, drop columns from the new dataset that are not in base matrix, mapping targets, or required columns')

    args = parser.parse_args()

    job = NormalizeAndMerge(
        gb_matrix=args.gb_matrix,
        input_tsv=args.input_tsv,
        input_fasta=args.input_fasta,
        mapping_file=args.mapping_file,
        output=args.output,
        key=args.key,
        dataset_source=args.dataset_source,
        prefer=args.prefer,
        drop_unmapped=args.drop_unmapped,
        log_file=args.log_file
    )
    job.process()
