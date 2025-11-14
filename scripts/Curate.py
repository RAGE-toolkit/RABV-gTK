#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
from argparse import ArgumentParser
from typing import List, Dict, Tuple


class BlastAlignment:
    KEY = "gi_number"
    CURATOR = "curator"

    def __init__(self, gb_matrix: str, curated_file: str,
                 base_dir: str = "tmp", output_dir: str = "Curated",
                 output_file: str = "gB_matrix.tsv") -> None:
        self.gb_matrix = gb_matrix
        self.curated_file = curated_file
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.output_file = output_file

    def process(self) -> None:
        # Read original matrix (preserve header order)
        gb_header, gb_rows = self._read_tsv(self.gb_matrix, "GenBank matrix")

        # Read curated file
        cur_header, cur_rows = self._read_tsv(self.curated_file, "curated file")

        # Basic checks
        self._require_column(gb_header, self.KEY, "GenBank matrix")
        self._require_column(cur_header, self.KEY, "curated file")
        self._require_column(cur_header, self.CURATOR, "curated file")

        # curator must be non-empty per curated row
        empties = [r.get(self.KEY, "") for r in cur_rows
                   if (r.get(self.CURATOR, "") or "").strip() == ""]
        if empties:
            raise ValueError(
                f"'curator' is empty for {len(empties)} curated rows. "
                f"First few gi_accession: {empties[:10]}"
            )

        # Build index: gi_accession -> row dict (first occurrence kept)
        gb_index: Dict[str, Dict[str, str]] = {}
        for row in gb_rows:
            k = (row.get(self.KEY, "") or "").strip()
            if k and k not in gb_index:
                gb_index[k] = row

        # Replace fields
        updated_rows = 0
        updated_cells = 0
        missing = 0
        ignored_cols = set()

        # columns we will try to apply (exclude key and curator)
        candidate_cols = [c for c in cur_header if c not in (self.KEY, self.CURATOR)]

        for crow in cur_rows:
            key = (crow.get(self.KEY, "") or "").strip()
            if not key:
                continue
            target = gb_index.get(key)
            if target is None:
                missing += 1
                continue

            row_changes = 0
            for col in candidate_cols:
                if col not in gb_header:
                    ignored_cols.add(col)
                    continue
                new_val = (crow.get(col, "") or "").strip()
                if new_val == "":
                    continue  # do not overwrite with blank
                old_val = target.get(col, "")
                if new_val != old_val:
                    target[col] = new_val
                    updated_cells += 1
                    row_changes += 1
            if row_changes:
                updated_rows += 1

        # Write output
        out_dir = os.path.join(self.base_dir, self.output_dir)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, self.output_file)
        self._write_tsv(out_path, gb_header, gb_rows)

        # Small summary
        print("\n=== Curation Summary ===")
        print(f"Original rows        : {len(gb_rows)}")
        print(f"Curated rows         : {len(cur_rows)}")
        print(f"Rows updated         : {updated_rows}")
        print(f"Cells updated        : {updated_cells}")
        print(f"Curated keys missing : {missing}")
        if ignored_cols:
            print(f"Ignored curated cols : {sorted(ignored_cols)} (not in matrix header)")
        print(f"Output               : {out_path}")
        print("========================\n")

    # ---------- helpers ----------

    def _read_tsv(self, path: str, label: str) -> Tuple[List[str], List[Dict[str, str]]]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")
        with open(path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            header = reader.fieldnames or []
            rows: List[Dict[str, str]] = []
            for rec in reader:
                # Normalize None to empty strings
                for k in list(rec.keys()):
                    rec[k] = "" if rec[k] is None else str(rec[k])
                rows.append(rec)
        return header, rows

    def _write_tsv(self, path: str, header: List[str], rows: List[Dict[str, str]]) -> None:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t",
                                    quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
            writer.writeheader()
            for rec in rows:
                writer.writerow({c: rec.get(c, "") for c in header})

    def _require_column(self, header: List[str], col: str, label: str) -> None:
        if col not in header:
            raise ValueError(f"Required column '{col}' missing in {label}.")


def _build_argparser() -> ArgumentParser:
    p = ArgumentParser(description="Replace fields in gB_matrix from curated TSV by gi_accession.")
    p.add_argument('-g', '--gb_matrix', help='GenBank matrix file (TSV)',
                   default='tmp/test.txt')#tmp/GenBank-matrix/gB_matrix_raw.tsv')
    p.add_argument('-c', '--curated_file', help='Curated file (TSV)',
                   default='generic/curation.tsv')
    p.add_argument('-b', '--base_dir', help='Base directory for outputs', default='tmp')
    p.add_argument('-t', '--output_dir', help='Output subdirectory', default='Curated')
    p.add_argument('-o', '--output_file', help='Output filename (TSV)', default='gB_matrix.tsv')
    return p


def main() -> None:
    args = _build_argparser().parse_args()
    BlastAlignment(
        args.gb_matrix, args.curated_file, args.base_dir, args.output_dir, args.output_file
    ).process()


if __name__ == "__main__":
    main()

