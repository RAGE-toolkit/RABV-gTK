#!/usr/bin/env python3
import os
import csv
import re
from argparse import ArgumentParser
from typing import List, Dict, Tuple


class BlastAlignment:
    KEY = "gi_number"
    COMMENT_COL = "comment"   # the free-text comment column name (both files)
    _BLANKISH = {"", "-", "na", "n/a", "NA", "Na", "N/A"}

    def __init__(self, gb_matrix: str, curated_file: str,
                 base_dir: str = "tmp", output_dir: str = "Curated",
                 output_file: str = "gB_matrix.tsv") -> None:
        self.gb_matrix = gb_matrix
        self.curated_file = curated_file
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.output_file = output_file

    # ---------------------- Public API ----------------------

    def process(self) -> None:
        gb_header, gb_rows = self._read_tsv(self.gb_matrix, "GenBank matrix")
        cur_header, cur_rows = self._read_tsv(self.curated_file, "curated file")

        # Required columns
        self._require_column(gb_header, self.KEY, "GenBank matrix")
        self._require_column(cur_header, self.KEY, "curated file")

        # Optional columns
        comment_col = self.COMMENT_COL if self.COMMENT_COL in gb_header else None
        curator_col = "curator" if "curator" in gb_header else None

        # Index GenBank matrix by accession
        gb_index: Dict[str, Dict[str, str]] = {}
        for row in gb_rows:
            k = (row.get(self.KEY, "") or "").strip()
            if k and k not in gb_index:
                gb_index[k] = row

        # Columns to update (excluding comment and curator)
        candidate_cols = [c for c in cur_header if c not in (self.KEY, self.COMMENT_COL, "curator")]

        updated_rows = 0
        updated_cells = 0
        missing = 0
        ignored_cols = set()

        for crow in cur_rows:
            key = (crow.get(self.KEY, "") or "").strip()
            if not key:
                continue
            target = gb_index.get(key)
            if target is None:
                missing += 1
                continue

            row_changes = 0
            change_notes: List[str] = []

            # --- 1) Update metadata fields ---
            for col in candidate_cols:
                if col not in gb_header:
                    ignored_cols.add(col)
                    continue

                new_val_raw = crow.get(col, "")
                new_val = (new_val_raw if new_val_raw is not None else "").strip()
                if self._is_blankish(new_val):
                    continue

                old_val_raw = target.get(col, "")
                old_val = (old_val_raw if old_val_raw is not None else "").strip()

                if new_val != old_val:
                    target[col] = new_val
                    updated_cells += 1
                    row_changes += 1
                    change_notes.append(f"{col}: '{old_val}' -> '{new_val}'")

            # --- 2) Only append comment/curator if real metadata changed ---
            if row_changes > 0:
                updated_rows += 1

                # Handle comment
                if comment_col:
                    existing_comment = (target.get(comment_col, "") or "").strip()
                    curated_comment = (crow.get(comment_col, "") or "").strip()
                    existing_comment = self._strip_unwanted_fragments(existing_comment)

                    if not self._is_blankish(curated_comment):
                        existing_comment = self._append_unique(existing_comment, curated_comment)

                    if change_notes:
                        change_blob = "; ".join(change_notes)
                        existing_comment = self._append_unique(existing_comment, change_blob)

                    target[comment_col] = self._clean_separators(existing_comment)

                # Handle curator
                if curator_col:
                    existing_curator = (target.get(curator_col, "") or "").strip()
                    curated_curator = (crow.get(curator_col, "") or "").strip()
                    if not self._is_blankish(curated_curator):
                        target[curator_col] = self._append_unique(existing_curator, curated_curator)

        # --- Write output ---
        out_dir = os.path.join(self.base_dir, self.output_dir)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, self.output_file)
        self._write_tsv(out_path, gb_header, gb_rows)

        # --- Summary ---
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

    # ---------------------- Utilities ----------------------

    @staticmethod
    def _is_blankish(val: str) -> bool:
        return (val is None) or (val.strip() == "") or (val.strip() in BlastAlignment._BLANKISH) or (val.strip().lower() in {"na", "n/a"})

    @staticmethod
    def _split_comment_items(s: str) -> List[str]:
        parts = [p.strip() for p in s.split(";") if p.strip()]
        return parts

    def _append_unique(self, existing: str, new_item: str) -> str:
        if self._is_blankish(new_item):
            return existing
        existing_items = self._split_comment_items(existing)
        existing_norm = {re.sub(r"\s+", " ", it.lower()) for it in existing_items}
        new_norm = re.sub(r"\s+", " ", new_item.strip().lower())
        if new_norm not in existing_norm:
            existing_items.append(new_item.strip())
        return "; ".join(existing_items)

    @staticmethod
    def _clean_separators(s: str) -> str:
        s = re.sub(r"(^\s*NA\s*;?\s*|\s*;?\s*NA\s*$)", "", s, flags=re.IGNORECASE)
        s = re.sub(r"\s*;\s*;\s*", "; ", s)
        return s.strip(" ;")

    @staticmethod
    def _strip_unwanted_fragments(s: str) -> str:
        if not s:
            return s
        items = [p for p in BlastAlignment._split_comment_items(s) if not p.lower().strip().startswith("curator:")]
        return "; ".join(items)

    # ---------- I/O helpers ----------

    def _read_tsv(self, path: str, label: str) -> Tuple[List[str], List[Dict[str, str]]]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")
        with open(path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            header = reader.fieldnames or []
            rows: List[Dict[str, str]] = []
            for rec in reader:
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

    @staticmethod
    def _require_column(header: List[str], col: str, label: str) -> None:
        if col not in header:
            raise ValueError(f"Required column '{col}' missing in {label}.")


# ---------------------- CLI Wrapper ----------------------

def _build_argparser() -> ArgumentParser:
    p = ArgumentParser(description=(
        "Replace fields in gB_matrix from curated TSV by gi_accession. "
        "Curator/comment fields are appended only if other metadata fields differ."
    ))
    p.add_argument('-g', '--gb_matrix', help='GenBank matrix file (TSV)',
                   default='tmp/GenBank-matrix/gB_matrix_raw.tsv')
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

