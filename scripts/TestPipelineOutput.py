#!/usr/bin/env python3

import argparse
import csv
import os
import sqlite3
import sys
from typing import Callable, List, Optional, Tuple


def count_columns_tsv(path: str) -> int:
    with open(path, "r", encoding="utf-8") as handle:
        first = handle.readline().rstrip("\n")
    if not first:
        return 0
    return len(first.split("\t"))


def read_tsv_header(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as handle:
        first = handle.readline().rstrip("\n")
    return first.split("\t") if first else []


def iter_tsv_rows(path: str):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            yield row


def scalar(cur: sqlite3.Cursor, sql: str) -> int:
    cur.execute(sql)
    row = cur.fetchone()
    return row[0] if row else 0


def db_checks(db_path: str, segmented: bool) -> Tuple[List[str], List[str]]:
    lines: List[str] = []
    errors: List[str] = []

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    alignment_rows = scalar(cur, "SELECT COUNT(*) FROM sequence_alignment")
    distinct_seq_ids = scalar(cur, "SELECT COUNT(DISTINCT sequence_id) FROM sequence_alignment")
    meta_rows = scalar(cur, "SELECT COUNT(*) FROM meta_data")
    excluded_rows = scalar(cur, "SELECT COUNT(*) FROM excluded_accessions")

    excluded_unprojectable = scalar(
        cur,
        """
        SELECT COUNT(*)
        FROM excluded_accessions
        WHERE LOWER(COALESCE(reason, '')) LIKE '%cannot be projected into merged segment alignment%'
        """,
    )

    excluded_not_enough_matches = scalar(
        cur,
        """
        SELECT COUNT(*)
        FROM excluded_accessions
        WHERE LOWER(COALESCE(reason, '')) LIKE '%unable to align: not enough matches%'
        """,
    )

    excluded_other = excluded_rows - excluded_unprojectable - excluded_not_enough_matches

    missing_from_alignment = scalar(
        cur,
        """
        SELECT COUNT(*) FROM meta_data m
        WHERE m.primary_accession IS NOT NULL
          AND TRIM(m.primary_accession) <> ''
          AND NOT EXISTS (
              SELECT 1 FROM sequence_alignment s
              WHERE s.sequence_id = m.primary_accession
          )
          AND NOT EXISTS (
              SELECT 1 FROM excluded_accessions e
              WHERE e.primary_accession = m.primary_accession
          )
        """,
    )

    missing_alignment_ref_meta = scalar(
        cur,
        """
        SELECT COUNT(*)
        FROM sequence_alignment s
        LEFT JOIN meta_data r ON r.primary_accession = s.alignment_name
        WHERE s.alignment_name IS NOT NULL
          AND TRIM(s.alignment_name) <> ''
          AND r.primary_accession IS NULL
        """,
    )

    lines.append(f"  - sequence_alignment rows: {alignment_rows}")
    lines.append(f"  - unique aligned sequence_id values: {distinct_seq_ids}")
    lines.append(f"  - meta_data rows: {meta_rows}")
    lines.append(f"  - excluded_accessions rows: {excluded_rows}")
    lines.append(f"    - excluded: unprojectable reference projection: {excluded_unprojectable}")
    lines.append(f"    - excluded: nextalign not-enough-matches: {excluded_not_enough_matches}")
    lines.append(f"    - excluded: other reasons: {excluded_other}")
    lines.append(f"  - non-excluded meta_data accessions missing from sequence_alignment: {missing_from_alignment}")
    lines.append(f"  - alignment_name values missing in meta_data: {missing_alignment_ref_meta}")

    segment_mismatch = 0
    if segmented:
        segment_mismatch = scalar(
            cur,
            """
            SELECT COUNT(*)
            FROM sequence_alignment s
            JOIN meta_data q ON q.primary_accession = s.sequence_id
            JOIN meta_data r ON r.primary_accession = s.alignment_name
            WHERE COALESCE(TRIM(q.segment), '') <> ''
              AND COALESCE(TRIM(r.segment), '') <> ''
              AND TRIM(q.segment) <> TRIM(r.segment)
            """,
        )
        lines.append(f"  - segment mismatches (query vs alignment_name): {segment_mismatch}")

    cur.execute(
        """
        SELECT reason, COUNT(*)
        FROM excluded_accessions
        GROUP BY reason
        ORDER BY COUNT(*) DESC
        LIMIT 5
        """
    )
    top_reasons = cur.fetchall()
    if top_reasons:
        lines.append("  - top exclusion reasons (count | reason):")
        for reason, count in top_reasons:
            reason_text = (reason or "").strip().replace("\n", " ")
            if len(reason_text) > 140:
                reason_text = reason_text[:137] + "..."
            lines.append(f"    - {count} | {reason_text}")

    if alignment_rows == 0:
        errors.append("sequence_alignment is empty")
    if missing_from_alignment != 0:
        errors.append(f"{missing_from_alignment} meta_data sequences missing alignment")
    if missing_alignment_ref_meta != 0:
        errors.append(f"{missing_alignment_ref_meta} alignment_name entries not present in meta_data")
    if segmented and segment_mismatch != 0:
        errors.append(f"{segment_mismatch} sequence_alignment rows map across different segments")

    conn.close()
    return lines, errors


def run_segmented(args) -> Tuple[List[str], bool]:
    out: List[str] = ["=== Testing Segmented Virus Pipeline Output ===", ""]

    out.append("Test 1: Checking annotated BLAST file structure...")
    cols = count_columns_tsv(args.annotated_blast)
    if cols == 5:
        out.append("✓ PASS: Annotated BLAST file has 5 columns (query, reference, score, strand, segment)")
    else:
        out.append(f"✗ FAIL: Annotated BLAST file has {cols} columns, expected 5")
        return out, False
    out.append("")

    out.append("Test 2: Checking segment_validated column in matrix...")
    header = read_tsv_header(args.validated_matrix)
    if "segment_validated" in header:
        out.append("✓ PASS: segment_validated column exists")
        seg_idx = header.index("segment_validated")
        total_count = 0
        segment_count = 0
        for row in iter_tsv_rows(args.validated_matrix):
            total_count += 1
            val = row[seg_idx].strip() if seg_idx < len(row) else ""
            if val and val != "not found":
                segment_count += 1
        out.append(f"  - Found {segment_count} records with valid segments out of {total_count} total")
        if segment_count > 0:
            out.append("✓ PASS: At least some records have segment assignments")
        else:
            out.append("⚠ WARNING: No records have valid segment assignments")
    else:
        out.append("✗ FAIL: segment_validated column not found in matrix")
        return out, False
    out.append("")

    if args.pivoted_matrix and os.path.exists(args.pivoted_matrix):
        out.append("Test 3: Checking pivoted segments matrix...")
        p_header = read_tsv_header(args.pivoted_matrix)
        if "Complete_status" in p_header:
            out.append("✓ PASS: Pivoted matrix has Complete_status column")
        else:
            out.append("✗ FAIL: Pivoted matrix missing Complete_status column")
            return out, False

        segment_cols = sum(1 for col in p_header if col in {str(x) for x in range(1, 9)})
        out.append(f"  - Found {segment_cols} segment columns")

        status_idx = p_header.index("Complete_status")
        complete = 0
        incomplete = 0
        for row in iter_tsv_rows(args.pivoted_matrix):
            status = row[status_idx].strip() if status_idx < len(row) else ""
            if status == "Complete":
                complete += 1
            elif status == "Incomplete":
                incomplete += 1

        out.append(f"  - Complete genomes: {complete}")
        out.append(f"  - Incomplete genomes: {incomplete}")
        if complete + incomplete > 0:
            out.append("✓ PASS: Pivoted matrix contains strain data")
        else:
            out.append("⚠ WARNING: Pivoted matrix is empty")
    else:
        out.append("Test 3: SKIPPED - Pivoted matrix not expected for this run")
    out.append("")

    out.append("Test 4: Checking sequence_alignment table in SQLite DB...")
    db_lines, db_errors = db_checks(args.sqlite_db, segmented=True)
    out.extend(db_lines)
    if db_errors:
        out.append("✗ FAIL: " + "; ".join(db_errors))
        return out, False
    out.append("✓ PASS: sequence_alignment covers all non-excluded metadata and segment-consistent mapping")
    out.append("")

    out.append("=== All segmented virus tests completed ===")
    return out, True


def run_non_segmented(args) -> Tuple[List[str], bool]:
    out: List[str] = ["=== Testing Non-Segmented Virus Pipeline Output ===", ""]

    out.append("Test 1: Checking BLAST file structure...")
    cols = count_columns_tsv(args.blast_hits)
    if cols == 4:
        out.append("✓ PASS: BLAST file has 4 columns (query, reference, score, strand)")
    else:
        out.append(f"✗ FAIL: BLAST file has {cols} columns, expected 4")
        return out, False
    out.append("")

    out.append("Test 2: Verifying no segment column for non-segmented virus...")
    gb_header = read_tsv_header(args.gb_matrix)
    if "segment_validated" in gb_header:
        out.append("⚠ WARNING: segment_validated column found (unexpected for non-segmented virus)")
    else:
        out.append("✓ PASS: No segment_validated column (correct for non-segmented virus)")
    out.append("")

    out.append("Test 3: Checking sequence counts...")
    blast_count = sum(1 for _ in open(args.blast_hits, "r", encoding="utf-8"))
    matrix_count = sum(1 for _ in iter_tsv_rows(args.gb_matrix))
    out.append(f"  - BLAST hits: {blast_count}")
    out.append(f"  - Matrix records: {matrix_count}")
    if blast_count > 0 and matrix_count > 0:
        out.append("✓ PASS: Pipeline processed sequences")
    else:
        out.append("⚠ WARNING: No sequences processed")
    out.append("")

    out.append("Test 4: Checking sequence_alignment table in SQLite DB...")
    db_lines, db_errors = db_checks(args.sqlite_db, segmented=False)
    out.extend(db_lines)
    if db_errors:
        out.append("✗ FAIL: " + "; ".join(db_errors))
        return out, False
    out.append("✓ PASS: sequence_alignment covers all non-excluded meta_data sequences")
    out.append("")

    out.append("=== All non-segmented virus tests completed ===")
    return out, True


def write_output(lines: List[str], output_path: str) -> None:
    with open(output_path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pipeline output validator for test profiles")
    parser.add_argument("--mode", choices=["segmented", "non_segmented"], required=True)
    parser.add_argument("--sqlite_db", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument("--annotated_blast")
    parser.add_argument("--validated_matrix")
    parser.add_argument("--pivoted_matrix")

    parser.add_argument("--blast_hits")
    parser.add_argument("--gb_matrix")
    return parser


def require_args(args, names: List[str]) -> Optional[str]:
    missing = [name for name in names if not getattr(args, name)]
    if missing:
        return "Missing required arguments for mode {}: {}".format(args.mode, ", ".join(missing))
    return None


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "segmented":
        err = require_args(args, ["annotated_blast", "validated_matrix"])
        if err:
            print(err, file=sys.stderr)
            return 2
        lines, ok = run_segmented(args)
    else:
        err = require_args(args, ["blast_hits", "gb_matrix"])
        if err:
            print(err, file=sys.stderr)
            return 2
        lines, ok = run_non_segmented(args)

    write_output(lines, args.output)
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
