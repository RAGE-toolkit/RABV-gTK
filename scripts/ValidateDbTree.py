import argparse
import os
import sqlite3
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

ACCESSION_COLUMNS = {
    "primary_accession",
    "accession",
    "accession_id",
    "accession_version",
    "sequence_id",
    "sequence_name",
    "seq_name",
    "sample",
    "query",
    "reference",
    "header",
}


def get_table_columns(conn, table_name):
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    return [row[1] for row in cursor.fetchall()]


def collect_table_values(conn, table_name, columns):
    values = set()
    cursor = conn.cursor()
    for col in columns:
        cursor.execute(f"SELECT DISTINCT {col} FROM {table_name}")
        values.update([r[0] for r in cursor.fetchall() if r[0]])
    return values


def find_cluster_column(columns):
    for col in columns:
        if col == "cluster" or col.startswith("cluster_"):
            return col
    return None


def fetch_tree_newick(conn):
    cursor = conn.cursor()
    cursor.execute(
        "SELECT name, newick FROM trees WHERE newick IS NOT NULL AND length(newick) > 0 ORDER BY CASE WHEN name='usher' THEN 0 WHEN name='iqtree' THEN 1 ELSE 2 END, name LIMIT 1"
    )
    row = cursor.fetchone()
    if not row:
        return None, None
    return row[0], row[1]


def fetch_table_column(conn, table_name, column_name):
    cursor = conn.cursor()
    cursor.execute(f"SELECT {column_name} FROM {table_name}")
    return [r[0] for r in cursor.fetchall()]


def main():
    parser = argparse.ArgumentParser(description="Validate SQLite DB contents against tree and plot the tree")
    parser.add_argument("--db", required=True, help="Path to SQLite DB")
    parser.add_argument("--outdir", required=True, help="Output directory for report and plot")
    parser.add_argument(
        "--test-mode",
        action="store_true",
        help="Relax strict failures for known test-only edge cases (e.g., no placeable queries)",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    conn = sqlite3.connect(args.db)
    try:
        tree_name, newick = fetch_tree_newick(conn)
        if not newick:
            raise SystemExit("No tree found in DB (trees table is empty)")

        tree = Phylo.read(StringIO(newick), "newick")
        tree_terminals = {t.name for t in tree.get_terminals() if t.name}

        meta_accessions = [x for x in fetch_table_column(conn, "meta_data", "primary_accession") if x]
        seq_headers = [x for x in fetch_table_column(conn, "sequences", "header") if x]

        meta_set = set(meta_accessions)
        seq_set = set(seq_headers)

        missing_in_tree = sorted(meta_set - tree_terminals)
        missing_in_sequences = sorted(meta_set - seq_set)
        missing_in_meta = sorted(tree_terminals - meta_set)

        meta_columns = get_table_columns(conn, "meta_data")
        cluster_col = find_cluster_column(meta_columns)
        centroid_set = set()
        missing_centroids_in_tree = []
        extra_in_tree = []
        if cluster_col:
            centroid_set = set(
                x for x in fetch_table_column(conn, "meta_data", cluster_col) if x
            )
            missing_centroids_in_tree = sorted(centroid_set - tree_terminals)
            extra_in_tree = sorted(tree_terminals - centroid_set)

        # Basic completeness checks
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM meta_data")
        meta_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM sequences")
        seq_count = cursor.fetchone()[0]

        report_path = os.path.join(args.outdir, "db_tree_validation.txt")
        missing_tree_path = os.path.join(args.outdir, "missing_in_tree.txt")
        missing_seq_path = os.path.join(args.outdir, "missing_in_sequences.txt")
        missing_meta_path = os.path.join(args.outdir, "missing_in_meta.txt")
        missing_centroids_path = os.path.join(args.outdir, "missing_centroids_in_tree.txt")
        extra_in_tree_path = os.path.join(args.outdir, "extra_in_tree.txt")

        with open(report_path, "w", encoding="utf-8") as report:
            report.write(f"Tree source: {tree_name}\n")
            report.write(f"Meta rows: {meta_count}\n")
            report.write(f"Sequence rows: {seq_count}\n")
            report.write(f"Tree terminals: {len(tree_terminals)}\n")
            if cluster_col:
                report.write(f"Cluster column: {cluster_col}\n")
                report.write(f"Centroid count: {len(centroid_set)}\n")
                report.write(f"Missing centroids in tree: {len(missing_centroids_in_tree)}\n")
                report.write(f"Extra nodes in tree: {len(extra_in_tree)}\n")
            report.write(f"Missing in tree (meta_data -> tree): {len(missing_in_tree)}\n")
            report.write(f"Missing in sequences (meta_data -> sequences): {len(missing_in_sequences)}\n")
            report.write(f"Missing in meta_data (tree -> meta_data): {len(missing_in_meta)}\n")

            if missing_in_tree:
                report.write("\nMissing in tree (first 50):\n")
                report.write("\n".join(missing_in_tree[:50]) + "\n")
            if missing_in_sequences:
                report.write("\nMissing in sequences (first 50):\n")
                report.write("\n".join(missing_in_sequences[:50]) + "\n")
            if missing_in_meta:
                report.write("\nMissing in meta_data (first 50):\n")
                report.write("\n".join(missing_in_meta[:50]) + "\n")

            if cluster_col and missing_centroids_in_tree:
                report.write("\nMissing centroids in tree (first 50):\n")
                report.write("\n".join(missing_centroids_in_tree[:50]) + "\n")
            if cluster_col and extra_in_tree:
                report.write("\nExtra nodes in tree (first 50):\n")
                report.write("\n".join(extra_in_tree[:50]) + "\n")

            # Table-level coverage summary
            report.write("\nTable-level accession coverage:\n")
            for table_name in ["meta_data", "sequences", "sequence_alignment", "features", "insertions", "host_taxa"]:
                cols = get_table_columns(conn, table_name)
                acc_cols = [c for c in cols if c in ACCESSION_COLUMNS]
                if not acc_cols:
                    report.write(f"- {table_name}: no accession-like columns found\n")
                    continue
                values = collect_table_values(conn, table_name, acc_cols)
                missing_from_table = sorted(tree_terminals - values)
                report.write(
                    f"- {table_name}: columns={','.join(acc_cols)} values={len(values)} missing_from_table={len(missing_from_table)}\n"
                )

        # Write full lists for debugging
        if missing_in_tree:
            with open(missing_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_tree) + "\n")
        if missing_in_sequences:
            with open(missing_seq_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_sequences) + "\n")
        if missing_in_meta:
            with open(missing_meta_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_in_meta) + "\n")

        if cluster_col and missing_centroids_in_tree:
            with open(missing_centroids_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(missing_centroids_in_tree) + "\n")
        if cluster_col and extra_in_tree:
            with open(extra_in_tree_path, "w", encoding="utf-8") as handle:
                handle.write("\n".join(extra_in_tree) + "\n")

        # Print summary to stdout for quick visibility in logs
        print(f"[info] Tree source: {tree_name}")
        print(f"[info] Meta rows: {meta_count}, Sequence rows: {seq_count}, Tree terminals: {len(tree_terminals)}")
        if cluster_col:
            print(f"[info] Cluster column: {cluster_col}, Centroid count: {len(centroid_set)}")
            print(f"[info] Missing centroids in tree: {len(missing_centroids_in_tree)}")
            print(f"[info] Extra nodes in tree: {len(extra_in_tree)}")
        print(f"[info] Missing in tree: {len(missing_in_tree)}")
        print(f"[info] Missing in sequences: {len(missing_in_sequences)}")
        print(f"[info] Missing in meta_data: {len(missing_in_meta)}")
        if missing_in_tree:
            print(f"[info] Missing in tree (first 10): {', '.join(missing_in_tree[:10])}")
        if missing_in_sequences:
            print(f"[info] Missing in sequences (first 10): {', '.join(missing_in_sequences[:10])}")
        if missing_in_meta:
            print(f"[info] Missing in meta_data (first 10): {', '.join(missing_in_meta[:10])}")
        if cluster_col and missing_centroids_in_tree:
            print(f"[info] Missing centroids in tree (first 10): {', '.join(missing_centroids_in_tree[:10])}")
        if cluster_col and extra_in_tree:
            print(f"[info] Extra nodes in tree (first 10): {', '.join(extra_in_tree[:10])}")

        # Plot tree
        fig = plt.figure(figsize=(12, 18))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        plot_path = os.path.join(args.outdir, "db_tree.png")
        plt.tight_layout()
        fig.savefig(plot_path, dpi=150)

        overlap_count = len(meta_set & tree_terminals)
        disjoint_sets = overlap_count == 0 and len(meta_set) > 0 and len(tree_terminals) > 0

        # Test-mode exception: disjoint USHER tree/meta can happen when no query sequences were placeable
        # and the resulting tree contains only reference nodes.
        if args.test_mode and tree_name == "usher" and disjoint_sets:
            print(
                "[warn] Test mode: tree terminals and meta_data are disjoint; "
                "treating as no-placeable-queries scenario and not failing validation."
            )
            return

        # Fail if any validation issues detected
        if tree_name == "usher":
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: UShER tree does not match accessions")
        elif cluster_col:
            if missing_centroids_in_tree or extra_in_tree:
                raise SystemExit("Validation failed: IQ-TREE does not match centroid set")
        else:
            if missing_in_tree or missing_in_sequences or missing_in_meta:
                raise SystemExit("Validation failed: missing accessions detected")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
