import csv
from pathlib import Path

from FilterAndExtractSequences import FilterAndExtractSequences


def _write_tsv(path: Path, rows, header):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def _read_tsv(path: Path):
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_filter_columns_writes_query_ref_and_exclusions(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [
            ["VRL", "REF1", "4", "0", "0", ""],
            ["VRL", "Q1", "4", "0", "0", ""],
            ["VRL", "Q2", "3", "0", "0", ""],
        ],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )

    seqs.write_text(
        ">REF1\nATGC\n>Q1\nAATT\n>Q2\nAAA\n",
        encoding="utf-8",
    )
    refs.write_text("REF1\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=4,
        real_length=4,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    query_fa = (tmp_path / "Sequences" / "query_seq.fa").read_text(encoding="utf-8")
    ref_fa = (tmp_path / "Sequences" / "ref_seq.fa").read_text(encoding="utf-8")

    assert ">Q1" in query_fa
    assert ">Q2" not in query_fa
    assert ">REF1" in ref_fa

    rows = _read_tsv(matrix)
    by_acc = {r["gi_number"]: r for r in rows}
    assert by_acc["Q2"]["exclusion_status"] == "1"


def test_segmented_mode_writes_exclusion_refs_file(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [["VRL", "REF1", "4", "0", "0", ""]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">REF1\nATGC\n", encoding="utf-8")
    refs.write_text("REF1\tmaster\t1\nEXCL1\texclusion_list\t1\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="Y",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    exclusion_refs = (tmp_path / "Sequences" / "exclusion_refs.txt").read_text(encoding="utf-8").splitlines()
    assert exclusion_refs == ["EXCL1"]


def test_existing_exclusion_status_is_preserved(tmp_path: Path):
    matrix = tmp_path / "gB_matrix_raw.tsv"
    seqs = tmp_path / "sequences.fa"
    refs = tmp_path / "refs.tsv"

    _write_tsv(
        matrix,
        [["VRL", "QX", "4", "0", "1", "already excluded"]],
        ["division", "gi_number", "length", "n", "exclusion_status", "exclusion_criteria"],
    )
    seqs.write_text(">QX\nATGC\n", encoding="utf-8")
    refs.write_text("REF1\treference\n", encoding="utf-8")

    processor = FilterAndExtractSequences(
        genbank_matrix=str(matrix),
        sequence_file=str(seqs),
        genbank_matrix_filtered=str(tmp_path),
        ref_file=str(refs),
        base_dir=str(tmp_path),
        output_dir="Sequences",
        total_length=1,
        real_length=1,
        prop_ambigious=None,
        segmented_virus="N",
        gb_division=None,
        valid_divisions=["VRL", "ENV"],
        seq_type=None,
    )
    processor.process()

    rows = _read_tsv(matrix)
    assert rows[0]["exclusion_status"] == "1"
    assert rows[0]["exclusion_criteria"] == "already excluded"
