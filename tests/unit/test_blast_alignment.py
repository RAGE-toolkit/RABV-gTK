import csv
from pathlib import Path

import pytest

from BlastAlignment import BlastAlignment


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "blast_alignment"


def read_tsv_as_dicts(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def make_processor(tmp_path: Path, gb_matrix_path: Path, is_segmented: str = "N"):
    return BlastAlignment(
        query_fasta=str(DATA_DIR / "query.fa"),
        db_fasta=str(DATA_DIR / "ref.fa"),
        base_dir=str(tmp_path),
        output_dir="blast_out",
        output_file="query_tophits.tsv",
        is_segmented_virus=is_segmented,
        master_acc=str(DATA_DIR / "master_refs.tsv"),
        is_update="N",
        keep_blast_tmp_dir="N",
        gb_matrix=str(gb_matrix_path),
        segment_file=str(DATA_DIR / "master_refs.tsv"),
    )


def test_get_exclusion_list_refs(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    refs = processor.get_exclusion_list_refs()
    assert refs == {"REF_EXCL"}


def test_get_master_list_from_file(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    assert processor.get_master_list() == ["REF_MASTER"]


def test_get_master_list_from_comma_string(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    processor.master_acc = "REF_MASTER,REF_OTHER"
    assert processor.get_master_list() == ["REF_MASTER", "REF_OTHER"]


def test_ref_segments_static():
    segment_map = BlastAlignment.ref_segments(str(DATA_DIR / "query_uniq_tophit_annotated.tsv"))
    assert segment_map == {
        "1": {"REF_MASTER": ["Q1"]},
        "2": {"REF_OTHER": ["Q4"]},
    }


def test_write_filtered_ref_fasta_excludes_exclusion_list(tmp_path: Path):
    processor = make_processor(tmp_path, DATA_DIR / "gb_matrix.tsv")
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    processor.write_filtered_ref_fasta(str(out_dir), {"REF_EXCL"})

    actual = (out_dir / "ref_seq_filtered.fa").read_text(encoding="utf-8")
    expected = (DATA_DIR / "expected_ref_seq_filtered.fa").read_text(encoding="utf-8")
    assert actual == expected


def test_update_gb_matrix_marks_no_hit_and_exclusion_refs(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    out_dir = tmp_path / "blast"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "query_tophits.tsv").write_text(
        (DATA_DIR / "query_tophits.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (out_dir / "query_uniq_tophits.tsv").write_text(
        (DATA_DIR / "query_uniq_tophits.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    processor = make_processor(tmp_path, gb_copy)
    processor.update_gB_matrix(
        query_fasta=str(DATA_DIR / "query.fa"),
        query_tophit_uniq=str(out_dir / "query_uniq_tophits.tsv"),
        gB_matrix_file=str(gb_copy),
    )

    actual = read_tsv_as_dicts(gb_copy)
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gb_matrix_after_update.tsv")
    assert actual == expected


def test_process_non_segmented_orchestration_without_external_tools(tmp_path: Path, monkeypatch):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy, is_segmented="N")
    out_dir = tmp_path / "blast_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(processor, "check_blast_exists", lambda _cmd: True)
    monkeypatch.setattr(processor, "run_makeblastdb", lambda _tmp: None)

    def fake_run_blastn(output_dir, _query_file):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        (Path(output_dir) / "query_tophits.tsv").write_text(
            (DATA_DIR / "query_tophits.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    def fake_process_non_segmented_virus(output_dir, _query_fasta):
        (Path(output_dir) / "query_uniq_tophits.tsv").write_text(
            (DATA_DIR / "query_uniq_tophits.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    monkeypatch.setattr(processor, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(processor, "process_non_segmented_virus", fake_process_non_segmented_virus)

    processor.process()

    actual_matrix = read_tsv_as_dicts(gb_copy)
    expected_matrix = read_tsv_as_dicts(DATA_DIR / "expected_gb_matrix_after_update.tsv")
    assert actual_matrix == expected_matrix

    filtered_fa = out_dir / "ref_seq_filtered.fa"
    assert filtered_fa.exists()
    assert "REF_EXCL" not in filtered_fa.read_text(encoding="utf-8")


def test_update_gb_matrix_raises_on_malformed_unique_hits_row(tmp_path: Path):
    gb_copy = tmp_path / "gb_matrix.tsv"
    gb_copy.write_text((DATA_DIR / "gb_matrix.tsv").read_text(encoding="utf-8"), encoding="utf-8")

    out_dir = tmp_path / "blast"
    out_dir.mkdir(parents=True, exist_ok=True)
    bad_uniq = out_dir / "query_uniq_tophits.tsv"
    bad_uniq.write_text("Q1\tREF_MASTER\t99.9\n", encoding="utf-8")

    processor = make_processor(tmp_path, gb_copy)
    with pytest.raises(ValueError, match="Malformed unique BLAST hit row"):
        processor.update_gB_matrix(
            query_fasta=str(DATA_DIR / "query.fa"),
            query_tophit_uniq=str(bad_uniq),
            gB_matrix_file=str(gb_copy),
        )


def test_process_raises_when_gb_matrix_missing_gi_number(tmp_path: Path):
    gb_bad = tmp_path / "gb_bad.tsv"
    gb_bad.write_text("primary_accession\tsequence\nQ1\tATGC\n", encoding="utf-8")

    processor = make_processor(tmp_path, gb_bad)
    with pytest.raises(ValueError, match="missing required columns: gi_number"):
        processor.process()
