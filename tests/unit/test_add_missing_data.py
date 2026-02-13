import csv
import subprocess
import sys
from pathlib import Path

import pytest

from AddMissingData import AddMissingData


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "add_missing_data"
SCRIPT_PATH = REPO_ROOT / "scripts" / "AddMissingData.py"


def read_tsv_as_dicts(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_get_header_and_header_lookup():
    processor = AddMissingData(tmp_dir=".", gb_matrix=str(DATA_DIR / "input_gb_matrix.tsv"))
    header = processor.get_header(str(DATA_DIR / "fillup.tsv"))
    assert header == ["primary_accession", "country", "host", "collection_date"]

    keys = processor.header_lookup(header)
    assert keys == ["primary_accession", "country", "host", "collection_date"]


def test_add_missing_values_matches_expected(tmp_path: Path):
    out_dir = tmp_path / "out"
    processor = AddMissingData(
        tmp_dir=str(out_dir),
        gb_matrix=str(DATA_DIR / "input_gb_matrix.tsv"),
        fillup_file=str(DATA_DIR / "fillup.tsv"),
    )
    processor.add_missing_values()

    actual = read_tsv_as_dicts(out_dir / "gB_matrix_updated.csv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gB_matrix_updated.tsv")
    assert actual == expected


def test_bulk_replace_matches_expected(tmp_path: Path):
    out_dir = tmp_path / "out"
    processor = AddMissingData(
        tmp_dir=str(out_dir),
        gb_matrix=str(DATA_DIR / "input_gb_matrix.tsv"),
        bulk_file=str(DATA_DIR / "bulk.tsv"),
    )
    processor.bulk_replace()

    actual = read_tsv_as_dicts(out_dir / "gB_matrix_replaced.tsv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gB_matrix_replaced.tsv")
    assert actual == expected


def test_process_requires_either_fillup_or_bulk(tmp_path: Path):
    processor = AddMissingData(
        tmp_dir=str(tmp_path / "out"),
        gb_matrix=str(DATA_DIR / "input_gb_matrix.tsv"),
    )
    with pytest.raises(ValueError, match="Either fillup_file or bulk_file must be provided"):
        processor.process()


def test_cli_fillup_mode(tmp_path: Path):
    out_dir = tmp_path / "cli_fillup"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-d",
            str(out_dir),
            "-g",
            str(DATA_DIR / "input_gb_matrix.tsv"),
            "-f",
            str(DATA_DIR / "fillup.tsv"),
        ],
        check=True,
    )

    actual = read_tsv_as_dicts(out_dir / "gB_matrix_updated.csv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gB_matrix_updated.tsv")
    assert actual == expected


def test_cli_bulk_mode(tmp_path: Path):
    out_dir = tmp_path / "cli_bulk"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-d",
            str(out_dir),
            "-g",
            str(DATA_DIR / "input_gb_matrix.tsv"),
            "-b",
            str(DATA_DIR / "bulk.tsv"),
        ],
        check=True,
    )

    actual = read_tsv_as_dicts(out_dir / "gB_matrix_replaced.tsv")
    expected = read_tsv_as_dicts(DATA_DIR / "expected_gB_matrix_replaced.tsv")
    assert actual == expected
