import csv
from pathlib import Path

import pytest

from FluPivotTable import pivot_data, required_segments


def write_tsv(path: Path, rows, header):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv_dicts(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_pivot_data_complete_and_incomplete(tmp_path: Path):
    input_tsv = tmp_path / "matrix.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    rows = []
    for seg in range(1, 9):
        rows.append([f"A{seg}", "StrainComplete", str(seg), ""])

    rows.extend(
        [
            ["B1", "StrainIncomplete", "1", ""],
            ["B2", "StrainIncomplete", "2.0", ""],
            ["EXCL", "StrainExcluded", "3", "manual exclusion"],
        ]
    )

    write_tsv(
        input_tsv,
        rows,
        header=["primary_accession", "Parsed_strain", "segment_validated", "exclusion"],
    )

    pivot_data(str(input_tsv), str(output_tsv), required_segments)

    out = read_tsv_dicts(output_tsv)
    by_strain = {r["Parsed_strain"]: r for r in out}

    assert by_strain["StrainComplete"]["Complete_status"] == "Complete"
    assert by_strain["StrainComplete"]["1"] == "A1"
    assert by_strain["StrainComplete"]["8"] == "A8"

    assert by_strain["StrainIncomplete"]["Complete_status"] == "Incomplete"
    assert by_strain["StrainIncomplete"]["1"] == "B1"
    assert by_strain["StrainIncomplete"]["2"] == "B2"

    assert "StrainExcluded" not in by_strain


def test_pivot_data_missing_required_columns_raises(tmp_path: Path):
    input_tsv = tmp_path / "bad.tsv"
    output_tsv = tmp_path / "pivot.tsv"

    write_tsv(
        input_tsv,
        [["A1", "StrainX", "1"]],
        header=["primary_accession", "Parsed_strain", "segment_validated"],
    )

    with pytest.raises(ValueError, match="missing one or more required columns"):
        pivot_data(str(input_tsv), str(output_tsv), required_segments)
