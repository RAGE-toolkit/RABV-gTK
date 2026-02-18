from pathlib import Path

import pytest

from CalcGenomeCords import CalculateGenomeCoordinates


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "calc_genome_cords"


def test_extract_alignment_coordinates_with_explicit_master():
    calc = CalculateGenomeCoordinates(
        paded_alignment=str(DATA_DIR / "alignment.fa"),
        master_accession="MASTER1",
    )
    coords = calc.extract_alignment_coordinates()

    assert coords["MASTER1"] == ["MASTER1", 1, 9]
    assert coords["Q_A"] == ["MASTER1", 1, 9]
    assert coords["Q_B"] == ["MASTER1", 3, 7]


def test_extract_alignment_coordinates_default_master_is_first_record():
    calc = CalculateGenomeCoordinates(paded_alignment=str(DATA_DIR / "alignment.fa"))
    coords = calc.extract_alignment_coordinates()

    assert calc.master_accession == "MASTER1"
    assert coords["Q_A"] == ["MASTER1", 1, 9]


def test_extract_alignment_coordinates_raises_when_master_missing():
    calc = CalculateGenomeCoordinates(
        paded_alignment=str(DATA_DIR / "alignment.fa"),
        master_accession="DOES_NOT_EXIST",
    )
    with pytest.raises(ValueError, match="Master sequence with ID"):
        calc.extract_alignment_coordinates()


def test_extract_alignment_coordinates_raises_when_file_missing(tmp_path: Path):
    calc = CalculateGenomeCoordinates(
        paded_alignment=str(tmp_path / "missing.fa"),
        master_accession="MASTER1",
    )
    with pytest.raises(FileNotFoundError, match="Padded alignment file not found"):
        calc.extract_alignment_coordinates()


def test_extract_alignment_coordinates_raises_when_alignment_empty(tmp_path: Path):
    empty = tmp_path / "empty.fa"
    empty.write_text("", encoding="utf-8")
    calc = CalculateGenomeCoordinates(paded_alignment=str(empty), master_accession="MASTER1")
    with pytest.raises(ValueError, match="No FASTA records found"):
        calc.extract_alignment_coordinates()
