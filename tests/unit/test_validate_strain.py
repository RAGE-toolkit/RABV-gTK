import csv
from pathlib import Path

from ValidateStrain import (
    clean_strain_name,
    curate_serotype,
    extract_serotype,
    parse_strain_from_definition,
    process_matrix,
)


def write_tsv(path: Path, rows, header):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv_dicts(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_parse_and_clean_strain():
    definition = "Influenza A virus (A/chicken/Chile/11(H5N1)) segment"
    parsed = parse_strain_from_definition(definition)
    assert parsed == "A/chicken/Chile/11"

    assert clean_strain_name("A/chicken Chile///11___") == "A/chicken_Chile///11"


def test_extract_serotype_and_curation_edge_cases():
    assert extract_serotype({"serotype": " H3N2 ", "definition": "x"}) == "H3N2"
    assert extract_serotype({"serotype": "", "definition": "virus (H10N7)"}) == "H10N7"
    assert extract_serotype({"definition": "virus without subtype"}) == ""

    replacements = {"HXNY": "unknown"}
    assert curate_serotype("A/H1N1", replacements) == "H1N1"
    assert curate_serotype("h5n1", replacements) == "H5N1"
    assert curate_serotype("H7", replacements) == "H7"
    assert curate_serotype("N2", replacements) == "N2"
    assert curate_serotype("H?N3", replacements) == "N3"
    assert curate_serotype("H9N?", replacements) == "H9"
    assert curate_serotype("HXNY", replacements) == "unknown"
    assert curate_serotype("???", replacements) == "unknown"


def test_process_matrix_end_to_end(tmp_path: Path):
    input_tsv = tmp_path / "gB_matrix.tsv"
    output_tsv = tmp_path / "gB_matrix_validated.tsv"

    write_tsv(
        input_tsv,
        [
            ["ACC1", "A/duck/CA/1", "", "Influenza A virus (A/duck/CA/1(H3N2))"],
            ["ACC2", "", "", "Influenza A virus (A/chicken/PE/99(H5N1))"],
            ["ACC3", "", "A/H1N1", "Influenza A virus (A/human/US/20(H1N1))"],
            ["ACC4", "", "", "no subtype info"],
        ],
        header=["primary_accession", "strain", "serotype", "definition"],
    )

    process_matrix(
        str(input_tsv),
        str(output_tsv),
        overwrite=False,
        serotype_replacements={"UNKNOWNVALUE": "unknown"},
    )

    rows = read_tsv_dicts(output_tsv)

    assert rows[0]["Parsed_strain"] == "A/duck/CA/1"
    assert rows[0]["serotype_validated"] == "H3N2"

    assert rows[1]["Parsed_strain"] == "A/chicken/PE/99"
    assert rows[1]["serotype_validated"] == "H5N1"

    assert rows[2]["Parsed_strain"] == "A/human/US/20"
    assert rows[2]["serotype_validated"] == "H1N1"

    assert rows[3]["Parsed_strain"] == ""
    assert rows[3]["serotype_validated"] == "unknown"


def test_process_matrix_overwrite_mode(tmp_path: Path):
    input_tsv = tmp_path / "gB_matrix.tsv"
    write_tsv(
        input_tsv,
        [["ACCX", "", "", "Influenza A virus (A/goose/CN/7(H9N2))"]],
        header=["primary_accession", "strain", "serotype", "definition"],
    )

    process_matrix(str(input_tsv), str(input_tsv), overwrite=True, serotype_replacements={})
    rows = read_tsv_dicts(input_tsv)

    assert rows[0]["Parsed_strain"] == "A/goose/CN/7"
    assert rows[0]["serotype_validated"] == "H9N2"
