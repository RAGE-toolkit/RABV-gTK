import csv
from pathlib import Path

import pandas as pd
import pytest

from GenBankParser import GenBankParser


def _gbseq_xml(primary_accession: str, sequence: str, country: str = "UK:London", collection_date: str = "12-Jan-2020") -> str:
    return f"""
<GBSeq>
  <GBSeq_locus>{primary_accession}</GBSeq_locus>
  <GBSeq_length>{len(sequence)}</GBSeq_length>
  <GBSeq_moltype>genomic RNA</GBSeq_moltype>
  <GBSeq_topology>linear</GBSeq_topology>
  <GBSeq_division>VRL</GBSeq_division>
  <GBSeq_update-date>01-JAN-2024</GBSeq_update-date>
  <GBSeq_create-date>01-JAN-2020</GBSeq_create-date>
  <GBSeq_definition>Example definition {primary_accession}</GBSeq_definition>
  <GBSeq_primary-accession>{primary_accession}</GBSeq_primary-accession>
  <GBSeq_accession-version>{primary_accession}.1</GBSeq_accession-version>
  <GBSeq_source>Virus</GBSeq_source>
  <GBSeq_organism>Virus organism</GBSeq_organism>
  <GBSeq_taxonomy>Viruses; Example</GBSeq_taxonomy>
  <GBSeq_feature-table>
    <GBFeature>
      <GBFeature_key>source</GBFeature_key>
      <GBFeature_quals>
        <GBQualifier><GBQualifier_name>country</GBQualifier_name><GBQualifier_value>{country}</GBQualifier_value></GBQualifier>
        <GBQualifier><GBQualifier_name>host</GBQualifier_name><GBQualifier_value>bat</GBQualifier_value></GBQualifier>
        <GBQualifier><GBQualifier_name>collection_date</GBQualifier_name><GBQualifier_value>{collection_date}</GBQualifier_value></GBQualifier>
      </GBFeature_quals>
    </GBFeature>
    <GBFeature>
      <GBFeature_key>gene</GBFeature_key>
      <GBFeature_location>1..3</GBFeature_location>
      <GBFeature_quals>
        <GBQualifier><GBQualifier_name>gene</GBQualifier_name><GBQualifier_value>N</GBQualifier_value></GBQualifier>
      </GBFeature_quals>
    </GBFeature>
  </GBSeq_feature-table>
  <GBSeq_references>
    <GBReference>
      <GBReference_reference>1</GBReference_reference>
      <GBReference_position>1..{len(sequence)}</GBReference_position>
      <GBReference_authors><GBAuthor>A. Author</GBAuthor></GBReference_authors>
      <GBReference_title>Test title</GBReference_title>
      <GBReference_journal>Test journal</GBReference_journal>
      <GBReference_pubmed>12345</GBReference_pubmed>
    </GBReference>
  </GBSeq_references>
  <GBSeq_sequence>{sequence.lower()}</GBSeq_sequence>
</GBSeq>
"""


def _write_xml(path: Path, gbseq_blocks):
    xml = "<GBSet>\n" + "\n".join(gbseq_blocks) + "\n</GBSet>\n"
    path.write_text(xml, encoding="utf-8")


def test_load_ref_list_supports_two_and_three_column_lines(tmp_path: Path):
    ref_file = tmp_path / "ref_list.tsv"
    ref_file.write_text("REF1\tmaster\t1\nREF2\treference\nBADLINE\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
    )

    refs = parser.load_ref_list(str(ref_file))
    assert refs == {"REF1": "master", "REF2": "reference"}


def test_xml_to_tsv_parses_and_applies_exclusion(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(
        xml_file,
        [
            _gbseq_xml("REF1", "ATGCNN"),
            _gbseq_xml("Q1", "AATTCG", country="France:Paris", collection_date="2021"),
        ],
    )

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
    )

    rows = parser.xml_to_tsv(str(xml_file), {"REF1": "master"}, ["Q1"])
    by_acc = {r["primary_accession"]: r for r in rows}

    assert by_acc["REF1"]["accession_type"] == "master"
    assert by_acc["Q1"]["accession_type"] == "excluded"
    assert by_acc["Q1"]["exclusion_status"] == "1"
    assert by_acc["Q1"]["country"] == "France"
    assert by_acc["Q1"]["geo_loc"] == "Paris"
    assert by_acc["REF1"]["n"] == 2


def test_select_xml_files_raises_on_missing_critical_refs(tmp_path: Path):
    xml_file = tmp_path / "batch-1.xml"
    _write_xml(xml_file, [_gbseq_xml("REF1", "ATGC")])

    parser = GenBankParser(
        input_dir=str(tmp_path),
        base_dir=str(tmp_path),
        output_dir="out",
        ref_list=None,
        exclusion_list=None,
        is_segmented_virus="N",
        test_run=True,
        test_xml_limit=2,
        require_refs=True,
    )

    with pytest.raises(ValueError, match="Missing CRITICAL reference accessions"):
        parser.select_xml_files({"REF1": "master", "REF_MISSING": "reference"}, ["batch-1.xml"])


def test_process_writes_matrix_and_fasta(tmp_path: Path):
    in_dir = tmp_path / "xml"
    in_dir.mkdir(parents=True, exist_ok=True)

    _write_xml(in_dir / "batch-1.xml", [_gbseq_xml("REF1", "ATGC")])
    _write_xml(in_dir / "batch-2.xml", [_gbseq_xml("Q1", "AATT")])

    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    parser = GenBankParser(
        input_dir=str(in_dir),
        base_dir=str(tmp_path),
        output_dir="GenBank-matrix",
        ref_list=str(ref_file),
        exclusion_list=None,
        is_segmented_virus="N",
    )
    parser.process()

    matrix_path = tmp_path / "GenBank-matrix" / "gB_matrix_raw.tsv"
    fasta_path = tmp_path / "GenBank-matrix" / "sequences.fa"

    assert matrix_path.exists()
    assert fasta_path.exists()

    df = pd.read_csv(matrix_path, sep="\t")
    assert set(df["primary_accession"].tolist()) == {"REF1", "Q1"}
    assert "sequence" not in df.columns

    fasta_text = fasta_path.read_text(encoding="utf-8")
    assert ">REF1" in fasta_text
    assert ">Q1" in fasta_text
