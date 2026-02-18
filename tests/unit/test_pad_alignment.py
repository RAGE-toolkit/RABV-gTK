import shutil
from pathlib import Path

from Bio import SeqIO

from PadAlignment import PadAlignment


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "test_data" / "unit" / "pad_alignment"


def _copy_fixture_tree(tmp_path: Path):
    input_dir = tmp_path / "query_aln"
    (input_dir / "REF_OK").mkdir(parents=True, exist_ok=True)
    (input_dir / "REF_ORPHAN").mkdir(parents=True, exist_ok=True)

    shutil.copyfile(DATA_DIR / "REF_OK.aligned.fasta", input_dir / "REF_OK" / "REF_OK.aligned.fasta")
    shutil.copyfile(
        DATA_DIR / "REF_ORPHAN.aligned.fasta",
        input_dir / "REF_ORPHAN" / "REF_ORPHAN.aligned.fasta",
    )

    return input_dir


def _make_processor(tmp_path: Path):
    return PadAlignment(
        reference_alignment=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(tmp_path / "query_aln"),
        base_dir=str(tmp_path),
        output_dir="pad_out",
        keep_intermediate_files=True,
        new_outputfile=False,
    )


def test_find_orphan_query_references_identifies_unprojectable_refs(tmp_path: Path):
    input_dir = _copy_fixture_tree(tmp_path)
    processor = _make_processor(tmp_path)

    orphan_refs = processor.find_orphan_query_references(
        str(DATA_DIR / "master_ref.aligned.fasta"),
        str(input_dir),
    )

    assert set(orphan_refs.keys()) == {"REF_ORPHAN"}
    assert sorted(orphan_refs["REF_ORPHAN"]) == ["Q_ORPHAN_1", "Q_ORPHAN_2"]


def test_process_master_alignment_warns_and_skips_orphan_refs(tmp_path: Path, capsys):
    input_dir = _copy_fixture_tree(tmp_path)
    processor = _make_processor(tmp_path)

    processor.process_master_alignment(
        reference_alignment_file=str(DATA_DIR / "master_ref.aligned.fasta"),
        input_dir=str(input_dir),
        base_dir=str(tmp_path),
        output_dir=str(tmp_path / "pad_out"),
        keep_intermediate_files=True,
    )

    out = capsys.readouterr().out
    assert "[warn]" in out
    assert "REF_ORPHAN" in out

    merged_path = tmp_path / "pad_out" / "master_ref.aligned_merged_MSA.fasta"
    assert merged_path.exists()

    merged_ids = {r.id for r in SeqIO.parse(str(merged_path), "fasta")}
    assert "Q_OK_1" in merged_ids
    assert "Q_OK_2" in merged_ids
    assert "Q_ORPHAN_1" not in merged_ids
    assert "Q_ORPHAN_2" not in merged_ids
