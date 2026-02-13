import csv
from pathlib import Path

from Bio import SeqIO

from MMseqsClustering import (
    apply_trimming,
    parse_ranges,
    run_mmseqs_clustering,
    trim_alignment,
)


def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def test_parse_ranges_with_open_bounds():
    # alignment length is used for open-ended ranges
    assert parse_ranges(["1:3", "5:"], alignment_length=8) == [(0, 3), (4, 8)]
    assert parse_ranges([":2"], alignment_length=8) == [(0, 2)]


def test_trim_alignment_and_apply_trimming(tmp_path: Path):
    input_dir = tmp_path / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    in_fasta = input_dir / "segA.fasta"
    write_fasta(in_fasta, [("S1", "ABCDEFGH"), ("S2", "12345678")])

    out_fasta = tmp_path / "trimmed.fasta"
    trim_alignment(str(in_fasta), str(out_fasta), ["1:2", "5:6"])

    recs = list(SeqIO.parse(str(out_fasta), "fasta"))
    assert str(recs[0].seq) == "ABEF"
    assert str(recs[1].seq) == "1256"

    cds = tmp_path / "cds.tsv"
    with cds.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["input_fasta", "ranges"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"input_fasta": "segA.fasta", "ranges": "2:3,7:8"})

    trimmed_map = apply_trimming(str(cds), str(input_dir), str(tmp_path / "trimmed_dir"))
    assert "segA.fasta" in trimmed_map

    trimmed_records = list(SeqIO.parse(trimmed_map["segA.fasta"], "fasta"))
    assert str(trimmed_records[0].seq) == "BCGH"


def test_run_mmseqs_clustering_builds_expected_command_chain(tmp_path: Path, monkeypatch):
    in_fasta = tmp_path / "input.fasta"
    write_fasta(in_fasta, [("S1", "AAAA"), ("S2", "AAAT")])

    calls = []

    def fake_run(cmd, check=True):
        calls.append(cmd)
        return 0

    monkeypatch.setattr("MMseqsClustering.subprocess.run", fake_run)

    run_mmseqs_clustering(str(in_fasta), str(tmp_path / "out"), min_seq_id=0.9, threads=3)

    assert len(calls) == 7
    assert calls[0][:2] == ["mmseqs", "createdb"]
    assert calls[1][:2] == ["mmseqs", "cluster"]
    assert "--min-seq-id" in calls[1]
    assert "--threads" in calls[1]
    assert calls[2][:2] == ["mmseqs", "createtsv"]
    assert calls[3][:2] == ["mmseqs", "createseqfiledb"]
    assert calls[4][:2] == ["mmseqs", "result2flat"]
    assert calls[5][:2] == ["mmseqs", "createsubdb"]
    assert calls[6][:2] == ["mmseqs", "convert2fasta"]
