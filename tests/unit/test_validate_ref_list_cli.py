import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ValidateRefList.py"


def test_validate_ref_list_non_segmented_ok(tmp_path: Path):
    ref = tmp_path / "ref.tsv"
    ref.write_text("ACC1\tmaster\nACC2\treference\nACC3\texclusion_list\n", encoding="utf-8")

    out = tmp_path / "validated.txt"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "-r",
            str(ref),
            "-s",
            "N",
            "-o",
            str(out),
        ],
        check=True,
    )

    assert out.read_text(encoding="utf-8").strip() == "Reference list validation passed."


def test_validate_ref_list_segmented_missing_column_fails(tmp_path: Path):
    ref = tmp_path / "bad.tsv"
    ref.write_text("ACC1\tmaster\nACC2\treference\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "reference list must contain three columns" in result.stderr.lower()


def test_validate_ref_list_invalid_indicator_fails(tmp_path: Path):
    ref = tmp_path / "bad_indicator.tsv"
    ref.write_text("ACC1\tmaster\t1\nACC2\tbadvalue\t1\n", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "invalid values" in result.stderr.lower()


def test_validate_ref_list_exclusion_only_segment_is_allowed(tmp_path: Path):
    ref = tmp_path / "ref.tsv"
    ref.write_text(
        "ACC1\tmaster\t1\n"
        "ACC2\treference\t1\n"
        "ACC3\texclusion_list\t9\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "-r", str(ref), "-s", "Y"],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0
    assert "validation passed" in result.stdout.lower()
