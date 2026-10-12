"""Regression: a staged input that is a broken symlink must say so.

A bare `nextflow run -resume` resumes the *last session*, not the last run of
the current profile. SOFTWARE_VERSION took no inputs, so it hashed identically
in every profile and a disposable `-profile segmented_test -w work_test_seg`
run could satisfy a production HCV run. When that test work dir was later
deleted, CREATE_SQLITE_DB was left staging software_info.tsv through a dangling
symlink and died a day into the run with a bare "file not found" for a name
that was visibly present in the task directory.
"""

import os
from pathlib import Path

import pytest

from CreateSqliteDB import CreateSqliteDB


def test_broken_staged_symlink_names_its_missing_target(tmp_path: Path):
    missing_source = tmp_path / "work_test_seg" / "c3" / "Software_info" / "software_info.tsv"
    staged = tmp_path / "software_info.tsv"
    os.symlink(missing_source, staged)

    assert staged.is_symlink() and not staged.exists()

    with pytest.raises(FileNotFoundError) as excinfo:
        CreateSqliteDB._require_file(str(staged), "proj_settings")

    message = str(excinfo.value)
    assert "broken symlink" in message
    assert str(missing_source) in message


def test_intact_symlink_is_accepted(tmp_path: Path):
    real = tmp_path / "real.tsv"
    real.write_text("software\tversion\n")
    staged = tmp_path / "software_info.tsv"
    os.symlink(real, staged)

    CreateSqliteDB._require_file(str(staged), "proj_settings")


def test_plain_missing_file_keeps_the_original_message(tmp_path: Path):
    absent = tmp_path / "nope.tsv"

    with pytest.raises(FileNotFoundError) as excinfo:
        CreateSqliteDB._require_file(str(absent), "proj_settings")

    message = str(excinfo.value)
    assert "proj_settings file not found" in message
    assert "broken symlink" not in message


@pytest.mark.parametrize("empty", [None, ""])
def test_unset_path_still_raises(empty):
    with pytest.raises(FileNotFoundError):
        CreateSqliteDB._require_file(empty, "proj_settings")


def test_no_cached_process_is_without_inputs():
    """Root-cause guard: an input-less process hashes the same in every profile,
    so -resume can satisfy it from another profile's (possibly deleted) work dir.
    Every process must take an input or opt out of caching."""
    import re

    nf = Path(__file__).resolve().parents[2] / "vgtk-init.nf"
    blocks = re.split(r"(?m)^process\s+", nf.read_text())[1:]
    offenders = []
    for block in blocks:
        name = re.match(r"\w+", block).group(0)
        header = re.split(r"(?m)^\s*(?:shell|script|exec):", block)[0]
        if not re.search(r"(?m)^\s*input:", header) and not re.search(r"(?m)^\s*cache\s+false", header):
            offenders.append(name)
    assert offenders == [], f"processes with no inputs and caching on: {offenders}"
