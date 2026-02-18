import csv
from pathlib import Path

import pytest

from GenBankFetcher import GenBankFetcher


class _FakeResponse:
    def __init__(self, json_data=None, text_data="", status_code=200):
        self._json_data = json_data
        self.text = text_data
        self.status_code = status_code

    def json(self):
        return self._json_data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise Exception(f"HTTP {self.status_code}")


def test_fetch_ids_paginates_and_strips_versions(tmp_path: Path, monkeypatch):
    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=2,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    def fake_get(url):
        if "retmax=0" in url:
            return _FakeResponse({"esearchresult": {"count": "3", "webenv": "W", "querykey": "1"}})
        if "retstart=0" in url:
            return _FakeResponse({"esearchresult": {"idlist": ["A.1", "B.2"]}})
        if "retstart=2" in url:
            return _FakeResponse({"esearchresult": {"idlist": ["C.9"]}})
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    ids = fetcher.fetch_ids()
    assert ids == ["A", "B", "C"]


def test_fetch_genbank_data_adds_ref_list_and_saves_batches(tmp_path: Path, monkeypatch):
    ref_file = tmp_path / "refs.tsv"
    ref_file.write_text("REF1\tmaster\n", encoding="utf-8")

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
        ref_list=str(ref_file),
    )
    fetcher.efetch_batch_size = 2

    calls = []

    def fake_get(url):
        calls.append(url)
        return _FakeResponse(json_data={}, text_data="<GBSet></GBSet>")

    saved = []

    def fake_save(data, marker):
        saved.append((data, marker))

    monkeypatch.setattr("GenBankFetcher.requests.get", fake_get)
    monkeypatch.setattr(fetcher, "save_data", fake_save)

    fetcher.fetch_genbank_data(["A", "A", "B"])

    assert len(saved) >= 1
    assert all(data == "<GBSet></GBSet>" for data, _ in saved)
    assert any("efetch.fcgi" in c for c in calls)


def test_update_requires_primary_accession_column(tmp_path: Path):
    bad_tsv = tmp_path / "bad.tsv"
    with bad_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["wrong_col"])
        writer.writerow(["A"])

    fetcher = GenBankFetcher(
        taxid="11292",
        base_url="https://example/",
        email="x@y.com",
        output_dir=str(tmp_path),
        batch_size=100,
        sleep_time=0,
        base_dir="GenBank-XML",
        update_file=None,
    )

    with pytest.raises(ValueError, match="primary_accession"):
        fetcher.update(str(bad_tsv))
