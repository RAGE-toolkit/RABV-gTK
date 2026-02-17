#!/usr/bin/env python3
import os
import csv
import time
import sqlite3
import requests
from time import sleep
from os.path import join
from argparse import ArgumentParser

'''
Command:
	python GenBankFetcher.py -t 11292 --update  --gb_db ./../tmp/SqliteDB/rabv-gDB.db --update_log update_accession.tsv --update --meta_acc_col accession_version --sleep_time 2 --batch_size 100000 --tmp_dir tmp --base_dir GenBank-XML


'''

class GenBankFetcher:
    def __init__(
        self,
        taxid,
        base_url,
        email,
        output_dir,
        batch_size,
        sleep_time,
        base_dir,
        update_log,
    ):
        self.taxid = taxid
        self.base_url = base_url
        self.email = email
        self.output_dir = output_dir
        self.batch_size = batch_size
        self.efetch_batch_size = 100
        self.sleep_time = sleep_time
        self.base_dir = base_dir
        self.update_log = update_log

        # internal flag to route downloads to tmp/update when in update mode
        self._in_update_mode = False

    # -----------------------------
    # Output routing (normal vs update)
    # -----------------------------
    def _ensure_update_output_dirs(self):
        update_root = join(self.output_dir, "update")
        os.makedirs(update_root, exist_ok=True)
        os.makedirs(join(update_root, self.base_dir), exist_ok=True)
        return update_root

    def _xml_root_dir(self):
        if self._in_update_mode:
            return join(self.output_dir, "update", self.base_dir)
        return join(self.output_dir, self.base_dir)

    def _log_path(self):
        if self._in_update_mode:
            os.makedirs(join(self.output_dir, "update"), exist_ok=True)
            return join(self.output_dir, "update", self.update_log)
        return join(self.output_dir, self.update_log)

    # -----------------------------
    # NCBI helpers
    # -----------------------------
    def get_record_count(self):
        search_url = (
            f"{self.base_url}esearch.fcgi?db=nucleotide"
            f"&term=txid{self.taxid}[Organism:exp]"
            f"&retmode=json&email={self.email}"
        )
        response = requests.get(search_url, timeout=60)
        response.raise_for_status()
        data = response.json()
        return int(data["esearchresult"]["count"])

    def fetch_ids(self):
        retmax = self.get_record_count()
        search_url = (
            f"{self.base_url}esearch.fcgi?db=nucleotide"
            f"&term=txid{self.taxid}[Organism:exp]"
            f"&retmax={retmax}&idtype=acc"
            f"&usehistory=y&email={self.email}&retmode=json"
        )
        response = requests.get(search_url, timeout=60)
        response.raise_for_status()
        data = response.json()
        return data["esearchresult"]["idlist"]

    def fetch_accs(self):
        start_all = time.time()

        t0 = time.time()
        hist_url = (
            f"{self.base_url}esearch.fcgi?db=nucleotide"
            f"&term=txid{self.taxid}[Organism:exp]"
            f"&retmax=0&idtype=acc"
            f"&usehistory=y&email={self.email}&retmode=json"
        )
        hist = requests.get(hist_url, timeout=60)
        hist.raise_for_status()
        hist = hist.json()["esearchresult"]

        webenv = hist["webenv"]
        querykey = hist["querykey"]
        count = int(hist["count"])
        print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

        all_accs = []
        for start in range(0, count, self.batch_size):
            t1 = time.time()
            page_url = (
                f"{self.base_url}esearch.fcgi?db=nucleotide"
                f"&WebEnv={webenv}&query_key={querykey}"
                f"&retstart={start}&retmax={self.batch_size}"
                f"&idtype=acc&retmode=json"
            )
            page_resp = requests.get(page_url, timeout=60)
            page_resp.raise_for_status()
            page = page_resp.json()["esearchresult"]["idlist"]

            all_accs.extend(page)
            print(
                f"[fetch_accs] chunk {start:,}-{min(start+self.batch_size, count):,} "
                f"({len(page):,} records) took {time.time() - t1:.1f}s"
            )
            sleep(self.sleep_time)

        print(f"[fetch_accs] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
        return all_accs

    # -----------------------------
    # EFetch + file IO
    # -----------------------------
    def fetch_genbank_data(self, ids):
        batch_n = self.efetch_batch_size
        for i in range(0, len(ids), batch_n):
            chunk = ids[i : i + batch_n]
            ids_str = ",".join(chunk)
            url = (
                f"{self.base_url}efetch.fcgi?db=nucleotide"
                f"&id={ids_str}"
                f"&retmode=xml&email={self.email}"
            )
            resp = requests.get(url, timeout=120)
            resp.raise_for_status()

            self.save_data(resp.text, i + batch_n)
            print(f" - Downloaded XML {i+1:,}–{min(i+batch_n, len(ids)):,}")
            sleep(self.sleep_time)

    def save_data(self, data, batch_size):
        xml_dir = self._xml_root_dir()
        os.makedirs(xml_dir, exist_ok=True)

        # NEW: file prefix changes in update mode
        prefix = "batch-update" if self._in_update_mode else "batch"

        base_filename = join(xml_dir, f"{prefix}-{batch_size}")
        filename = f"{base_filename}.xml"
        counter = 1
        while os.path.exists(filename):
            filename = f"{base_filename}_{counter}.xml"
            counter += 1

        with open(filename, "w", encoding="utf-8") as file:
            file.write(data)
        print(f"Data written to: {filename}")

    # -----------------------------
    # DB helpers (default accession column + fallbacks)
    # -----------------------------
    def _detect_meta_data_acc_col(self, conn, preferred="accession_version"):
        cur = conn.cursor()
        cur.execute("PRAGMA table_info(meta_data)")
        cols = [row[1] for row in cur.fetchall()]

        if preferred in cols:
            return preferred

        fallbacks = ["primary_accession", "locus", "accession"]
        for c in fallbacks:
            if c in cols:
                print(f"[db] preferred column {preferred!r} not found; using {c!r} instead")
                return c

        raise ValueError(
            "Could not find an accession column in meta_data. "
            f"Tried preferred={preferred!r} and fallbacks={fallbacks}. Found columns={cols}"
        )

    def _split_accession_version(self, accession_version):
        if not accession_version:
            return (None, None)
        av = str(accession_version).strip()
        if not av:
            return (None, None)

        base, sep, ver = av.rpartition(".")
        if not sep:
            return (av, None)
        try:
            return (base, int(ver))
        except ValueError:
            return (base, None)

    # -----------------------------
    # Modes
    # -----------------------------
    def download(self):
        self._in_update_mode = False
        os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

        ids = self.fetch_ids()
        print(f"Found {len(ids)} IDs")
        self.fetch_genbank_data(ids)

    def update_from_tsv(self, update_tsv_path):
        self._in_update_mode = True
        self._ensure_update_output_dirs()

        with open(update_tsv_path, "r", encoding="utf-8") as file:
            reader = csv.DictReader(file, delimiter="\t")
            if "accession_version" not in reader.fieldnames:
                raise ValueError("Expecting a column called 'accession_version'")
            accession_versions = [row["accession_version"] for row in reader]

        print("first 10 accession_version in TSV:", accession_versions[:10])
        self._update_from_accession_versions(accession_versions)

    def update_from_db(self, db_path, meta_acc_col="accession_version"):
        self._in_update_mode = True
        self._ensure_update_output_dirs()

        if not os.path.exists(db_path):
            raise FileNotFoundError(f"SQLite DB not found: {db_path}")

        conn = sqlite3.connect(db_path)
        try:
            acc_col = self._detect_meta_data_acc_col(conn, preferred=meta_acc_col)
            cursor = conn.cursor()
            cursor.execute(f"SELECT {acc_col} FROM meta_data WHERE {acc_col} IS NOT NULL")
            accession_versions = [row[0] for row in cursor.fetchall()]
        finally:
            conn.close()

        print(f"[db] using meta_data column: {acc_col}")
        print("first 10 values in DB:", accession_versions[:10])
        self._update_from_accession_versions(accession_versions)

    # -----------------------------
    # Update logic
    # -----------------------------
    def _write_update_log(self, updated_versions, new_accessions):
        if not self.update_log:
            return

        log_path = self._log_path()
        wrote_any = False

        with open(log_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["old_accession_version", "new_accession_version"])
            for old_acc, new_acc in updated_versions:
                w.writerow([old_acc, new_acc])
                wrote_any = True
            for new_acc in new_accessions:
                w.writerow(["NA", new_acc])
                wrote_any = True

        if wrote_any:
            print(f"[update] Wrote update log: {log_path}")
        else:
            print(f"[update] No updates/new accessions; log created (header only): {log_path}")

    def _update_from_accession_versions(self, accession_versions):
        print(f"Found {len(accession_versions)} local accessions (from DB/TSV)")
        ids = self.fetch_accs()
        print("first 10 IDs in NCBI:", ids[:10])

        local_set = set(str(x).strip() for x in accession_versions if x is not None)

        # For each base accession, store max local version
        local_versions = {}
        for acc_ver in accession_versions:
            base, ver = self._split_accession_version(acc_ver)
            if not base:
                continue
            if ver is None:
                local_versions.setdefault(base, 0)
                continue
            prev = local_versions.get(base)
            if prev is None or ver > prev:
                local_versions[base] = ver

        missing_ids = []
        updated_versions = []
        new_accessions = []

        for acc_ver in ids:
            acc_ver = str(acc_ver).strip()
            if not acc_ver:
                continue

            if acc_ver in local_set:
                continue

            base, ver = self._split_accession_version(acc_ver)

            if base in local_versions and ver is not None and ver > local_versions[base]:
                old_acc_ver = f"{base}.{local_versions[base]}" if local_versions[base] > 0 else base
                updated_versions.append((old_acc_ver, acc_ver))
                missing_ids.append(acc_ver)
            elif base not in local_versions:
                new_accessions.append(acc_ver)
                missing_ids.append(acc_ver)

        if updated_versions:
            print("Updated accessions found (old -> new):")
            for old_acc_ver, new_acc_ver in updated_versions[:50]:
                print(f"  {old_acc_ver} -> {new_acc_ver}")
            if len(updated_versions) > 50:
                print(f"  ... ({len(updated_versions)-50} more)")

        print(
            f"Found {len(missing_ids)} missing/new IDs total "
            f"(updates={len(updated_versions)}, brand-new={len(new_accessions)})"
        )

        self._write_update_log(updated_versions, new_accessions)

        if missing_ids:
            self.fetch_genbank_data(missing_ids)
        else:
            print("[update] Nothing to download.")


def build_argparser():
    p = ArgumentParser(description="Download and update GenBank XML files for a given TaxID")

    p.add_argument("-t", "--taxid", help="TaxID example: 11292", required=True)
    p.add_argument("-o", "--tmp_dir", help="Output directory where XML files are stored", default="tmp")
    p.add_argument("-b", "--batch_size", help="Max number of accessions per esearch page", default=100000, type=int)

    p.add_argument(
        "--update",
        action="store_true",
        help="Enable update mode. Outputs go under tmp/update/ (XML in tmp/update/GenBank-XML, logs in tmp/update/).",
    )
    p.add_argument("--gb_db", help="SQLite DB file containing meta_data table (preferred update source)")
    p.add_argument("--update_tsv", help="TSV file for update mode (must have column accession_version)")

    p.add_argument(
        "--meta_acc_col",
        default="accession_version",
        help="Column name in meta_data to use as accession source (default: accession_version).",
    )

    p.add_argument("-u", "--base_url", help="Base URL", default="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
    p.add_argument("-e", "--email", help="Email ID", default="your_email@example.com")
    p.add_argument("-s", "--sleep_time", help="Delay after each efetch batch", default=2, type=int)
    p.add_argument("-d", "--base_dir", help="Directory name for XML storage", default="GenBank-XML")
    p.add_argument("--update_log", help="TSV log of updated/new accessions", default="update_accessions.tsv")
    return p


if __name__ == "__main__":
    parser = build_argparser()
    args = parser.parse_args()

    fetcher = GenBankFetcher(
        taxid=args.taxid,
        base_url=args.base_url,
        email=args.email,
        output_dir=args.tmp_dir,
        batch_size=args.batch_size,
        sleep_time=args.sleep_time,
        base_dir=args.base_dir,
        update_log=args.update_log,
    )

    if args.update:
        if args.gb_db:
            fetcher.update_from_db(args.gb_db, meta_acc_col=args.meta_acc_col)
        elif args.update_tsv:
            fetcher.update_from_tsv(args.update_tsv)
        else:
            raise SystemExit("Update mode requires either --gb_db or --update_tsv")
    else:
        fetcher.download()
