#!/usr/bin/env python3
import os
import re
import pandas as pd
from time import sleep
from os.path import join
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from date_utils import split_date_components

"""
Running with update mode:
	python GenBankParser.py -b tmp -o GenBank-matrix -r ./../../RABV-gTK/generic/rabv/ref_list.txt -e ./../../RABV-gTK/generic/rabv/exclusion_list.txt --update --update_accession tmp/update/update_accession.tsv

No update mode:
	python GenBankParser.py -b tmp -o GenBank-matrix -r ./../../RABV-gTK/generic/rabv/ref_list.txt -e ./../../RABV-gTK/generic/rabv/exclusion_list.txt
"""

ACC_VERSION_RE = re.compile(r"^(.+?)\.(\d+)$")


def norm_acc(acc: str) -> str:
    """Normalize accession to unversioned, uppercase."""
    if acc is None:
        return ""
    acc = str(acc).strip()
    if not acc:
        return ""
    acc = acc.split()[0].upper()
    m = ACC_VERSION_RE.match(acc)
    if m:
        acc = m.group(1)
    return acc


class GenBankParser:
    def __init__(
        self,
        base_dir,
        output_dir,
        ref_list,
        exclusion_list,
        is_segmented_virus,
        update_mode=False,
        update_accessions_tsv=None,
        input_dir=None,
    ):
        # base_dir is tmp-like root (e.g. "tmp")
        # output_dir is subdir name (e.g. "GenBank-matrix")
        self.base_dir = base_dir
        self.output_dir = output_dir

        self.ref_list = ref_list
        self.exclusion_list = exclusion_list
        self.is_segmented_virus = is_segmented_virus

        self.update_mode = update_mode
        self.update_accessions_tsv = update_accessions_tsv

        # Input XML directory:
        # normal: <base_dir>/GenBank-XML
        # update : <base_dir>/update/GenBank-XML
        if input_dir:
            self.input_dir = input_dir
        else:
            self.input_dir = (
                join(self.base_dir, "update", "GenBank-XML")
                if self.update_mode
                else join(self.base_dir, "GenBank-XML")
            )

        # Output root directory:
        # normal: <base_dir>/<output_dir>
        # update : <base_dir>/update/<output_dir>
        self.out_root = (
            join(self.base_dir, "update", self.output_dir)
            if self.update_mode
            else join(self.base_dir, self.output_dir)
        )
        os.makedirs(self.out_root, exist_ok=True)

        # populated in process() when update_mode
        self._update_accver_set = None

    def count_ATGCN(self, sequence):
        nucl_dict = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        sequence = (sequence or "").upper()
        for each_nucl in sequence:
            if each_nucl in nucl_dict:
                nucl_dict[each_nucl] += 1
        return nucl_dict["A"], nucl_dict["T"], nucl_dict["G"], nucl_dict["C"], nucl_dict["N"]

    def load_ref_list(self, acc_list_file):
        if acc_list_file is None:
            return {}
        ref_dict = {}
        try:
            with open(acc_list_file) as file:
                for line in file:
                    if not line.strip():
                        continue
                    if self.is_segmented_virus == "Y":
                        accession, segment_type, accession_type = line.strip().split("\t")
                    else:
                        accession, accession_type = line.strip().split("\t")
                    if accession not in ref_dict:
                        ref_dict[accession] = accession_type
        except FileNotFoundError:
            print(f"Error!!!: File {acc_list_file} not found. Exiting the program")
            raise SystemExit(1)
        return ref_dict

    def load_exclusion_list(self, acc_list_file):
        if acc_list_file is None:
            return []
        ref_list = []
        try:
            with open(acc_list_file) as file:
                for line in file:
                    accession = line.strip()
                    if accession and accession not in ref_list:
                        ref_list.append(accession)
        except FileNotFoundError:
            print(f"Warning: File {acc_list_file} not found. Proceeding with all the available sequences")
            sleep(2)
        return ref_list

    def load_update_accessions(self, update_tsv_path):
        """
        Reads update_accessions.tsv produced by the fetcher.
        Expected columns: old_accession_version, new_accession_version
        We keep ONLY new_accession_version values (ignore NA).
        Used to filter XML records by accession_version EXACT match.
        """
        if not update_tsv_path:
            raise ValueError("update_accessions_tsv path not provided")
        if not os.path.exists(update_tsv_path):
            raise FileNotFoundError(f"update_accessions.tsv not found: {update_tsv_path}")

        df = pd.read_csv(update_tsv_path, sep="\t", dtype=str).fillna("")
        if "new_accession_version" not in df.columns:
            raise ValueError("update_accessions.tsv must contain column: new_accession_version")

        new_accs = set(
            x.strip()
            for x in df["new_accession_version"].tolist()
            if x and x.strip() and x.strip().upper() != "NA"
        )
        print(f"[update] Loaded {len(new_accs)} new_accession_version values from {update_tsv_path}")
        print("[update] first 10:", list(sorted(new_accs))[:10])
        return new_accs

    def xml_to_rows(self, xml_file, ref_seq_dict: dict, exclusion_acc_list: list):
        tree = ET.parse(xml_file)
        root = tree.getroot()
        rows = []

        for gbseq in root.findall("GBSeq"):
            content = {}
            content["locus"] = gbseq.findtext("GBSeq_locus")
            content["length"] = gbseq.findtext("GBSeq_length")
            content["data_source"] = "ncbi"
            content["strandedness"] = gbseq.findtext("GBSeq_strandedness")
            content["molecule_type"] = gbseq.findtext("GBSeq_moltype")
            content["topology"] = gbseq.findtext("GBSeq_topology")
            content["division"] = gbseq.findtext("GBSeq_division")
            content["update_date"] = gbseq.findtext("GBSeq_update-date")
            content["create_date"] = gbseq.findtext("GBSeq_create-date")
            content["definition"] = gbseq.findtext("GBSeq_definition")
            content["primary_accession"] = gbseq.findtext("GBSeq_primary-accession")
            content["accession_version"] = gbseq.findtext("GBSeq_accession-version")
            content["gi_number"] = content["primary_accession"]
            content["source"] = gbseq.findtext("GBSeq_source")
            content["organism"] = gbseq.findtext("GBSeq_organism")
            content["taxonomy"] = gbseq.findtext("GBSeq_taxonomy")

            # UPDATE FILTER: only keep records listed in update_accessions.tsv
            if self.update_mode and self._update_accver_set is not None:
                if content["accession_version"] not in self._update_accver_set:
                    continue

            # accession_type
            if content["primary_accession"] in ref_seq_dict:
                content["accession_type"] = ref_seq_dict[content["primary_accession"]]
            else:
                content["accession_type"] = "query"

            # exclusion
            if content["primary_accession"] in exclusion_acc_list:
                content["accession_type"] = "excluded"
                content["exclusion_criteria"] = "excluded by the user"
                content["exclusion_status"] = "1"
            else:
                content["exclusion_criteria"] = ""
                content["exclusion_status"] = "0"

            mol_type = ""
            strain = ""
            isolate = ""
            isolation_source = ""
            db_xref = ""
            country = ""
            host = ""
            collection_date = ""
            segment = ""
            serotype = ""
            genes = []
            cds = []

            for gb_feature in gbseq.findall("GBSeq_feature-table/GBFeature"):
                key = gb_feature.findtext("GBFeature_key")
                if key == "source":
                    for qualifier in gb_feature.findall("GBFeature_quals/GBQualifier"):
                        name = qualifier.findtext("GBQualifier_name")
                        value = qualifier.findtext("GBQualifier_value")
                        if name == "mol_type":
                            mol_type = value or ""
                        elif name == "strain":
                            strain = value or ""
                        elif name == "isolate":
                            isolate = value or ""
                        elif name == "isolation_source":
                            isolation_source = value or ""
                        elif name == "db_xref":
                            db_xref = value or ""
                        elif name in ("country", "geo_loc_name"):
                            country = value or ""
                        elif name == "host":
                            host = value or ""
                        elif name == "collection_date":
                            collection_date = value or ""
                        elif name == "segment":
                            segment = value or ""
                        elif name == "serotype":
                            serotype = value or ""
                elif key == "gene":
                    gene_info = {"gene_location": gb_feature.findtext("GBFeature_location")}
                    for qualifier in gb_feature.findall("GBFeature_quals/GBQualifier"):
                        name = qualifier.findtext("GBQualifier_name")
                        value = qualifier.findtext("GBQualifier_value")
                        if name == "gene":
                            gene_info["gene_name"] = value
                            genes.append(gene_info)
                elif key == "CDS":
                    cds_info = {"cds_location": gb_feature.findtext("GBFeature_location")}
                    for qualifier in gb_feature.findall("GBFeature_quals/GBQualifier"):
                        name = qualifier.findtext("GBQualifier_name")
                        value = qualifier.findtext("GBQualifier_value")
                        if name:
                            cds_info[name] = value
                    cds.append(cds_info)

            pubmed_ids = []
            for reference in gbseq.findall("GBSeq_references/GBReference"):
                pubmed_tag = reference.find("GBReference_pubmed")
                if pubmed_tag is not None and pubmed_tag.text:
                    pubmed_ids.append(pubmed_tag.text)

            content["pubmed_id"] = "; ".join(pubmed_ids)
            content["mol_type"] = mol_type
            content["strain"] = strain
            content["isolate"] = isolate
            content["isolation_source"] = isolation_source
            content["db_xref"] = db_xref

            if country and ":" in country:
                tmp_country = country.split(":", 1)
                content["country"] = tmp_country[0]
                content["geo_loc"] = tmp_country[1] if len(tmp_country) > 1 else ""
            else:
                content["country"] = country
                content["geo_loc"] = ""

            content["host"] = host
            content["collection_date"] = collection_date

            split_collection_date = split_date_components(content["collection_date"])
            content["collection_day"] = split_collection_date.get("day")
            content["collection_mon"] = split_collection_date.get("month")
            content["collection_year"] = split_collection_date.get("year")

            content["segment"] = segment
            content["serotype"] = serotype

            sequence = gbseq.find("GBSeq_sequence")
            content["sequence"] = sequence.text if sequence is not None else ""

            content["genes"] = "; ".join(
                [f"{g.get('gene_name','')}({g.get('gene_location','')})" for g in genes if g.get("gene_name")]
            )
            content["cds_info"] = "; ".join([f"{k}: {v}" for each_cds in cds for k, v in each_cds.items()])

            a, t, g, c, n = self.count_ATGCN(content["sequence"])
            content["a"] = a
            content["t"] = t
            content["g"] = g
            content["c"] = c
            content["n"] = n
            content["real_length"] = int(a) + int(t) + int(g) + int(c)
            content["comment"] = "NA"

            rows.append(content)

        return rows

    def count_xml_files(self):
        if not os.path.isdir(self.input_dir):
            print(f"Error: '{self.input_dir}' is not a valid directory.")
            return 0
        return len([f for f in os.listdir(self.input_dir) if f.lower().endswith(".xml")])

    def process(self):
        ref_seq_dict = self.load_ref_list(self.ref_list)
        exclusion_acc_list = self.load_exclusion_list(self.exclusion_list)

        if self.update_mode:
            update_path = self.update_accessions_tsv
            if not update_path:
                update_path = join(self.base_dir, "update", "update_accessions.tsv")
            self._update_accver_set = self.load_update_accessions(update_path)

        total_xml = self.count_xml_files()
        if total_xml == 0:
            print(f"No XML files found in: {self.input_dir}")
            return

        merged_rows = []
        count = 1
        for each_xml in os.listdir(self.input_dir):
            if not each_xml.lower().endswith(".xml"):
                continue
            print(f"Parsing: {count} of {total_xml}: {each_xml}")
            rows = self.xml_to_rows(join(self.input_dir, each_xml), ref_seq_dict, exclusion_acc_list)
            merged_rows.extend(rows)
            count += 1

        if not merged_rows:
            if self.update_mode:
                print("[update] No matching accession_version records found in XMLs (nothing to write).")
            else:
                print("No records parsed.")
            return

        df = pd.DataFrame(merged_rows)

        out_matrix = join(self.out_root, "gB_matrix_raw.tsv")
        out_fasta = join(self.out_root, "sequences.fa")

        df_dedup = df.drop_duplicates(subset="locus", keep="last")

        with open(out_fasta, "w") as fasta_file:
            for _, row in df_dedup.iterrows():
                fasta_file.write(f">{row['primary_accession']}\n{row['sequence']}\n")

        df_out = df_dedup.drop(columns=["sequence"], errors="ignore")
        df_out.to_csv(out_matrix, sep="\t", index=False)

        print(f"Wrote FASTA: {out_fasta}")
        print(f"Wrote TSV:   {out_matrix}")
        print(f"Input XML dir used: {self.input_dir}")
        print(f"Output dir used:    {self.out_root}")


if __name__ == "__main__":
    ap = ArgumentParser(description="Extract GenBank XML files to a TSV table + FASTA")
    ap.add_argument("-b", "--base_dir", help="Base directory (tmp-like)", default="tmp")
    ap.add_argument("-o", "--output_dir", help="Output directory name under base_dir", default="GenBank-matrix")
    ap.add_argument("-r", "--ref_list", help="Set of reference accessions", required=True)
    ap.add_argument("-e", "--exclusion_list", help="Set of sequence accessions to be excluded")
    ap.add_argument("-s", "--is_segmented_virus", help="Is segmented virus (Y/N)", default="N")

    # optional override (otherwise auto-chosen based on --update)
    ap.add_argument("-d", "--input_dir", help="Input directory (optional override)", default=None)

    # update mode
    ap.add_argument(
        "--update",
        action="store_true",
        help="If set, read XMLs from <base_dir>/update/GenBank-XML and write outputs to <base_dir>/update/<output_dir>.",
    )
    ap.add_argument(
        "--update_accessions",
        default=None,
        help="Path to update_accessions.tsv (default: <base_dir>/update/update_accessions.tsv)",
    )

    args = ap.parse_args()

    parser = GenBankParser(
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        ref_list=args.ref_list,
        exclusion_list=args.exclusion_list,
        is_segmented_virus=args.is_segmented_virus,
        update_mode=args.update,
        update_accessions_tsv=args.update_accessions,
        input_dir=args.input_dir,
    )
    parser.process()
