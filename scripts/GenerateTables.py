import os
import csv
import time
import urllib.error
import read_file
from Bio import Entrez
from os.path import join
from itertools import islice
from argparse import ArgumentParser

'''
	python GenerateTables.py -g tmp/update/GenBank-matrix/gB_matrix_raw.tsv -bh tmp/update/Blast/query_uniq_tophits.tsv -p tmp/update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa -n tmp/update/Nextalign --update
'''

class GenerateTables:
    def __init__(
        self,
        genbank_matrix,
        base_dir,
        output_dir,
        blast_hits,
        paded_aln,
        nextalign_dir,
        email,
        update=False,
    ):
        self.genbank_matrix = genbank_matrix
        self.base_dir = base_dir
        self.output_dir = output_dir  # e.g. "Tables"
        self.blast_hits = blast_hits
        self.paded_aln = paded_aln
        self.nextalign_dir = nextalign_dir
        self.email = email
        self.update = update

        # NEW:
        # Normal mode -> <base_dir>/<output_dir>   e.g. tmp/Tables
        # Update mode  -> <base_dir>/update/<output_dir>  e.g. tmp/update/Tables
        if self.update:
            self.effective_output_dir = join("update", self.output_dir)
        else:
            self.effective_output_dir = self.output_dir

        os.makedirs(join(self.base_dir, self.effective_output_dir), exist_ok=True)

    def fetch_taxonomy_details(self, tax_id, max_retries=5, delay=2):
        """
        Kept for potential future use, but not called currently.
        """
        Entrez.email = self.email
        for attempt in range(1, max_retries + 1):
            try:
                handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                records = Entrez.read(handle)
                time.sleep(1)
                handle.close()

                if records:
                    tax_record = records[0]
                    taxonomy_info = {
                        "Scientific Name": tax_record.get("ScientificName", "N/A"),
                        "Taxonomy ID": tax_record.get("TaxId", "N/A"),
                        "Rank": tax_record.get("Rank", "N/A"),
                        "Lineage": tax_record.get("Lineage", "N/A"),
                        "Other Names": tax_record.get("OtherNames", {}).get("Synonym", []),
                    }
                    return taxonomy_info
                return "No taxonomy details found for the given ID."

            except urllib.error.HTTPError as e:
                print(f"HTTPError on attempt {attempt} for TaxID {tax_id}: {e}")
                if attempt == max_retries:
                    print("Max retries reached. Skipping this TaxID.")
                    return {
                        "Scientific Name": "N/A",
                        "Taxonomy ID": tax_id,
                        "Rank": "N/A",
                        "Lineage": "N/A",
                        "Other Names": [],
                    }
                time.sleep(delay)

    @staticmethod
    def load_blast_hits(blast_hit_file):
        acc_dict = {}
        with open(blast_hit_file) as f:
            for line in f:
                query, ref, score, strand = line.strip().split("\t")
                acc_dict[query] = ref
        return acc_dict

    def created_alignment_table(self, blast_dict):
        """
        Writes:
          - Normal mode: <base_dir>/<output_dir>/sequence_alignment.tsv
          - Update mode : <base_dir>/update/<output_dir>/sequence_alignment.tsv
        """
        accessions = {}
        missing_accs = []
        seqs = {}

        out_path = join(self.base_dir, self.effective_output_dir, "sequence_alignment.tsv")

        header = ["sequence_id", "alignment_name", "alignment"]
        with open(out_path, "w") as write_file:
            write_file.write("\t".join(header) + "\n")

            # 1) From padded alignment
            rds = read_file.fasta(self.paded_aln)
            for rows in rds:
                seq_id = rows[0].strip()
                seq = rows[1]

                if seq_id in accessions:
                    continue

                seqs[seq_id] = seq
                accessions[seq_id] = 1

                if seq_id in blast_dict:
                    write_file.write(seq_id + "\t" + blast_dict[seq_id] + "\t" + seq + "\n")
                else:
                    missing_accs.append(seq_id)

            # 2) Fill missing ones from nextalign results
            for each_ref_aln in os.listdir(self.nextalign_dir):
                ref_dir = join(self.nextalign_dir, each_ref_aln)
                if not os.path.isdir(ref_dir):
                    continue

                for each_ref_aln_file in os.listdir(ref_dir):
                    aln_dir = join(ref_dir, each_ref_aln_file)
                    if not os.path.isdir(aln_dir):
                        continue

                    aln_fa = join(aln_dir, each_ref_aln_file + ".aligned.fasta")
                    if not os.path.exists(aln_fa):
                        continue

                    rds2 = read_file.fasta(aln_fa)
                    for rows2 in rds2:
                        seq_id2 = rows2[0].strip()
                        if seq_id2 in missing_accs:
                            if seq_id2 in seqs:
                                write_file.write(
                                    seq_id2 + "\t" + each_ref_aln_file + "\t" + seqs[seq_id2] + "\n"
                                )
                                accessions[seq_id2] = 1

        print("Removing the Sequence redundancy")
        self.remove_redundancy_from_alignment(out_path)

    def remove_redundancy_from_alignment(self, file_path, delimiter="\t"):
        """
        Removes duplicate sequence_id (column 1). Keeps first occurrence. Preserves header.
        """
        seen = set()
        lines_to_write = []

        with open(file_path, "r") as infile:
            header = infile.readline()
            if header:
                lines_to_write.append(header)

            for line in infile:
                if not line.strip():
                    continue
                key = line.strip().split(delimiter)[0]
                if key not in seen:
                    seen.add(key)
                    lines_to_write.append(line)

        with open(file_path, "w") as outfile:
            outfile.writelines(lines_to_write)

        print(f"Redundancy removed. File updated: {file_path}")

    def create_insertion_table(self):
        """
        Writes:
          - Normal mode: <base_dir>/<output_dir>/insertions.tsv
          - Update mode : <base_dir>/update/<output_dir>/insertions.tsv
        """
        out_path = join(self.base_dir, self.effective_output_dir, "insertions.tsv")

        header = ["accession", "reference", "insertion"]
        with open(out_path, "w") as write_file:
            write_file.write("\t".join(header) + "\n")

            for aln_dir in os.listdir(self.nextalign_dir):
                aln_dir_path = join(self.nextalign_dir, aln_dir)
                if not os.path.isdir(aln_dir_path):
                    continue

                for each_aln_dir in os.listdir(aln_dir_path):
                    each_aln_path = join(aln_dir_path, each_aln_dir)
                    if not os.path.isdir(each_aln_path):
                        continue

                    ins_csv = join(each_aln_path, each_aln_dir + ".insertions.csv")
                    if not os.path.exists(ins_csv):
                        continue

                    with open(ins_csv) as f:
                        for each_line in islice(f, 1, None):  # skip header
                            parts = each_line.strip().split(",")
                            if len(parts) < 2:
                                continue
                            accession = parts[0].strip()
                            insertion = parts[1].strip()
                            if insertion:
                                data = [accession, each_aln_dir, insertion]
                                write_file.write("\t".join(data) + "\n")

    def process(self):
        out_dir_abs = join(self.base_dir, self.effective_output_dir)
        print(f"[mode] update={self.update} -> writing outputs to: {out_dir_abs}")

        blast_dictionary = self.load_blast_hits(self.blast_hits)
        self.created_alignment_table(blast_dictionary)
        self.create_insertion_table()


if __name__ == "__main__":
    parser = ArgumentParser(description="Generating tables sqlite DB")

    parser.add_argument(
        "-g",
        "--genbank_matrix",
        help="Genbank matrix table",
        default="tmp/GenBank-matrix/gB_matrix_raw.tsv",
    )
    parser.add_argument("-b", "--base_dir", help="base directory", default="tmp")
    parser.add_argument(
        "-o",
        "--output_dir",
        help='output directory name (default: "Tables"). '
             'In update mode outputs go to: <base_dir>/update/<output_dir>/',
        default="Tables",
    )
    parser.add_argument(
        "-bh",
        "--blast_hits",
        help="BLASTN unique hits",
        default="tmp/Blast/query_uniq_tophits.tsv",
    )
    parser.add_argument(
        "-p",
        "--paded_aln",
        help="Paded alignment file",
        default="tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta",
    )
    parser.add_argument(
        "-n",
        "--nextalign_dir",
        help="Nextalign aligned directory",
        default="tmp/Nextalign/",
    )
    parser.add_argument("-e", "--email", help="Email id", default="your-email@example.com")

    parser.add_argument(
        "--update",
        action="store_true",
        help="When set, write outputs into <base_dir>/update/<output_dir>/ (e.g. tmp/update/Tables/).",
    )

    args = parser.parse_args()

    processor = GenerateTables(
        args.genbank_matrix,
        args.base_dir,
        args.output_dir,
        args.blast_hits,
        args.paded_aln,
        args.nextalign_dir,
        args.email,
        update=args.update,
    )
    processor.process()
