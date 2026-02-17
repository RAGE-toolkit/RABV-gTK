#!/usr/bin/env python3
import os
import csv
import shutil
import tempfile
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser

'''
Running the script with update:
    python FilterAndExtractSequences.py -g tmp/update/GenBank-matrix/gB_matrix_raw.tsv -r ./../generic/rabv/ref_list.txt -b tmp -o Sequences --update

Runnning the script without update:
	python FilterAndExtractSequences.py -g tmp/GenBank-matrix/genbank_matrix.tsv -r tmp/ref_list.txt
'''

class FilterAndExtractSequences:
    def __init__(
        self,
        genbank_matrix,
        sequence_file,
        genbank_matrix_filtered,
        ref_file,
        base_dir,
        output_dir,
        total_length,
        real_length,
        prop_ambigious,
        segmented_virus,
        gb_division,
        valid_divisions,
        seq_type,
        update_mode=False,   # NEW
    ):
        self.genbank_matrix = genbank_matrix
        self.sequence_file = sequence_file
        self.genbank_matrix_filtered = genbank_matrix_filtered
        self.ref_file = ref_file
        self.base_dir = base_dir
        self.output_dir = output_dir
        self.total_length = total_length
        self.real_length = real_length
        self.prop_ambigious = prop_ambigious
        self.segmented_virus = segmented_virus
        self.gb_division = gb_division
        self.valid_divisions = valid_divisions
        self.seq_type = seq_type
        self.update_mode = update_mode  # NEW

        # NEW: route output under tmp/update when update mode is on
        if self.update_mode:
            self.out_root = join(self.base_dir, "update", self.output_dir)
        else:
            self.out_root = join(self.base_dir, self.output_dir)

        os.makedirs(self.out_root, exist_ok=True)

    def fasta_to_dict(self):
        seq_dict = {}
        for record in SeqIO.parse(self.sequence_file, "fasta"):
            seq_dict[record.id] = str(record.seq)
        return seq_dict

    def read_ref_list(self):
        ref_list = {}
        with open(self.ref_file) as f:
            for each_ref in f:
                split_col = each_ref.strip().split("\t")
                if len(split_col) < 2:
                    continue
                accs = split_col[0].strip()
                acc_type = split_col[1].strip()
                ref_list[accs] = acc_type
        return ref_list

    def read_ref_list_segmented_virus(self):
        ref_list = []
        with open(self.ref_file) as f:
            for each_ref_line in f:
                split_ref_line = each_ref_line.strip().split("|")
                if split_ref_line and split_ref_line[0].strip():
                    ref_list.append(split_ref_line[0].strip())
        return ref_list

    def filter_columns(self):
        sequence_dict = self.fasta_to_dict()
        exclusion_dict = {}

        if self.segmented_virus == "Y":
            ref_list = self.read_ref_list_segmented_virus()
        else:
            ref_list = self.read_ref_list()

        # NEW: write query/ref fasta into out_root
        query_seq_file = join(self.out_root, "query_seq.fa")
        ref_seq_file = join(self.out_root, "ref_seq.fa")

        with open(query_seq_file, "w") as write_query_seq, open(ref_seq_file, "w") as write_ref_seq:
            with open(self.genbank_matrix) as file, tempfile.NamedTemporaryFile("w", delete=False, newline="") as tmpfile:
                csv_reader = csv.DictReader(file, delimiter="\t")
                fieldnames = csv_reader.fieldnames or []

                if "exclusion_status" not in fieldnames:
                    fieldnames += ["exclusion_status"]
                if "exclusion_criteria" not in fieldnames:
                    fieldnames += ["exclusion_criteria"]

                writer = csv.DictWriter(tmpfile, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()

                for row in csv_reader:
                    if str(row.get("exclusion_status", "")).strip() == "1":
                        writer.writerow(row)
                        continue

                    gb_division = row.get("division", "")
                    seq_len = int(row.get("length", 0) or 0)
                    n_count = int(row.get("n", 0) or 0)
                    seq_len_without_n = seq_len - n_count

                    # NOTE: original code used gi_number; keep that behavior
                    accession = row.get("gi_number", "")
                    if accession not in sequence_dict:
                        # skip if sequence missing from fasta dict
                        row["exclusion_status"] = 1
                        row["exclusion_criteria"] = (row.get("exclusion_criteria", "").strip() + "; " if row.get("exclusion_criteria") else "") + "Sequence missing in FASTA"
                        writer.writerow(row)
                        continue

                    sequence = sequence_dict[accession]

                    exclusion_status = 0
                    exclusion_criteria = ""

                    if self.gb_division is not None and gb_division not in self.gb_division:
                        exclusion_status = 1
                        exclusion_criteria = f"GenBank division {gb_division} not in user-defined exclusion list"
                    elif gb_division not in self.valid_divisions:
                        exclusion_status = 1
                        exclusion_criteria = f"GenBank division {gb_division} not in valid division list"

                    if seq_len < self.total_length:
                        exclusion_status = 1
                        exclusion_criteria = "Sequence length less than total length"

                    if seq_len_without_n < self.real_length:
                        exclusion_status = 1
                        exclusion_criteria = "Real sequence length (length - Ns) less than the threshold"

                    existing_criteria = row.get("exclusion_criteria", "").strip()
                    if existing_criteria:
                        row["exclusion_criteria"] = f"{existing_criteria}; {exclusion_criteria}" if exclusion_criteria else existing_criteria
                    else:
                        row["exclusion_criteria"] = exclusion_criteria

                    row["exclusion_status"] = exclusion_status
                    writer.writerow(row)

                    if exclusion_status == 1:
                        exclusion_dict[accession] = row["exclusion_criteria"]
                        continue

                    if accession in ref_list:
                        write_ref_seq.write(f">{accession}\n{sequence}\n")
                    else:
                        write_query_seq.write(f">{accession}\n{sequence}\n")

        # Replace original matrix file with updated temp file
        # (keep original behavior: overwrites the input matrix)
        shutil.move(tmpfile.name, self.genbank_matrix)

        print(f"Wrote: {query_seq_file}")
        print(f"Wrote: {ref_seq_file}")
        print(f"Updated matrix in-place: {self.genbank_matrix}")

    def process(self):
        self.filter_columns()


if __name__ == "__main__":
    parser = ArgumentParser(description="Filter sequences and prepare sequences for BLAST alignment")
    parser.add_argument("-g", "--genbank_matrix", help="GenBank matrix file", required=True)
    parser.add_argument("-sf", "--sequence_file", help="FASTA sequence file", default="tmp/GenBank-matrix/sequences.fa")
    parser.add_argument("-f", "--genbank_matrix_filtered", help="Filtered GenBank matrix file storage directory", default="tmp/GenBank-matrix")
    parser.add_argument("-r", "--ref_file", help="Text file containing list of reference sequence accessions", required=True)
    parser.add_argument("-b", "--base_dir", help="Base directory", default="tmp")
    parser.add_argument("-o", "--output_dir", help="Output directory to process the data", default="Sequences")
    parser.add_argument("-l", "--total_length", help="Total length of the sequence should be >= length provided", default=1, type=int)
    parser.add_argument("-n", "--real_length", help="Length of sequences without N", default=1, type=int)
    parser.add_argument("-a", "--prop_ambigious_data", help="Proportion of ambiguous data to be excluded.", nargs="+", default=None, type=str)
    parser.add_argument("-v", "--segmented_virus", help="Is segmented virus", default="N")
    parser.add_argument(
        "-d",
        "--genbank_division",
        help=(
            "GenBank division to be excluded. General divisions are:\n"
            '    "VRL" = Viral sequences\n'
            '    "ENV" = Environmental sequences\n'
            '    "PAT" = Patented sequences\n'
            '    "SYN" = Synthetic and chimeric sequences\n'
        ),
        nargs="+",
        default=["VRL", "PAT", "SYN", "ENV"],
        type=str,
    )
    parser.add_argument("-vd", "--valid_divisions", help="Valid GenBank divisions to be considered", nargs="+", default=["VRL", "ENV"], type=str)
    parser.add_argument("-s", "--seq_type", help="Sequence type", default=None, type=str)

    # NEW: update flag
    parser.add_argument(
        "--update",
        action="store_true",
        help="If set, write outputs (query_seq.fa/ref_seq.fa) under tmp/update/<output_dir> instead of tmp/<output_dir>.",
    )

    args = parser.parse_args()

    processor = FilterAndExtractSequences(
        genbank_matrix=args.genbank_matrix,
        sequence_file=args.sequence_file,
        genbank_matrix_filtered=args.genbank_matrix_filtered,
        ref_file=args.ref_file,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        total_length=args.total_length,
        real_length=args.real_length,
        prop_ambigious=args.prop_ambigious_data,
        segmented_virus=args.segmented_virus,
        gb_division=args.genbank_division,
        valid_divisions=args.valid_divisions,
        seq_type=args.seq_type,
        update_mode=args.update,  # NEW
    )
    processor.process()
