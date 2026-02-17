#!/usr/bin/env python3
import os
import re
import csv
import read_file
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
from argparse import ArgumentParser
from TextFileHandler import TextFileLoader
from FastaHandler import RemoveRedundantSequence

'''
Running command with update flag:
	python NextalignAlignment.py --gb_matrix tmp/update/GenBank-matrix/gB_matrix_raw.tsv --query_dir tmp/update/Blast/grouped_fasta --ref_dir tmp/update/Blast/ref_seqs --ref_fa_file tmp/update/Sequences/ref_seq.fa --master_seq_dir tmp/update/Blast/master_seq --master_ref NC_001542 --update

Running command without update flag:
	python NextalignAlignment.py --gb_matrix tmp/update/GenBank-matrix/gB_matrix_raw.tsv --query_dir tmp/update/Blast/grouped_fasta --ref_dir tmp/update/Blast/ref_seqs --ref_fa_file tmp/update/Sequences/ref_seq.fa --master_seq_dir tmp/update/Blast/master_seq --master_ref NC_001542
'''

class NextalignAlignment:
	def __init__(
		self,
		gb_matrix,
		query_dir,
		ref_dir,
		ref_fa_file,
		master_seq_dir,
		tmp_dir,
		master_ref,
		nextalign_dir,
		reference_alignment,
		update_mode=False,          # NEW
		update_root="update",       # NEW
	):
		self.gb_matrix = gb_matrix
		self.query_dir = query_dir
		self.ref_dir = ref_dir
		self.ref_fa_file = ref_fa_file
		self.master_seq_dir = master_seq_dir
		self.master_ref = master_ref
		self.tmp_dir = tmp_dir
		self.min_seed = "44"
		self.seed_spacing = "50"
		self.min_match_rate = "0.1"
		self.reference_alignment = reference_alignment
		self.nextalign_dir = nextalign_dir

		# NEW
		self.update_mode = update_mode
		self.update_root = update_root

	@staticmethod
	def path_to_basename(file_path):
		path = os.path.basename(file_path)
		return path.split('.')[0]

	def _out_base(self):
		"""
		Base output directory:
		  normal: <tmp_dir>/<nextalign_dir>/
		  update : <tmp_dir>/update/<nextalign_dir>/
		"""
		if self.update_mode:
			return join(self.tmp_dir, self.update_root, self.nextalign_dir)
		return join(self.tmp_dir, self.nextalign_dir)

	def _ensure_outdirs(self, *dirs):
		for d in dirs:
			os.makedirs(d, exist_ok=True)

	def nextalign_master(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(ref_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]

		command_str = " ".join(command)
		print(f"Executing command: {command_str}")

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def nextalign_query(self, query_acc_path, ref_acc_path, query_aln_op):
		accession = self.path_to_basename(query_acc_path)
		command = [
			'nextalign', 'run',
			'--min-seeds', f'{self.min_seed}',
			'--seed-spacing', f'{self.seed_spacing}',
			'--min-match-rate', f'{self.min_match_rate}',
			'--input-ref', ref_acc_path,
			'--output-all', join(query_aln_op, f'{accession}'),
			'--output-basename', f'{accession}',
			'--include-reference',
			query_acc_path
		]

		command_str = " ".join(command)

		return_code = os.system(command_str)
		if return_code == 0:
			print(f"{accession} completed successfully.")
		else:
			print(f"{accession} failed with return code {return_code}")

	def update_gb_matrix(self, alignment_dir: list, gB_matrix_file):
		failed_accessions = {}

		for each_aln_type in alignment_dir:
			if not os.path.isdir(each_aln_type):
				continue
			for each_aln in os.listdir(each_aln_type):
				error_file_path = join(each_aln_type, each_aln, each_aln + ".errors.csv")
				if os.path.exists(error_file_path):
					with open(error_file_path) as f:
						for line in f:
							if "In sequence" in line:
								acc = line.split(",")[0].strip()
								error = line.split(",")[1].strip()
								failed_accessions.setdefault(acc, []).append(error)

		updated_rows = []
		with open(gB_matrix_file, newline='') as csvfile:
			reader = csv.DictReader(csvfile, delimiter='\t')
			fieldnames = reader.fieldnames or []

			if 'exclusion_status' not in fieldnames:
				fieldnames.append('exclusion_status')
			if 'exclusion_criteria' not in fieldnames:
				fieldnames.append('exclusion_criteria')

			for row in reader:
				gi = row.get('gi_number')
				if gi in failed_accessions:
					row['exclusion_status'] = '1'
					existing_criteria = row.get('exclusion_criteria', '')
					new_criteria = '; '.join(failed_accessions[gi])
					if existing_criteria:
						row['exclusion_criteria'] = existing_criteria + '; ' + new_criteria
					else:
						row['exclusion_criteria'] = new_criteria
				else:
					row['exclusion_status'] = row.get('exclusion_status', '0')
					row['exclusion_criteria'] = row.get('exclusion_criteria', '')
				updated_rows.append(row)

		with open(gB_matrix_file, 'w', newline='') as csvfile:
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
			writer.writeheader()
			writer.writerows(updated_rows)

	def process(self):
		# NEW: output under tmp/update/Nextalign/... when update_mode
		out_base = self._out_base()
		query_aln_output_dir = join(out_base, "query_aln")
		ref_aln_output_dir = join(out_base, "reference_aln")

		self._ensure_outdirs(out_base, query_aln_output_dir, ref_aln_output_dir)

		if self.reference_alignment:
			# align each query against its corresponding ref fasta
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file),
					query_aln_output_dir
				)
			self.update_gb_matrix([query_aln_output_dir], self.gb_matrix)

		else:
			# align query against reference sequence (per-file)
			for each_query_file in os.listdir(self.query_dir):
				ref_file = each_query_file
				self.nextalign_query(
					join(self.query_dir, each_query_file),
					join(self.ref_dir, ref_file),
					query_aln_output_dir
				)

			# align master/ref alignment
			for each_ref in os.listdir(self.master_seq_dir):
				self.nextalign_master(
					self.ref_fa_file,
					join(self.master_seq_dir, each_ref),
					ref_aln_output_dir
				)

			input_seq = join(ref_aln_output_dir, self.master_ref, self.master_ref + ".aligned.fasta")
			output_seq = input_seq
			unique_seqs = RemoveRedundantSequence(input_seq, output_seq)
			unique_seqs.remove_redundant_fasta()

			self.update_gb_matrix([query_aln_output_dir, ref_aln_output_dir], self.gb_matrix)


if __name__ == "__main__":
	parser = ArgumentParser(description='Performs the nextalign of each sequence')

	# NEW: update flag to route outputs to tmp/update/
	parser.add_argument('--update', action='store_true',
	                    help='If set, write nextalign outputs under tmp/update/<nextalign_dir>/...')

	parser.add_argument('-g', '--gB_matrix', help='GenBank matrix (meta data) file.',
	                    default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-q', '--query_dir', help='Query file directory.',
	                    default="tmp/Blast/grouped_fasta")
	parser.add_argument('-r', '--ref_dir', help='Reference fasta directory',
	                    default="tmp/Blast/ref_seqs")
	parser.add_argument('-f', '--ref_fa_file', help='Reference fasta sequences',
	                    default="tmp/Sequences/ref_seq.fa")
	parser.add_argument('-ms', '--master_seq_dir', help='Master sequence directory',
	                    default="tmp/Blast/master_seq")
	parser.add_argument('-t', '--tmp_dir', help='Temp directory to process the data',
	                    default="tmp")
	parser.add_argument('-m', '--master_ref', help='Master reference accession. e.g. NC_001542', required=True)
	parser.add_argument('-n', '--nextalign_dir', help='Nextalign output directory name',
	                    default="Nextalign")
	parser.add_argument('-ra', '--ref_alignment_file',
	                    help='Use your own reference alignment file instead of Nextalign perfoms the alignment of reference against the master reference sequence')

	args = parser.parse_args()

	processor = NextalignAlignment(
		gb_matrix=args.gB_matrix,
		query_dir=args.query_dir,
		ref_dir=args.ref_dir,
		ref_fa_file=args.ref_fa_file,
		master_seq_dir=args.master_seq_dir,
		tmp_dir=args.tmp_dir,
		master_ref=args.master_ref,
		nextalign_dir=args.nextalign_dir,
		reference_alignment=args.ref_alignment_file,
		update_mode=args.update,   # NEW
		update_root="update",      # NEW
	)
	processor.process()
