import os
import csv
import shutil
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
import read_file

'''
	
	python BlastAlignment.py -b tmp -t Blast -q tmp/update/Sequences/query_seq.fa -m NC_001542 --update -r clade_files/ref_seq.fa -g tmp/update/GenBank-matrix/gB_matrix_raw.tsv


'''

class BlastAlignment:
	def __init__(
		self,
		query_fasta,
		db_fasta,
		base_dir,
		output_dir,
		output_file,
		is_segmented_virus,
		master_acc,
		is_update,
		keep_blast_tmp_dir,
		gb_matrix,
		segment_file=None,
		update_mode=False,          # NEW (directory routing)
		update_root_dir="update",   # NEW
	):
		self.query_fasta = query_fasta
		self.db_fasta = db_fasta
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.is_segmented_virus = is_segmented_virus
		self.segment_file = segment_file
		self.master_acc = master_acc
		self.gb_matrix = gb_matrix
		self.is_update = is_update
		self.keep_blast_tmp_dir = keep_blast_tmp_dir
		self.db_file_name = os.path.basename(db_fasta)

		# NEW: if update_mode, write everything under tmp/update/<output_dir>
		self.update_mode = update_mode
		self.update_root_dir = update_root_dir

		# Base output path used by the pipeline
		self.out_root = self._resolve_out_root()

	def _resolve_out_root(self):
		if self.update_mode:
			return join(self.base_dir, self.update_root_dir, self.output_dir)
		return join(self.base_dir, self.output_dir)

	def _ensure_out_root(self):
		os.makedirs(self.out_root, exist_ok=True)

	def read_gb_matrix(self):
		accessions = {}
		with open(self.gb_matrix, newline="", encoding="utf-8") as file:
			reader = csv.DictReader(file, delimiter="\t")
			for row in reader:
				# NOTE: gb_matrix used here assumes a 'sequence' column exists
				accessions[row["gi_number"]] = row["sequence"]
		return accessions

	def update(self, query_tmp_dir):
		if os.path.exists(self.query_fasta):
			query_accession = [record.id for record in SeqIO.parse(self.query_fasta, "fasta")]
			db_accession = [record.id for record in SeqIO.parse(self.db_fasta, "fasta")]
			all_accessions = self.read_gb_matrix()
			missing_accessions = [acc for acc in all_accessions.keys() if acc not in query_accession and acc not in db_accession]
			missing_accessions_count = len(missing_accessions)
			if missing_accessions_count > 0:
				print(f"{missing_accessions_count} new sequences to process")
				with open(join(query_tmp_dir, "query.fa"), "w") as write_file:
					for each_missing_acc in missing_accessions:
						write_file.write(">" + each_missing_acc + "\n")
						write_file.write(all_accessions[each_missing_acc] + "\n")
			else:
				print("No new accession id's to process")
		else:
			print("No query file exists for update, you need to run it without update to perform blast on existing sequences")

	@staticmethod
	def check_blast_exists(command):
		try:
			subprocess.run([command, "-version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			print(f"{command} is working.")
			return True
		except subprocess.CalledProcessError:
			print(f"{command} is not installed correctly.")
			return False
		except FileNotFoundError:
			print(f"{command} is not found on the system.")
			return False

	@staticmethod
	def delete_directory(dir_path):
		if os.path.exists(dir_path):
			if os.path.isdir(dir_path):
				shutil.rmtree(dir_path)
				print(f"Directory '{dir_path}' has been deleted. Ignore the message")
			else:
				print(f"'{dir_path}' exists but is not a directory. Ignore the message")
		else:
			print(f"Directory '{dir_path}' does not exist. Ignore the message")

	@staticmethod
	def ref_segments(query_tophits_annotated):
		segment_dict = {}
		for each_line in open(query_tophits_annotated):
			query, ref, score, strand, segment = each_line.strip().split("\t")
			if segment not in segment_dict:
				segment_dict[segment] = {}
			if ref not in segment_dict[segment]:
				segment_dict[segment][ref] = []
			segment_dict[segment][ref].append(query)
		return segment_dict

	def run_makeblastdb(self, tmp_dir):
		os.makedirs(join(tmp_dir, "DB"), exist_ok=True)
		command = [
			"makeblastdb",
			"-in", self.db_fasta,
			"-out", join(tmp_dir, "DB", self.db_file_name),
			"-title", "alignment",
			"-dbtype", "nucl",
		]
		try:
			subprocess.run(command, check=True)
			print(f"makeblastdb ran successfully on {self.db_fasta}")
		except subprocess.CalledProcessError as e:
			print(f"Error running makeblastdb: {e}")

	def run_blastn(self, output_dir, query_file):
		command = [
			"blastn",
			"-query", query_file,
			"-db", join(output_dir, "DB", self.db_file_name),
			"-task", "blastn",
			"-max_target_seqs", "1",
			"-max_hsps", "1",
			"-out", join(output_dir, self.output_file),
			"-outfmt", "6 qacc sacc pident sstrand",
		]
		try:
			subprocess.run(command, check=True)
			print(f"blastn ran successfully. Results saved in {self.output_file}")
		except subprocess.CalledProcessError as e:
			print(f"Error running blastn: {e}")

	def write_master_seq(self, output_dir):
		with open(self.db_fasta, "r") as infile:
			records = SeqIO.parse(infile, "fasta")
			selected_records = [record for record in records if record.id == self.master_acc]

		if selected_records:
			with open(join(output_dir, self.master_acc + ".fasta"), "w") as outfile:
				SeqIO.write(selected_records, outfile, "fasta")
				print(f"Sequence '{self.master_acc}' has been saved to {join(output_dir, self.master_acc)}")
		else:
			print(f"Sequence ID '{self.master_acc}' not found in {self.db_fasta}")

	def process_non_segmented_virus(self, output_dir, query_fasta):
		input_file = join(output_dir, "query_tophits.tsv")
		query_tophit_uniq = join(output_dir, "query_uniq_tophits.tsv")
		grouped_fasta = join(output_dir, "grouped_fasta")
		sorted_fasta = join(output_dir, "sorted_fasta")
		merged_fasta = join(output_dir, "merged_fasta")
		sorted_all = join(output_dir, "sorted_all")
		ref_seq_dir = join(output_dir, "ref_seqs")
		master_seq = join(output_dir, "master_seq")

		os.makedirs(grouped_fasta, exist_ok=True)
		os.makedirs(ref_seq_dir, exist_ok=True)
		os.makedirs(sorted_fasta, exist_ok=True)
		os.makedirs(merged_fasta, exist_ok=True)
		os.makedirs(sorted_all, exist_ok=True)
		os.makedirs(master_seq, exist_ok=True)

		records = {}
		values = {}

		self.write_master_seq(master_seq)

		with open(input_file, newline="") as file:
			reader = csv.reader(file, delimiter="\t")
			for row in reader:
				col1, col2, col3, col4 = row[0], row[1], float(row[2]), row[3]
				if col1 in values:
					existing_value = values[col1]
					if col3 > existing_value:
						records[col1] = [col1, col2, col3, col4]
						values[col1] = col3
				else:
					records[col1] = [col1, col2, col3, col4]
					values[col1] = col3

		with open(query_tophit_uniq, "w", newline="") as file:
			writer = csv.writer(file, delimiter="\t")
			for record in records.values():
				writer.writerow(record)

		seq_dicts = {}
		query_seqs = read_file.fasta(query_fasta)
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		# Separate plus and minus strand sequences
		for each_line in open(query_tophit_uniq, "r"):
			query_acc, ref_acc, identity, strand = each_line.strip().split("\t")
			if strand == "plus":
				with open(join(sorted_fasta, "plus.fa"), "a") as file_plus:
					file_plus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")
			else:
				with open(join(sorted_fasta, "minus.fa"), "a") as file_minus:
					file_minus.write(">" + query_acc + "\n" + seq_dicts[query_acc] + "\n")

		for each_file in os.listdir(sorted_fasta):
			if "minus" in each_file:
				command = ["seqkit", "seq", "-r", "-p", "-v", "-t", "dna", join(sorted_fasta, each_file), ">", join(merged_fasta, each_file)]
			else:
				command = ["cp", join(sorted_fasta, each_file), join(merged_fasta, each_file)]

			try:
				print(" ".join(command))
				os.system(" ".join(command))
				print(f"seqkit ran successfully for {each_file}")
			except subprocess.CalledProcessError as e:
				print(f"Error running seqkit: {e}")

		file_list = []
		for each_file in os.listdir(merged_fasta):
			prefix = merged_fasta + "/"
			file_list.append(prefix + each_file)

		command = ["cat", " ".join(file_list), ">", join(sorted_all, "query_seq.fa")]
		print(" ".join(command))
		try:
			os.system(" ".join(command))
			print("concatenation successful")
		except subprocess.CalledProcessError as e:
			print(f"Error in concatenation: {e}")

		grouped_dict = {}
		for each_line in open(query_tophit_uniq, "r"):
			query_acc, ref_acc, identity, strand = each_line.strip().split("\t")
			if ref_acc not in grouped_dict:
				grouped_dict[ref_acc] = [query_acc]
			else:
				grouped_dict[ref_acc].append(query_acc)

		#print("\n".join(grouped_dict.keys()))
		seq_dicts = {}
		query_seqs = read_file.fasta(join(sorted_all, "query_seq.fa"))
		for rows in query_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for each_ref_acc, list_of_query_acc in grouped_dict.items():
			with open(join(grouped_fasta, each_ref_acc + ".fasta"), "a") as write_file:
				for each_query_acc in list_of_query_acc:
					seqs = seq_dicts[each_query_acc]
					write_file.write(">" + each_query_acc + "\n")
					for i in range(0, len(seqs), 80):
						write_file.write(seqs[i:i + 80] + "\n")

		# writing reference sequences into individual fasta files
		ref_seqs = read_file.fasta(self.db_fasta)
		for rows in ref_seqs:
			seq_dicts[rows[0].strip()] = rows[1].strip()

		for ref_accs in grouped_dict.keys():
			seqs = seq_dicts[ref_accs]
			with open(join(ref_seq_dir, ref_accs + ".fasta"), "w") as write_file:
				write_file.write(">" + ref_accs + "\n")
				for i in range(0, len(seqs), 80):
					write_file.write(seqs[i:i + 80] + "\n")

	def update_gB_matrix(self, query_fasta, query_tophit_uniq, gB_matrix_file):
		uniq_blast_acc = {}
		with open(query_tophit_uniq) as f:
			for line in f:
				query_acc, ref_acc, score, strand = line.strip().split("\t")
				uniq_blast_acc[query_acc] = ref_acc

		query_acc_status = {}
		read_query_obj = read_file.fasta(query_fasta)
		for each_seq_obj in read_query_obj:
			header = each_seq_obj[0]
			if header not in uniq_blast_acc:
				query_acc_status[header] = 1  # excluded
			else:
				query_acc_status[header] = 0  # not excluded

		updated_rows = []
		with open(gB_matrix_file, newline="") as infile:
			reader = csv.DictReader(infile, delimiter="\t")
			fieldnames = reader.fieldnames or []

			if "exclusion_status" not in fieldnames:
				fieldnames.append("exclusion_status")
			if "exclusion_criteria" not in fieldnames:
				fieldnames.append("exclusion_criteria")

			for row in reader:
				gi = row.get("gi_number")
				status = query_acc_status.get(gi, 0)
				existing_criteria = row.get("exclusion_criteria", "")
				if status == 1 and gi in query_acc_status:
					new_criteria = "excluded due to no hit"
					if existing_criteria:
						row["exclusion_criteria"] = f"{existing_criteria}; {new_criteria}"
					else:
						row["exclusion_criteria"] = new_criteria
					row["exclusion_status"] = "1"
				else:
					row["exclusion_criteria"] = existing_criteria

				updated_rows.append(row)

		with open(gB_matrix_file, "w", newline="") as outfile:
			writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
			writer.writeheader()
			writer.writerows(updated_rows)

	def process(self):
		print(f"Using {self.gb_matrix} GenBank matrix file")
		self.check_blast_exists("blastn")

		# ensure base output root exists
		self._ensure_out_root()

		if self.is_segmented_virus == "Y" and not self.segment_file:
			raise ValueError("Missing segment file with accession and segment information")

		if self.is_segmented_virus == "Y":
			self.run_makeblastdb(self.out_root)
			self.run_blastn(self.out_root, self.query_fasta)

			for each_segment_dir in [
				join(self.out_root, "segment_sorted"),
				join(self.out_root, "segment_sorted_all"),
				join(self.out_root, "segment_merged_fasta"),
			]:
				self.delete_directory(each_segment_dir)

			raise NotImplementedError(
				"Segmented-virus path in the provided script appears incomplete (e.g., write_master_file). "
				"If you paste the complete working version, I’ll add update-mode routing safely."
			)

		else:
			if self.is_update == "Y":
				blast_tmp_dir = join(self.out_root, "tmp_dir")
				os.makedirs(blast_tmp_dir, exist_ok=True)

				self.update(blast_tmp_dir)

				qfa = join(blast_tmp_dir, "query.fa")
				if not os.path.exists(qfa) or os.path.getsize(qfa) == 0:
					print("[update] No new sequences written to query.fa; skipping BLAST.")
					return

				self.run_makeblastdb(blast_tmp_dir)
				self.run_blastn(blast_tmp_dir, qfa)
				self.process_non_segmented_virus(blast_tmp_dir, qfa)
				self.update_gB_matrix(qfa, join(blast_tmp_dir, "query_uniq_tophits.tsv"), self.gb_matrix)

				if self.keep_blast_tmp_dir == "N":
					shutil.rmtree(blast_tmp_dir)

			else:
				self.run_makeblastdb(self.out_root)
				self.run_blastn(self.out_root, self.query_fasta)
				self.process_non_segmented_virus(self.out_root, self.query_fasta)
				self.update_gB_matrix(self.query_fasta, join(self.out_root, "query_uniq_tophits.tsv"), self.gb_matrix)


if __name__ == "__main__":
	parser = ArgumentParser(description="Performs the BLAST alignment of query sequences against the given reference sequences")

	# NOTE: when --update_dir is set, defaults can be automatically redirected
	parser.add_argument("-b", "--base_dir", help="Base directory", default="tmp")
	parser.add_argument("-t", "--output_dir", help="Output directory", default="Blast")
	parser.add_argument("--update", action="store_true",
	                    help="If set, write outputs under tmp/update/<output_dir> instead of tmp/<output_dir>.")

	# These defaults stay as you had them (normal-mode). If --update_dir is set and user does not override,
	# we remap them below to tmp/update/Sequences/...
	parser.add_argument("-q", "--query_fa", help="Query fasta file", default="tmp/Sequences/query_seq.fa")
	parser.add_argument("-r", "--ref_fa", help="Blast DB fasta file", default="tmp/Sequences/ref_seq.fa")

	parser.add_argument("-o", "--output_file", help="Output file", default="query_tophits.tsv")
	parser.add_argument("-s", "--is_segmented_virus", help="Type Y for segmented virus else N", default="N")
	parser.add_argument("-f", "--segment_file", help="File containing information about the segments")
	parser.add_argument("-m", "--master_acc", help="Master accession (e.g. NC_001542)", required=True)
	parser.add_argument("-u", "--is_update", help="Only BLAST new sequences (Y/N)", default="N")
	parser.add_argument("-k", "--keep_blast_tmp_dir", help="Retain BLAST temp directory for debugging (Y/N)", default="N")
	parser.add_argument("-g", "--gb_matrix", help="GenBank matrix file", default="tmp/GenBank-matrix/gB_matrix_raw.tsv")

	args = parser.parse_args()

	# Review it and delete. If user enabled update_dir and didn't override query/ref defaults, redirect them.
	#if args.update:
	#	if args.query_fa == "tmp/Sequences/query_seq.fa":
	#		args.query_fa = join(args.base_dir, "update", "Sequences", "query_seq.fa")
	#	if args.ref_fa == "tmp/Sequences/ref_seq.fa":
	#		args.ref_fa = join(args.base_dir, "update", "Sequences", "ref_seq.fa")

	processor = BlastAlignment(
		query_fasta=args.query_fa,
		db_fasta=args.ref_fa,
		base_dir=args.base_dir,
		output_dir=args.output_dir,
		output_file=args.output_file,
		is_segmented_virus=args.is_segmented_virus,
		master_acc=args.master_acc,
		is_update=args.is_update,
		keep_blast_tmp_dir=args.keep_blast_tmp_dir,
		gb_matrix=args.gb_matrix,
		segment_file=args.segment_file,
		update_mode=args.update,   # NEW
		update_root_dir="update",      # NEW
	)
	processor.process()
