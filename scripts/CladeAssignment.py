'''
Clade_assignments
	- read the meta data table and separate the accession based on reference and query
	- read the clade files (accession, major, minor clade)
	- read the padded alignment file and separate the sequences based on query and reference
	- perform $iqtree2 -s alUnc509RefseqsMafftHandModified.fa -m GTR+G -nt 7 -pre ref_tree #on reference sequences
	- perform $epa-ng on query_sequence_sets $epa-ng --redo -m GTR+G -t ref_tree.treefile -s alUnc509RefseqsMafftHandModified.fa -q query_seq_aln.fa  
	- perform $gappa to assign major clades gappa examine assign   --jplace-path epa_result.jplace   --taxon-file major_clade_and_accession.tsv   --out-dir gappa_assign_output   --per-query-results #to assign major clades

	- perform $gappa examine assign   --jplace-path epa_result.jplace   --taxon-file minor_clade_and_accession_noNull.tsv   --out-dir gappa_assign_output   --per-query-results #to assign minor clades
	- add the clade assignments to meta_data sheet
	
'''

import os
import sys
import csv
import shutil
import read_file
import subprocess
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser

class CladeAssignment:
	def __init__(self, major_clade, minor_clade, padded_alignment, base_dir, output_dir, gb_matrix, threads, iqtree_model):
		self.major_clade = major_clade
		self.minor_clade = minor_clade
		self.padded_alignment = padded_alignment
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.gb_matrix = gb_matrix
		self.threads = threads
		self.iqtree_model = iqtree_model
	
	def sequence_accessions(self):
		accessions = {}
		exc_type = {}
		with open(self.gb_matrix, newline='', encoding='utf-8') as file:
			reader = csv.DictReader(file, delimiter='\t') 
			for row in reader:
				
				if row["accession_type"] not in accessions:
					if row["exclusion_status"] != "1":
						accessions[row["accession_type"]] = [row["gi_number"]]
				else:
					if row["exclusion_status"] != "1":
						accessions[row["accession_type"]].append(row["gi_number"])

		return accessions

	def alignment(self):
		output_path = join(self.basce_dir, self.output_dir)
		os.makedirs(output_path, exist_ok=True)

		query_alignment = open(join(output_path, "query_aln.fa"), "w")
		reference_alignment = open(join(output_path, "reference_aln.fa"), "w")
  
		acc_type = self.sequence_accessions()
		alignment_dict = {}
		for rows in read_file.fasta(self.padded_alignment):
			if rows[0] in acc_type["query"]:
				query_alignment.write(">" + rows[0] + "\n" + rows[1].strip() + "\n")
			else:
				reference_alignment.write(">" + rows[0] + "\n" + rows[1].strip() + "\n")

	def iqtree(self, ref_fa, query_fa, prefix):
		command = [
			'iqtree2',
			'-s', ref_fa,
			'-m', self.iqtree_model,
			'-nt', self.threads,
			'-pre', prefix,
			]
		try:
			subprocess.run(command, check=True)
			print(f"iqtree ran successfully. Results saved in {prefix}")
		except subprocess.CalledProcessError as e:
			print(f"Error running iqtree: {e}")
		
	def epa_ng(self, ref_fa, query_fa, prefix):
		command = [
			'epa-ng',
			'--redo',
			'-m', self.iqtree_model,
			'-t', prefix + '.treefile',
			'-s', ref_fa,
			'-q', query_fa,
			'--outdir', join(self.base_dir, self.output_dir)
			]
		try:
			subprocess.run(command, check=True)
			print(f"EPA-ng ran successfully. ")
		except subprocess.CalledProcessError as e:
			print(f"Error running EPA-ng: {e}")

	def major_clades(self, prefix):
		command = [
			'gappa', 
			'examine', 
			'assign',
			'--jplace-path', join(prefix, 'epa_result.jplace'),
			'--taxon-file', self.major_clade,
			'--out-dir', join(prefix, 'gappa_major_clades'),
			'--per-query-results'
			]
		try:
			print(command)
			subprocess.run(command, check=True)
			print(f"Gappa major clades assigned successfully. ")
		except subprocess.CalledProcessError as e:
			print(f"Error running Gappa: {e}")

	def minor_clades(self, prefix):
		command = [
			'gappa', 
			'examine',
			'assign',
			'--jplace-path', join(prefix, 'epa_result.jplace'),
			'--taxon-file', self.minor_clade,
			'--out-dir', join(prefix, 'gappa_minor_clades'),
			'--per-query-results'
			]
		try:
			subprocess.run(command, check=True)
			print(f"Gappa minor clades assigned successfully.")
		except subprocess.CalledProcessError as e:
			print(f"Error running Gappa: {e}")
				
	def update_gb_matrix_with_clades(self, major_clades_file, minor_clades_file):
		major_clades = {}
		gb_matrix_file = self.gb_matrix

		with open(major_clades_file, 'r', newline='') as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				major_clades[row['name']] = row['taxopath']

		minor_clades = {}
		with open(minor_clades_file, 'r', newline='') as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				minor_clades[row['name']] = row['taxopath']

		temp_file = gb_matrix_file + '.tmp'
		with open(gb_matrix_file, 'r', newline='') as infile, \
			open(temp_file, 'w', newline='') as outfile:

				reader = csv.DictReader(infile, delimiter='\t')
				#fieldnames = reader.fieldnames + ['major_clade', 'minor_clade']
				base_fields = [f for f in reader.fieldnames if f not in ('major_clade', 'minor_clade')]
				new_fields = base_fields + ['major_clade', 'minor_clade']
				writer = csv.DictWriter(outfile, fieldnames=new_fields, delimiter='\t')
				writer.writeheader()

				for row in reader:
					row = {k: v for k, v in row.items() if k in base_fields}
					locus = row['gi_number']
					row['major_clade'] = major_clades.get(locus, '')
					row['minor_clade'] = minor_clades.get(locus, '')
					writer.writerow(row)

		os.replace(temp_file, gb_matrix_file)
		print(f"Updated file in place: {gb_matrix_file}")

	def process(self):
		ref_fa = join(self.base_dir, self.output_dir, 'reference_aln.fa')
		query_fa = join(self.base_dir, self.output_dir, 'query_aln.fa')
		prefix = join(self.base_dir, self.output_dir, 'ref_tree')
		output_dir = join(self.base_dir, self.output_dir)
		major_clade_output = join(self.base_dir, self.output_dir, 'gappa_major_clades', 'per_query.tsv')
		minor_clade_output = join(self.base_dir, self.output_dir, 'gappa_minor_clades', 'per_query.tsv')
		
		seperate_alignment = self.alignment()
		execute_iqtree = self.iqtree(ref_fa, query_fa, prefix)
		execute_epa_ng = self.epa_ng(ref_fa, query_fa, prefix)
		assign_major_clades = self.major_clades(output_dir)
		assign_minor_clades = self.minor_clades(output_dir)
		update_gb_matrix_with_clades = self.update_gb_matrix_with_clades(major_clade_output, minor_clade_output)	
	
if __name__ == "__main__":
	parser = ArgumentParser(description='Assignes the clades to the query sequences')
	parser.add_argument('-c', '--major_clade', help='Major clade file', default='generic/rabv/major_clades.tsv')
	parser.add_argument('-s', '--minor_clade', help='Minor clade file', default='generic/rabv/minor_clades.tsv')
	parser.add_argument('-p', '--padded_alignment', help='Padded alignment file', default='tmp/Pad-alignment/NC_001542.aligned_merged_MSA.fasta') 
	parser.add_argument('-b', '--base_dir', help='Base directory', default='tmp')
	parser.add_argument('-o', '--output_dir', help='Output directory', default='CladeAssignment')
	parser.add_argument('-g', '--gb_matrix', help='GenBank matrix file', default='tmp/GenBank-matrix/gB_matrix_raw.tsv')
	parser.add_argument('-t', '--threads', help='Threads to run iqtree', default='7')
	parser.add_argument('-m', '--iqtree_model', help='iqtree model to use', default='GTR+G')
	args = parser.parse_args()

	processor = CladeAssignment(
		args.major_clade,
		args.minor_clade,
		args.padded_alignment,
		args.base_dir,
		args.output_dir,
		args.gb_matrix,
		args.threads,
		args.iqtree_model
		)
	processor.process()
