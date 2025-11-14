import os
import csv
import sys
import argparse
from os.path import join

class HostTaxaTable:
	def __init__(self, gb_matrix, output_dir, names_file, nodes_file, base_dir, output_file):
		self.gb_matrix = gb_matrix
		self.output_dir = output_dir
		self.names_file = names_file
		self.nodes_file = nodes_file
		self.base_dir = base_dir
		self.output_file = output_file

	def load_taxa_ids_from_tsv(self):
		column="host_taxa_id"
		unique=True
		tsv_file = self.gb_matrix
		taxa_ids = []
		with open(tsv_file, newline="", encoding="utf-8") as f:
			reader = csv.DictReader(f, delimiter="\t")
			if column not in reader.fieldnames:
				raise ValueError(f"Column '{column}' not found in TSV file")
			for row in reader:
				val = row[column].strip()
				if val and val.upper() != "NA":
					taxa_ids.append(int(val))
		if unique:
			return sorted(set(taxa_ids))
		return taxa_ids

	def load_names(self):
		names = {}
		with open(self.names_file, encoding="utf-8") as f:
			for line in f:
				parts = [p.strip() for p in line.split("|")]
				taxid, name_txt, _, class_name = parts[:4]
				if class_name == "scientific name":
					names[int(taxid)] = name_txt
		return names

	def load_nodes(self):
		children_map = {}
		with open(self.nodes_file, encoding="utf-8") as f:
			for line in f:
				parts = [p.strip() for p in line.split("|")]
				taxid, parent_taxid = int(parts[0]), int(parts[1])
				children_map.setdefault(parent_taxid, []).append(taxid)
		return children_map

	def write_children(self):
		os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)
		children_map = self.load_nodes()
		names = self.load_names()

		taxa_list = self.load_taxa_ids_from_tsv()
		op_file = join(self.base_dir, self.output_dir, self.output_file)
		with open(op_file, "w", encoding="utf-8") as out:
			out.write("taxa_id\tscientific_name\tchild_taxa_id\tchild_scientific_name\n")
			for each_taxa in taxa_list:
				sci_name = names.get(each_taxa, "Unknown")
				children = children_map.get(each_taxa, [])
				if not children:
					out.write(f"{each_taxa}\t{sci_name}\tNA\tNA\n")
				else:
					for child in children:
						out.write(f"{each_taxa}\t{sci_name}\t{child}\t{names.get(child,'Unknown')}\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract children taxa for a given NCBI taxonomy ID")
	parser.add_argument("-g", "--gb_matrix", default="tmp-0/GenBank-matrix/gB_matrix_raw.tsv", help="GenBank matrix file")
	parser.add_argument("-o", "--output_dir", default="HostTaxa", help="Output directory")
	parser.add_argument("-n", "--names", default="tmp/Taxa/names.dmp",help="Path to names.dmp file (default: names.dmp)")
	parser.add_argument("-s", "--nodes", default="tmp/Taxa/nodes.dmp", help="Path to nodes.dmp file (default: nodes.dmp)")
	parser.add_argument("-b", "--base_dir", default="tmp", help="Path to base directory")
	parser.add_argument("-f", "--output_file", default="Host_taxa.tsv", help="Output TSV file")
	args = parser.parse_args()

	write_taxa = HostTaxaTable(args.gb_matrix, args.output_dir, args.names, args.nodes, args.base_dir, args.output_file)
	write_taxa.write_children() 

	#self.write_children(args.gb_matrix, args.output_dir, args.names, args.nodes, args.base_dir, args.output_file)
	print(f"Results written to {args.output_file}")
