import os
import re
import csv
import time
import sqlite3
import read_file
import subprocess
import pandas as pd
from Bio import SeqIO
from os.path import join
from datetime import datetime
from argparse import ArgumentParser
from collections import defaultdict

class CreateSqliteDB:
	def __init__(self, meta_data, features, pad_aln, gene_info, m49_countries, m49_interm_region, m49_regions, m49_sub_regions, proj_settings, fasta_sequence_file, insertions, host_taxa_file, base_dir, output_dir, db_name, db_status, tree_file=None, iqtree_file=None, usher_tree=None, cluster_tsv=None, cluster_min_seq_id=None, filtered_ids_file=None):
		self.meta_data = meta_data
		self.features = features
		self.pad_aln = pad_aln
		self.gene_info = gene_info
		self.m49_countries = m49_countries
		self.m49_interm_region = m49_interm_region
		self.m49_regions = m49_regions
		self.m49_sub_regions = m49_sub_regions
		self.proj_settings = proj_settings
		self.fasta_sequence_file = fasta_sequence_file
		self.insertions = insertions
		self.host_taxa_file = host_taxa_file
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.db_name = db_name
		self.db_status = db_status
		self.tree_file = tree_file
		self.iqtree_file = iqtree_file
		self.usher_tree = usher_tree
		self.cluster_tsv = cluster_tsv
		self.cluster_min_seq_id = cluster_min_seq_id
		self.filtered_ids_file = filtered_ids_file

	@staticmethod
	def _read_tree_file(tree_path: str):
		if not tree_path:
			return None
		try:
			with open(tree_path, "r", encoding="utf-8") as handle:
				return handle.read().strip()
		except FileNotFoundError:
			return None

	def _load_filtered_ids(self) -> set:
		"""Load the set of sequence IDs that were filtered during alignment."""
		if not self.filtered_ids_file:
			return set()
		try:
			with open(self.filtered_ids_file, "r", encoding="utf-8") as f:
				return {line.strip() for line in f if line.strip()}
		except FileNotFoundError:
			return set()

	def _add_cluster_column(self, df_meta_data: pd.DataFrame):
		if not self.cluster_tsv:
			return df_meta_data
		try:
			cluster_df = pd.read_csv(self.cluster_tsv, sep="\t", header=None, dtype=str)
		except FileNotFoundError:
			return df_meta_data

		if cluster_df.shape[1] < 2:
			return df_meta_data

		cluster_df = cluster_df.iloc[:, :2]
		cluster_df.columns = ["cluster_rep", "member"]
		cluster_map = dict(zip(cluster_df["member"], cluster_df["cluster_rep"]))

		try:
			min_id = float(self.cluster_min_seq_id) if self.cluster_min_seq_id is not None else None
		except (TypeError, ValueError):
			min_id = None

		if min_id is not None:
			pct = int(round(min_id * 100))
			col_name = f"cluster_{pct}pct"
		else:
			col_name = "cluster"

		if "primary_accession" in df_meta_data.columns:
			df_meta_data[col_name] = df_meta_data["primary_accession"].map(cluster_map)
		return df_meta_data

	def load_fasta(self):
		fasta_data = []
		for record in SeqIO.parse(self.fasta_sequence_file, "fasta"):
			fasta_data.append({"header": record.id, "sequence": str(record.seq)})

		return pd.DataFrame(fasta_data)

	@staticmethod
	def _normalize_db_status(db_status: str):
		s = (db_status or "").strip().lower()
		if s in {"new", "new db", "create", "created", "fresh"}:
			return "new db"
		if s in {"modified", "update", "updated", "changed"}:
			return "last updated"
		return db_status


	def create_db(self):
		output_dir = join(self.base_dir, self.output_dir)
		os.makedirs(output_dir, exist_ok=True)
		
		# Load filtered sequence IDs to exclude
		filtered_ids = self._load_filtered_ids()
		if filtered_ids:
			print(f"[CreateSqliteDB] Excluding {len(filtered_ids)} filtered sequences from DB")
		
		df_meta_data = pd.read_csv(join(self.meta_data), sep="\t", dtype=str)
		
		# Exclude filtered sequences from meta_data
		if filtered_ids and "primary_accession" in df_meta_data.columns:
			before_count = len(df_meta_data)
			df_meta_data = df_meta_data[~df_meta_data["primary_accession"].isin(filtered_ids)]
			after_count = len(df_meta_data)
			if before_count != after_count:
				print(f"[CreateSqliteDB] Removed {before_count - after_count} filtered sequences from meta_data")
		
		df_meta_data = self._add_cluster_column(df_meta_data)
		df_features = pd.read_csv(join(self.features), sep="\t")
		df_aln = pd.read_csv(join(self.pad_aln), sep="\t")
		df_gene = pd.read_csv(join(self.gene_info), sep="\t")
		df_m49_country = pd.read_csv(join(self.m49_countries), dtype={'m49_code': str}, sep=",")
		df_m49_interm = pd.read_csv(join(self.m49_interm_region), sep=",")
		df_m49_region = pd.read_csv(join(self.m49_regions), sep=",")
		df_m49_sub_region = pd.read_csv(join(self.m49_sub_regions), sep=",")
		df_proj_setting = pd.read_csv(join(self.proj_settings), sep="\t")
		df_insertions = pd.read_csv(join(self.insertions), sep="\t")
		df_host_taxa = pd.read_csv(join(self.host_taxa_file), sep="\t", dtype=str)
		df_fasta_sequences = self.load_fasta()
		conn = sqlite3.connect(join(output_dir, self.db_name + ".db"))
		cursor = conn.cursor()

		df_meta_data.to_sql("meta_data", conn, if_exists="replace", index=False)
		df_features.to_sql("features", conn, if_exists="replace", index=False)
		df_aln.to_sql("sequence_alignment", conn, if_exists="replace", index=False)
		df_gene.to_sql("genes", conn, if_exists="replace", index=False)
		df_m49_country.to_sql("m49_country", conn, if_exists="replace", index=False)
		df_m49_interm.to_sql("m49_intermediate", conn, if_exists="replace", index=False)
		df_m49_region.to_sql("m49_regions", conn, if_exists="replace", index=False)
		df_m49_sub_region.to_sql("m49_sub_regions", conn, if_exists="replace", index=False)
		df_proj_setting.to_sql("project_settings", conn, if_exists="replace", index=False)
		df_fasta_sequences.to_sql("sequences", conn, if_exists="replace", index=False)
		df_insertions.to_sql("insertions", conn, if_exists="replace", index=False)
		df_host_taxa.to_sql("host_taxa", conn, if_exists="replace", index=False)
		cursor.execute("PRAGMA foreign_keys = ON;")

		cursor.execute("""CREATE TABLE IF NOT EXISTS meta_data AS SELECT * FROM meta_data;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS features AS SELECT * FROM features;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS sequence_alignment AS SELECT * FROM sequence_alignment;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS genes AS SELECT * FROM genes;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_country AS SELECT * FROM m49_country;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_intermediate AS SELECT * FROM m49_intermediate;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_regions AS SELECT * FROM m49_regions;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS m49_sub_regions AS SELECT * FROM m49_sub_regions;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS project_settings AS SELECT * FROM project_settings;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS sequences AS SELECT * FROM sequences;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS insertions AS SELECT * FROM insertions;""")
		cursor.execute("""CREATE TABLE IF NOT EXISTS host_taxa AS SELECT * FROM host_taxa;""")

		cursor.execute("""CREATE TABLE IF NOT EXISTS trees (name TEXT, source TEXT, newick TEXT, created_at TEXT);""")
		now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
		for name, source, tree_path in [
			("veryfasttree", "veryfasttree", self.tree_file),
			("iqtree", "iqtree", self.iqtree_file),
			("usher", "usher", self.usher_tree),
		]:
			newick = self._read_tree_file(tree_path)
			if newick:
				cursor.execute(
					"INSERT INTO trees (name, source, newick, created_at) VALUES (?, ?, ?, ?)",
					(name, source, newick, now_str),
				)

		cursor.execute("""CREATE TABLE IF NOT EXISTS info (creation_type TEXT,date TEXT);""")
		creation_type = self._normalize_db_status(self.db_status)

		info_df = pd.DataFrame([{"creation_type": creation_type, "date": now_str}])
		info_df.to_sql("info", conn, if_exists="append", index=False)

		conn.commit()
		conn.close()

def process(args):
	db_creator = CreateSqliteDB(
			args.meta_data,
			args.features,
			args.pad_aln,
			args.gene_info,
			args.m49_countries,
			args.m49_interm_region,
			args.m49_regions,
			args.m49_sub_regions,
			args.proj_settings,
			args.fasta_sequences,
			args.insertion_file,
			args.host_taxa_file,
			args.base_dir,
			args.output_dir,
			args.db_name,
			args.db_status,
			args.tree_file,
			args.iqtree_file,
			args.usher_tree,
			args.cluster_tsv,
			args.cluster_min_seq_id,
			args.filtered_ids,
		)
	db_creator.create_db()


if __name__ == "__main__":
	parser = ArgumentParser(description='Creating sqlite DB')
	parser.add_argument('-m', '--meta_data', help='Meta data table', default="tmp/GenBank-matrix/gB_matrix_raw.tsv")
	parser.add_argument('-b', '--base_dir', help='Base directory', default="tmp")
	parser.add_argument('-o', '--output_dir', help='tmp directory where the database is stored', default="SqliteDB")
	parser.add_argument('-rf', '--features', help='Features table', default="tmp/Tables/features.tsv")
	#parser.add_argument('-qf', '--query_features', help='Query feature table', default="tmp/Tables/query_features.tsv")
	parser.add_argument('-p', '--pad_aln', help='Padded alignment file', default="tmp/Tables/sequence_alignment.tsv")
	parser.add_argument('-g', '--gene_info', help='Gene table', default="generic/rabv/Tables/gene_info.csv")
	parser.add_argument('-mc', '--m49_countries', help='M49 countries', default="assets/m49_country.csv")
	parser.add_argument('-mir', '--m49_interm_region', help='M49 intermediate regions', default="assets/m49_intermediate_region.csv")
	parser.add_argument('-mr', '--m49_regions', help='M49 regions', default="assets/m49_region.csv")
	parser.add_argument('-msr', '--m49_sub_regions', help='M49 sub-regions', default="assets/m49_sub_region.csv")
	parser.add_argument('-s', '--proj_settings', help='Project settings', default="tmp/Software_info/software_info.tsv")
	parser.add_argument('-fa', '--fasta_sequences', help='Fasta sequences', default="tmp/GenBank-matrix/sequences.fa")
	parser.add_argument('-i', '--insertion_file', help='Nextalign insertion file', default="tmp/Tables/insertions.tsv")
	parser.add_argument('-ht', '--host_taxa_file', help='Host Taxanomy file', default="tmp/HostTaxa/Host_taxa.tsv")
	parser.add_argument('-d', '--db_name', help='Name of the Sqlite database', default="gdb")
	parser.add_argument('-ds', '--db_status', help='Database status: "new db" (default) or "last modified"/"last updated". Determines info.creation_type.',default="new db")
	parser.add_argument('-t', '--tree_file', help='VeryFastTree Newick file', default=None)
	parser.add_argument('-it', '--iqtree_file', help='IQ-TREE Newick file', default=None)
	parser.add_argument('-ut', '--usher_tree', help='UShER output Newick file', default=None)
	parser.add_argument('-ct', '--cluster_tsv', help='MMseqs clustering TSV (rep\tmember)', default=None)
	parser.add_argument('-ci', '--cluster_min_seq_id', help='MMseqs min sequence identity used for clustering', default=None)
	parser.add_argument('-fi', '--filtered_ids', help='File with filtered sequence IDs (one per line) to exclude from DB', default=None)

	args = parser.parse_args()

	process(args)

