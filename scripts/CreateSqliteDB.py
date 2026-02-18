import os
import re
import csv
import time
import sqlite3
import sys
import read_file
import subprocess
import pandas as pd
from Bio import SeqIO
from os.path import join
from datetime import datetime
from argparse import ArgumentParser
from collections import defaultdict

class CreateSqliteDB:
	def __init__(self, meta_data, features, pad_aln, gene_info, m49_countries, m49_interm_region, m49_regions, m49_sub_regions, proj_settings, fasta_sequence_file, insertions, host_taxa_file, base_dir, output_dir, db_name, db_status, tree_file=None, iqtree_file=None, usher_tree=None, cluster_tsv=None, cluster_min_seq_id=None, filtered_ids_file=None, filtered_details_file=None, tree_manifest=None):
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
		self.filtered_details_file = filtered_details_file
		self.tree_manifest = tree_manifest

	@staticmethod
	def _read_tree_file(tree_path):
		if not tree_path:
			return None
		try:
			with open(tree_path, "r", encoding="utf-8") as handle:
				return handle.read().strip()
		except FileNotFoundError:
			return None

	@staticmethod
	def _load_tree_manifest(manifest_path):
		if not manifest_path or not os.path.isfile(manifest_path):
			return []
		rows = []
		with open(manifest_path, "r", encoding="utf-8") as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				path_val = (row.get("path") or "").strip()
				if not path_val:
					continue
				rows.append({
					"source": (row.get("source") or "").strip() or "unknown",
					"name": (row.get("name") or "").strip() or None,
					"segment_key": (row.get("segment_key") or "").strip() or None,
					"path": path_val,
				})
		return rows

	@staticmethod
	def _segment_from_key(segment_key):
		if not segment_key:
			return None
		key = str(segment_key).strip()
		if not key:
			return None
		if key.isdigit():
			return key

		patterns = [
			r"(?:^|[_-])segment[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])seg[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])refset[_-]?(\d+)(?:$|[_-])",
			r"(?:^|[_-])(\d+)(?:$|[_-])",
		]
		for pat in patterns:
			m = re.search(pat, key, flags=re.IGNORECASE)
			if m:
				return m.group(1)
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

	@staticmethod
	def _require_file(path: str, label: str):
		if not path or not os.path.isfile(path):
			raise FileNotFoundError(f"{label} file not found: {path}")

	@staticmethod
	def _read_tsv_required(path: str, required_columns, label: str, dtype=None):
		df = pd.read_csv(path, sep="\t", dtype=dtype)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	@staticmethod
	def _read_csv_required(path: str, required_columns, label: str, dtype=None):
		df = pd.read_csv(path, sep=",", dtype=dtype)
		missing = [c for c in required_columns if c not in df.columns]
		if missing:
			raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")
		return df

	def _load_filtered_details(self) -> dict:
		"""Load filtered sequence reasons from filtered_sequences.tsv if available."""
		reasons = {}
		if not self.filtered_details_file or not os.path.isfile(self.filtered_details_file):
			return reasons
		try:
			df = pd.read_csv(self.filtered_details_file, sep="\t", dtype=str).fillna("")
		except Exception:
			return reasons

		if "seq_name" not in df.columns:
			return reasons

		for _, row in df.iterrows():
			seq = str(row.get("seq_name", "")).strip()
			if not seq:
				continue
			err = str(row.get("error", "")).strip()
			ref = str(row.get("reference", "")).strip()
			if err:
				reason = f"alignment_filtering: {err}"
			elif ref:
				reason = f"alignment_filtering: filtered against reference {ref}"
			else:
				reason = "alignment_filtering"
			reasons[seq] = reason
		return reasons

	def _add_cluster_column(self, df_meta_data):
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

		self._require_file(self.meta_data, "meta_data")
		self._require_file(self.features, "features")
		self._require_file(self.pad_aln, "pad_aln")
		self._require_file(self.gene_info, "gene_info")
		self._require_file(self.m49_countries, "m49_countries")
		self._require_file(self.m49_interm_region, "m49_interm_region")
		self._require_file(self.m49_regions, "m49_regions")
		self._require_file(self.m49_sub_regions, "m49_sub_regions")
		self._require_file(self.proj_settings, "proj_settings")
		self._require_file(self.fasta_sequence_file, "fasta_sequences")
		self._require_file(self.insertions, "insertions")
		self._require_file(self.host_taxa_file, "host_taxa_file")
		
		excluded_records = []

		# Load filtered sequence IDs to exclude
		filtered_ids = self._load_filtered_ids()
		filtered_details = self._load_filtered_details()
		if filtered_ids:
			print(f"[CreateSqliteDB] Excluding {len(filtered_ids)} filtered sequences from DB")
			for fid in filtered_ids:
				excluded_records.append({"primary_accession": fid, "reason": filtered_details.get(fid, "alignment_filtering")})
		
		df_meta_data = self._read_tsv_required(
			join(self.meta_data),
			["primary_accession"],
			"meta_data",
			dtype=str,
		)

		# Track reference/master rows so they are always retained in meta_data
		acc_type_col = "accession_type" if "accession_type" in df_meta_data.columns else None
		if acc_type_col:
			acc_type_norm = df_meta_data[acc_type_col].fillna("").str.strip().str.lower()
			is_ref_or_master = acc_type_norm.isin(["reference", "master"])
		else:
			is_ref_or_master = pd.Series(False, index=df_meta_data.index)
		
		# Exclude filtered sequences from meta_data, but never remove reference/master rows
		if filtered_ids and "primary_accession" in df_meta_data.columns:
			before_count = len(df_meta_data)
			remove_mask = df_meta_data["primary_accession"].isin(filtered_ids) & (~is_ref_or_master)
			df_meta_data = df_meta_data[~remove_mask]
			after_count = len(df_meta_data)
			if before_count != after_count:
				print(f"[CreateSqliteDB] Removed {before_count - after_count} filtered non-reference sequences from meta_data")

		is_ref_or_master = is_ref_or_master.reindex(df_meta_data.index, fill_value=False)
		
		# Collect exclusions from meta_data (e.g. invalid division)
		if "exclusion" in df_meta_data.columns:
			# Find rows with non-empty exclusion
			exclusion_mask = df_meta_data["exclusion"].notna() & (df_meta_data["exclusion"] != "")
			excluded_rows = df_meta_data[exclusion_mask]
			
			if not excluded_rows.empty:
				print(f"[CreateSqliteDB] Found {len(excluded_rows)} rows with exclusions in meta_data")
				for _, row in excluded_rows.iterrows():
					acc = row.get("primary_accession", "")
					if acc:
						excluded_records.append({"primary_accession": acc, "reason": row["exclusion"]})

				# Remove excluded rows from main meta_data, but never remove reference/master rows
				remove_mask = exclusion_mask & (~is_ref_or_master)
				df_meta_data = df_meta_data[~remove_mask]
				retained_refs = (exclusion_mask & is_ref_or_master).sum()
				if retained_refs:
					print(f"[CreateSqliteDB] Retained {retained_refs} reference/master rows despite exclusion flags")

		df_meta_data = self._add_cluster_column(df_meta_data)
		df_features = self._read_tsv_required(join(self.features), [], "features")
		df_aln = self._read_tsv_required(join(self.pad_aln), ["primary_accession"], "pad_aln")
		df_gene = self._read_tsv_required(join(self.gene_info), [], "gene_info")
		df_m49_country = self._read_csv_required(join(self.m49_countries), ["m49_code"], "m49_countries", dtype={'m49_code': str})
		df_m49_interm = self._read_csv_required(join(self.m49_interm_region), [], "m49_interm_region")
		df_m49_region = self._read_csv_required(join(self.m49_regions), [], "m49_regions")
		df_m49_sub_region = self._read_csv_required(join(self.m49_sub_regions), [], "m49_sub_regions")
		df_proj_setting = self._read_tsv_required(join(self.proj_settings), [], "proj_settings")
		df_insertions = self._read_tsv_required(join(self.insertions), ["primary_accession"], "insertions")
		df_host_taxa = self._read_tsv_required(join(self.host_taxa_file), ["primary_accession"], "host_taxa_file", dtype=str)
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
		
		if excluded_records:
			df_excluded = pd.DataFrame(excluded_records)
			df_excluded = df_excluded.drop_duplicates(subset=["primary_accession"])
			df_excluded.to_sql("excluded_accessions", conn, if_exists="replace", index=False)
			print(f"[CreateSqliteDB] Created excluded_accessions table with {len(df_excluded)} records")
		else:
			cursor.execute("CREATE TABLE IF NOT EXISTS excluded_accessions (primary_accession TEXT, reason TEXT)")

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
		cursor.execute("""CREATE TABLE IF NOT EXISTS excluded_accessions AS SELECT * FROM excluded_accessions;""")

		cursor.execute("""CREATE TABLE IF NOT EXISTS trees (name TEXT, source TEXT, segment_key TEXT, segment TEXT, newick TEXT, created_at TEXT);""")
		now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

		accession_to_segment = {}
		if "primary_accession" in df_meta_data.columns and "segment" in df_meta_data.columns:
			for _, row in df_meta_data[["primary_accession", "segment"]].dropna().iterrows():
				acc = str(row["primary_accession"]).strip()
				seg = str(row["segment"]).strip()
				if acc and seg and acc not in accession_to_segment:
					accession_to_segment[acc] = seg

		tree_records = []
		for name, source, tree_path in [
			("veryfasttree", "veryfasttree", self.tree_file),
			("iqtree", "iqtree", self.iqtree_file),
			("usher", "usher", self.usher_tree),
		]:
			newick = self._read_tree_file(tree_path)
			if newick:
				tree_records.append({
					"name": name,
					"source": source,
					"segment_key": None,
					"segment": None,
					"newick": newick,
				})

		for entry in self._load_tree_manifest(self.tree_manifest):
			newick = self._read_tree_file(entry["path"])
			if not newick:
				continue
			seg_key = entry.get("segment_key")
			seg_num = accession_to_segment.get(seg_key) if seg_key else None
			if seg_num is None:
				seg_num = self._segment_from_key(seg_key)
			name = entry.get("name") or f"{entry['source']}_{seg_key or 'tree'}"
			tree_records.append({
				"name": name,
				"source": entry["source"],
				"segment_key": seg_key,
				"segment": seg_num,
				"newick": newick,
			})

		seen_tree_names = set()
		for tr in tree_records:
			key = (tr["source"], tr["name"])
			if key in seen_tree_names:
				continue
			seen_tree_names.add(key)
			cursor.execute(
				"INSERT INTO trees (name, source, segment_key, segment, newick, created_at) VALUES (?, ?, ?, ?, ?, ?)",
				(tr["name"], tr["source"], tr["segment_key"], tr["segment"], tr["newick"], now_str),
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
			args.filtered_details,
			args.tree_manifest,
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
	parser.add_argument('--tree_manifest', help='TSV manifest with columns: source, name, segment_key, path', default=None)
	parser.add_argument('-ct', '--cluster_tsv', help='MMseqs clustering TSV (rep\tmember)', default=None)
	parser.add_argument('-ci', '--cluster_min_seq_id', help='MMseqs min sequence identity used for clustering', default=None)
	parser.add_argument('-fi', '--filtered_ids', help='File with filtered sequence IDs (one per line) to exclude from DB', default=None)
	parser.add_argument('-fd', '--filtered_details', help='TSV with filtered sequence details (seq_name, reference, error, warnings)', default=None)

	args = parser.parse_args()
	try:
		process(args)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		sys.exit(2)

