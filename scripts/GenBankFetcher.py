import os
import csv
import requests
import time
from time import sleep
from os.path import join
from argparse import ArgumentParser

class GenBankFetcher:
	def __init__(self, taxid, base_url, email, output_dir, batch_size, sleep_time, base_dir, update_file, test_run=False, ref_list=None):
		self.taxid = taxid
		self.base_url = base_url
		self.email = email
		self.output_dir = output_dir
		self.batch_size = batch_size
		self.efetch_batch_size = 100
		self.sleep_time = sleep_time
		self.base_dir = base_dir
		self.update_file = update_file
		self.test_run = test_run
		self.ref_list = ref_list

	def get_record_count(self):
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmode=json&email={self.email}"
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return int(data['esearchresult']['count'])

	def fetch_ids(self):
		retmax = self.get_record_count()
		#search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmax={retmax}&usehistory=y&email={self.email}&retmode=json"
		search_url = (
        f"{self.base_url}esearch.fcgi?db=nucleotide"
        f"&term=txid{self.taxid}[Organism:exp]"
        f"&retmax={retmax}&idtype=acc"
        f"&usehistory=y&email={self.email}&retmode=json"
    	)
		
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return data["esearchresult"]["idlist"]

	def fetch_accs(self):
		retmax = self.get_record_count()
		search_url = f"{self.base_url}esearch.fcgi?db=nucleotide&term=txid{self.taxid}[Organism:exp]&retmax={retmax}&idtype=acc&usehistory=y&email={self.email}&retmode=json"
		response = requests.get(search_url)
		response.raise_for_status()
		data = response.json()
		return data["esearchresult"]["idlist"]


	def fetch_accs(self):
		start_all = time.time()

		# 1) inicializa la historia y recoge WebEnv, QueryKey y count
		t0 = time.time()
		hist_url = (
			f"{self.base_url}esearch.fcgi?db=nucleotide"
			f"&term=txid{self.taxid}[Organism:exp]"
			f"&retmax=0&idtype=acc"
			f"&usehistory=y&email={self.email}&retmode=json"
		)
		hist = requests.get(hist_url).json()["esearchresult"]
		webenv   = hist["webenv"]
		querykey = hist["querykey"]
		count    = int(hist["count"])
		print(f"[fetch_accs] ESearch history → count={count:,} took {time.time() - t0:.1f}s")

		# 2) paginación en bloques de self.batch_size
		all_accs = []
		for start in range(0, count, self.batch_size):
			t1 = time.time()
			page_url = (
				f"{self.base_url}esearch.fcgi?db=nucleotide"
				f"&WebEnv={webenv}&query_key={querykey}"
				f"&retstart={start}&retmax={self.batch_size}"
				f"&idtype=acc&retmode=json"
			)
			page = requests.get(page_url).json()["esearchresult"]["idlist"]
			page = [acc.split('.', 1)[0] for acc in page]
			all_accs.extend(page)
			print(f"[fetch_accs] chunk {start:,}-{min(start+self.batch_size, count):,} "
				f"({len(page):,} records) took {time.time() - t1:.1f}s")

		print(f"[fetch_accs] total accs fetched: {len(all_accs):,} in {time.time() - start_all:.1f}s")
		return all_accs


	def fetch_genbank_data(self, ids):
		# usa un batch interno distinto sólo para efetch
		batch_n = self.efetch_batch_size
		if self.test_run:
			ids=ids[:50]
			# read as tsv self.ref_list add to ids, refs are first column
			ref_list = []
			with open(self.ref_list, 'r') as f:
				for line in f:
					ref_list.append(line.strip().split('\t')[0])
			ids.extend(ref_list)
					
   

			batch_n=10
	
		max_ids=len(ids)
		for i in range(0, max_ids, batch_n):
			chunk = ids[i:i+batch_n]
			ids_str = ",".join(chunk)
			url = (
				f"{self.base_url}efetch.fcgi?db=nucleotide"
				f"&id={ids_str}"
				f"&retmode=xml&email={self.email}"
			)
			resp = requests.get(url)
			resp.raise_for_status()

			# aquí pasamos batch_n para nombrar el fichero igual que antes
			self.save_data(resp.text, i + batch_n)

			print(f"Downloaded XML {i+1:,}–{min(i+batch_n, len(ids)):,}")
			sleep(self.sleep_time)



	def save_data(self, data, batch_size):
		os.makedirs(self.output_dir, exist_ok=True)
		os.makedirs(join(self.output_dir, self.base_dir), exist_ok=True)

		base_filename = join(self.output_dir, self.base_dir, f"batch-{batch_size}")
		filename = f"{base_filename}.xml"
		counter = 1
		while os.path.exists(filename):
			filename = f"{base_filename}_{counter}.xml"
			counter += 1

		with open(filename, "w") as file:
				file.write(data)
		print(f"Data written to: {filename}")

	def download(self):
		ids = self.fetch_ids()
		print(f"Found {len(ids)} IDs")
		self.fetch_genbank_data(ids)

	def update(self, update_file):
		with open(update_file, 'r') as file:
			reader = csv.DictReader(file, delimiter='\t')
			if 'primary_accession' not in reader.fieldnames:
					raise ValueError("Expecting a column called 'primary_accession'")
			primary_accession = [row['primary_accession'] for row in reader]
			# check accession_version format
			print("first 10 primary_accession in TSV:", primary_accession[:10])

			print(f"Found {len(primary_accession)} primary_accession")
			ids = self.fetch_accs()
			
			# check NCBI IDs format
			print("first 10 IDs in NCBI:", ids[:10])
			primary_accession_set = set(primary_accession)
			missing_ids = [id for id in ids if id not in primary_accession_set]
			print(f"Found {len(missing_ids)} missing IDs")
			self.fetch_genbank_data(missing_ids)

if __name__ == "__main__":
	parser = ArgumentParser(description='This script downloads and updates GenBank XML files for a given taxonomic group identified by an NCBI taxonomic ID, usually a species')
	parser.add_argument('-t', '--taxid', help='TaxID example: 11292', required=True)
	parser.add_argument('-o', '--tmp_dir', help='Output directory where all the XML files are stored', default='tmp')
	parser.add_argument('-b', '--batch_size', help='Max number of XML files to pull and merge in a single file', default=100000, type=int)
	parser.add_argument('--update', help='Run the script in update mode to download only new sequences, expecting a TSV with a column called Accession Version')
	parser.add_argument('-u', '--base_url', help='Base URL to download the XML files', default='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/')
	parser.add_argument('-e', '--email', help='Email ID', default='your_email@example.com')
	parser.add_argument('-s', '--sleep_time', help='Delay after each set of information fetch', default=2, type=int)
	parser.add_argument('-d', '--base_dir', help='Directory where all the XML files are stored', default='GenBank-XML')
	parser.add_argument("--test_run", action="store_true", help="Run a test fetching only a few records for quick testing")
	parser.add_argument('--ref_list', help='Reference accession list for test run', default='generic/rabv/ref_list.txt')
	args = parser.parse_args()

	fetcher = GenBankFetcher(
		taxid=args.taxid,
		base_url=args.base_url,
		email=args.email,
		output_dir=args.tmp_dir,
		batch_size=args.batch_size,
		sleep_time=args.sleep_time,
		base_dir = args.base_dir,
		update_file = args.update,
		test_run = args.test_run,
		ref_list = args.ref_list
	)

	if args.update:
		fetcher.update(args.update)
	else:
		fetcher.download()

