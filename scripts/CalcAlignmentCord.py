#!/usr/bin/env python3
import os
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from GffToDictionary import GffDictionary
from CalcGenomeCords import CalculateGenomeCoordinates

'''
	python CalcAlignmentCord.py -i tmp/update/Pad-alignment/ -m NC_001542 -bh tmp/update/Blast/query_uniq_tophits.tsv -g tmp/update/Gff/NC_001542.gff3 --update
'''

class CalculateAlignmentCoordinates:
	def __init__(
		self,
		paded_alignment,
		master_gff,
		tmp_dir,
		output_dir,
		output_file,
		master_accession,
		blast_uniq_hits,
		update_mode=False,      # NEW
		update_root="update",   # NEW
	):
		self.paded_alignment = paded_alignment
		self.master_gff = master_gff
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_accession = master_accession
		self.blast_uniq_hits = blast_uniq_hits

		# NEW
		self.update_mode = update_mode
		self.update_root = update_root

	def _out_dir_abs(self):
		"""
		normal: <tmp_dir>/<output_dir>
		update : <tmp_dir>/update/<output_dir>
		"""
		if self.update_mode:
			return join(self.tmp_dir, self.update_root, self.output_dir)
		return join(self.tmp_dir, self.output_dir)

	def get_gap_ranges(self, sequence):
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1  # 1-based
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges

	def count_gaps_before_position(self, gap_ranges, position):
		count = 0
		for start, end in gap_ranges:
			if end < position:
				count += (end - start + 1)
			elif start <= position <= end:
				count += (position - start + 1)
		return count

	def recalculate_cds_coordinates(self, sequence_id, gap_ranges, cds_list, start_offset):
		adjusted_coords = []

		for cds in cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])

			gaps_before_start = self.count_gaps_before_position(gap_ranges, cds_start)
			gaps_before_end = self.count_gaps_before_position(gap_ranges, cds_end)

			adj_start = cds_start - gaps_before_start + (start_offset - 1)
			adj_end = cds_end - gaps_before_end + (start_offset - 1)

			if [adj_start, adj_end] not in adjusted_coords and adj_start != adj_end:
				adjusted_coords.append([adj_start, adj_end])

		return adjusted_coords

	def get_products_for_range(self, gff_cds_list, coord_range):
		query_start, query_end = int(coord_range[0]), int(coord_range[1])
		results = []

		for cds in gff_cds_list:
			cds_start = int(cds['start'])
			cds_end = int(cds['end'])

			if query_end >= cds_start and query_start <= cds_end:
				overlap_start = max(query_start, cds_start)
				overlap_end = min(query_end, cds_end)
				results.append({
					'start': overlap_start,
					'end': overlap_end,
					'product': cds['product']
				})

		return results

	def load_blast_hits(self):
		acc_dict = {}
		if not self.blast_uniq_hits:
			return acc_dict
		if not os.path.exists(self.blast_uniq_hits):
			print(f"[warn] blast_uniq_hits not found: {self.blast_uniq_hits} (continuing with master as reference)")
			return acc_dict

		for i in open(self.blast_uniq_hits):
			query, ref, score, strand = i.strip().split('\t')
			acc_dict[query] = ref
		return acc_dict

	def find_gaps_in_fasta(self):
		outdir = self._out_dir_abs()
		os.makedirs(outdir, exist_ok=True)

		gff_dict = GffDictionary(self.master_gff).gff_dict
		cds_list = gff_dict['CDS']

		fasta_file_dir = self.paded_alignment
		blast_dict = self.load_blast_hits()

		header = [
			"accession", "master_ref_accession", "reference_accession",
			"aln_start", "aln_end", "cds_start", "cds_end", "product"
		]

		out_path = join(outdir, self.output_file)
		with open(out_path, "w") as out_f:
			out_f.write("\t".join(header) + "\n")

			for fasta_file in os.listdir(fasta_file_dir):
				if not (fasta_file.endswith(".fa") or fasta_file.endswith(".fasta")):
					continue

				calc = CalculateGenomeCoordinates(join(fasta_file_dir, fasta_file), self.master_accession)
				genome_coords = calc.extract_alignment_coordinates()

				for record in SeqIO.parse(join(fasta_file_dir, fasta_file), "fasta"):
					sequence = str(record.seq)
					gaps = self.get_gap_ranges(sequence)

					# start offset: just after the first leading gap block
					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates(record.id, gaps, cds_list, start_offset)

					for each_cords in adjusted:
						product_hits = self.get_products_for_range(cds_list, each_cords)

						if record.id in genome_coords:
							_, genome_cord_start, genome_cord_end = genome_coords[record.id]
						else:
							# fallback if coord mapper didn't return this ID
							genome_cord_start, genome_cord_end = ("NA", "NA")

						ref_acc = blast_dict.get(record.id, self.master_accession)

						for overlap_product in product_hits:
							data = [
								record.id,
								self.master_accession,
								ref_acc,
								str(genome_cord_start),
								str(genome_cord_end),
								str(each_cords[0]),
								str(each_cords[1]),
								overlap_product['product']
							]
							out_f.write("\t".join(data) + "\n")

		print(f"[coords] Wrote: {out_path}")


if __name__ == "__main__":
	parser = ArgumentParser(description='Calculates the genome and CDS coordinates for sequences in padded alignments')
	parser.add_argument('-i', '--paded_alignment', help='Directory containing one or more fasta files.', required=True)
	parser.add_argument('-b', '--tmp_dir', help='Base directory', default="tmp")
	parser.add_argument('-d', '--output_dir', help='Output directory where results are stored', default='Tables')
	parser.add_argument('-o', '--output_file', help='Output file name', default='features.tsv')
	parser.add_argument('-m', '--master_accession', help='Master accession', required=True)
	parser.add_argument('-bh', '--blast_uniq_hits', help='Blast unique hits file', default='tmp/Blast/query_uniq_tophits.tsv')
	parser.add_argument('-g', '--master_gff', help='Master GFF3 file', required=True)

	# NEW: update flag (route outputs under tmp/update/)
	parser.add_argument('--update', action='store_true',
	                    help='If set, write outputs under <tmp_dir>/update/<output_dir> (e.g. tmp/update/Tables).')

	args = parser.parse_args()

	processor = CalculateAlignmentCoordinates(
		paded_alignment=args.paded_alignment,
		master_gff=args.master_gff,
		tmp_dir=args.tmp_dir,
		output_dir=args.output_dir,
		output_file=args.output_file,
		master_accession=args.master_accession,
		blast_uniq_hits=args.blast_uniq_hits,
		update_mode=args.update,   # NEW
		update_root="update",      # NEW
	)
	processor.find_gaps_in_fasta()
