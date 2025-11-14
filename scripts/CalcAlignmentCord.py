import os
from Bio import SeqIO
from os.path import join
from argparse import ArgumentParser
from GffToDictionary import GffDictionary
from CalcGenomeCords import CalculateGenomeCoordinates 

class CalculateAlignmentCoordinates:

	def __init__(self, paded_alignment, master_gff, tmp_dir, output_dir, output_file, master_accession, blast_uniq_hits):
		self.paded_alignment = paded_alignment
		self.master_gff = master_gff
		self.tmp_dir = tmp_dir
		self.output_dir = output_dir
		self.output_file = output_file
		self.master_accession = master_accession
		self.blast_uniq_hits = blast_uniq_hits

	def get_gap_ranges(self, sequence):
		gap_ranges = []
		start = None

		for i, char in enumerate(sequence):
			if char == '-':
				if start is None:
					start = i + 1  # Convert to 1-based indexing
			else:
				if start is not None:
					gap_ranges.append([start, i])
					start = None

		if start is not None:
			gap_ranges.append([start, len(sequence)])

		return gap_ranges

	def count_gaps_before_position(self, gap_ranges, position):
		"""Count how many positions are removed before a given alignment position."""
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

			adj_start = cds_start - gaps_before_start
			adj_end = cds_end - gaps_before_end

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
		for i in open(self.blast_uniq_hits):
			query, ref, score, strand = i.strip().split('\t')
			acc_dict[query] = ref
		return acc_dict

	def find_gaps_in_fasta(self): #, fasta_file_dir, gff_file):
		os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)

		gff_dict = GffDictionary(self.master_gff).gff_dict
		cds_list = gff_dict['CDS']
		fasta_file_dir = self.paded_alignment

		blast_dict = self.load_blast_hits()

		header = ["accession", "master_ref_accession", "reference_accession", "aln_start", "aln_end", "cds_start", "cds_end", "product"]
		with open(join(self.tmp_dir, self.output_dir, self.output_file), "w") as out_f:

			out_f.write("\t".join(header))
			out_f.write("\n")

			for fasta_file in os.listdir(fasta_file_dir):

				calc = CalculateGenomeCoordinates(join(fasta_file_dir, fasta_file), self.master_accession)
				genome_coords = calc.extract_alignment_coordinates()
				for record in SeqIO.parse(join(fasta_file_dir, fasta_file), "fasta"):

					sequence = str(record.seq)
					gaps = self.get_gap_ranges(sequence)
					aligned_length = len(sequence.replace('-', ''))

					# Calculate start offset: just after the first gap
					if gaps and gaps[0][0] == 1:
						start_offset = gaps[0][1] + 1
					else:
						start_offset = 1

					adjusted = self.recalculate_cds_coordinates(record.id, gaps, cds_list, start_offset)

					#print(f">{record.id}", adjusted)
					for each_cords in adjusted:
						product = self.get_products_for_range(cds_list, each_cords)
						master_acc, genome_cord_start, genome_cord_end = genome_coords[record.id]
						for overlap_product in product:
							if record.id in blast_dict:
								data = [record.id, self.master_accession, blast_dict[record.id], str(genome_cord_start), str(genome_cord_end), str(each_cords[0]), str(each_cords[1]), overlap_product['product']]
								out_f.write('\t'.join(data))
								out_f.write("\n")
							else:
								data = [record.id, self.master_accession, self.master_accession, str(genome_cord_start), str(genome_cord_end), str(each_cords[0]), str(each_cords[1]), overlap_product['product']]
								out_f.write('\t'.join(data))	
								out_f.write("\n")				
if __name__ == "__main__":
	parser = ArgumentParser(description='Calculates the genome and cds coordinates for a given sequences')
	parser.add_argument('-i', '--paded_alignment', help='Sequence file directory, it can be single or multiple fasta sequencce files.', required=True)
	parser.add_argument('-b', '--tmp_dir', help='Base directory', default="tmp")
	parser.add_argument('-d', '--output_dir', help='Output directory where processed data and results are stored', default='Tables')
	parser.add_argument('-o', '--output_file', help='Output file name', default='features.tsv')
	parser.add_argument('-m', '--master_accession', help='Master accession', required=True)
	parser.add_argument('-bh', '--blast_uniq_hits', help='Blast unique hits file', default='tmp/Blast/query_uniq_tophits.tsv')
	parser.add_argument('-g', '--master_gff', help='Master GFF3 file', required=True)
	args = parser.parse_args()

	processor = CalculateAlignmentCoordinates(args.paded_alignment, args.master_gff, args.tmp_dir, args.output_dir, args.output_file, args.master_accession, args.blast_uniq_hits)
	processor.find_gaps_in_fasta()

# Example usage:
#find_gaps_in_fasta("NC_001542.aligned_merged_MSA.fasta", "NC_001542.gff3")

