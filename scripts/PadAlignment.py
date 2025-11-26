import os
import shutil
import argparse
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq

class PadAlignment:
	def __init__(self, reference_alignment, input_dir, base_dir, output_dir, keep_intermediate_files, new_outputfile=False):
		self.reference_alignment = reference_alignment
		self.input_dir = input_dir
		self.base_dir = base_dir
		self.output_dir = output_dir
		self.keep_intermediate_files = keep_intermediate_files
		self.new_outputfile = new_outputfile

	def insert_gaps(self, reference_aligned, subalignment_seqs):
		ref_with_gaps_list = list(reference_aligned)
		updated_sequences = []
		for seq_record in subalignment_seqs:
			sequence = list(str(seq_record.seq))
			gapped_sequence = []
			seq_idx = 0
			for char in ref_with_gaps_list:
				if char == '-':
					gapped_sequence.append('-')
				else:
					if seq_idx < len(sequence):
						gapped_sequence.append(sequence[seq_idx])
						seq_idx += 1
			gapped_seq_str = ''.join(gapped_sequence)
			seq_record.seq = Seq(gapped_seq_str)
			updated_sequences.append(seq_record)
		return updated_sequences

	def process_master_alignment(self, reference_alignment_file, input_dir, base_dir, output_dir, keep_intermediate_files=False):
		master_alignment = SeqIO.to_dict(SeqIO.parse(reference_alignment_file, "fasta"))
		merged_sequences = []
		for ref_id, ref_record in master_alignment.items():
			ref_aligned = ref_record.seq
			ref_id = ref_id.split('.')[0]
			subalignment_file = os.path.join(input_dir, f"{ref_id}/{ref_id}.aligned.fasta")
			if os.path.exists(subalignment_file):
				print(f"Processing subalignment for {ref_id} using {subalignment_file}")
				subalignment_seqs = list(SeqIO.parse(subalignment_file, "fasta"))
				updated_seqs = self.insert_gaps(ref_aligned, subalignment_seqs)
				os.makedirs(join(output_dir), exist_ok=True)
				output_file = os.path.join(output_dir, f"{ref_id}_aligned_padded.fasta")
				with open(join(output_file), "w") as output_handle:
					SeqIO.write(updated_seqs, output_handle, "fasta")
					print(f"Saved updated alignment to {output_file}")
				merged_sequences.extend(updated_seqs)
			else:
				print(f"Subalignment file {subalignment_file} not found. Skipping {ref_id}.")

		if merged_sequences:
			merged_output_file = os.path.join(base_dir, output_dir, os.path.basename(reference_alignment_file).replace(".fasta", "_merged_MSA.fasta"))
			os.makedirs(join(base_dir, output_dir), exist_ok=True)
			with open(merged_output_file, "w") as merged_output_handle:
				SeqIO.write(merged_sequences, merged_output_handle, "fasta")
				print(f"Saved merged alignment to {merged_output_file}")
		else:
			print(f"No sequences found for {reference_alignment_file}, skipping merged file creation.")

		if not keep_intermediate_files:
			for ref_id in master_alignment:
				padded_file = os.path.join(output_dir, f"{ref_id.split('.')[0]}_aligned_padded.fasta")
				if os.path.exists(padded_file):
					os.remove(padded_file)
					print(f"Deleted intermediate file {padded_file}")
			shutil.rmtree(output_dir)

	def find_fasta_file(self, input_dir,new_outputfile=False):
		directory = join(self.base_dir, self.output_dir)
		if new_outputfile:
			return os.path.join(directory, "new_output.fasta")
		else:
			for file in os.listdir(directory):
				if file.endswith(".fasta") or file.endswith(".fa"):
					return os.path.join(directory, file)
			return None

	def remove_redundant_sequences(self):
		
		input_file = self.find_fasta_file(join(self.base_dir, self.output_dir),self.new_outputfile) 
		unique_records = {}
		for record in SeqIO.parse(input_file, "fasta"):
			accession = record.id.split('|')[0] if '|' in record.id else record.id
			if accession not in unique_records:
				unique_records[accession] = record
		with open(input_file, "w") as output_handle:
			SeqIO.write(unique_records.values(), output_handle, "fasta")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Insert gaps from master alignment into corresponding subalignments.")
	parser.add_argument("-r", "--reference_alignment", help="Path to master alignment file (FASTA format).", required=True)
	parser.add_argument("-i", "--input_dir", help="Directory containing subalignment files (Nextalign output).", default="tmp/Nextalign/query_aln")
	parser.add_argument("-d", "--base_dir", help="Base directory.", default="tmp")
	parser.add_argument("-o", "--output_dir", help="Directory to save padded subalignments and merged files.", default="Pad-alignment")
	parser.add_argument("--keep_intermediate_files", action="store_true", help="Keep intermediate files (padded subalignment). Default: disabled (files will be removed).")
	parser.add_argument("-n","--new_outputfile", action="store_true", help="New output file name for the final merged alignment.")
 
	args = parser.parse_args()

	processor = PadAlignment(args.reference_alignment, args.input_dir, args.base_dir, args.output_dir, args.keep_intermediate_files)
	processor.process_master_alignment(
		reference_alignment_file=args.reference_alignment,
		input_dir=args.input_dir,
		base_dir=args.base_dir,
		output_dir=args.output_dir,
		keep_intermediate_files=args.keep_intermediate_files
	)
	processor.remove_redundant_sequences()

