from Bio import SeqIO

class CalculateGenomeCoordinates:
    def __init__(self, paded_alignment, master_accession=None):
        self.paded_alignment = paded_alignment
        self.master_accession = master_accession

    def extract_alignment_coordinates(self):
        records = list(SeqIO.parse(self.paded_alignment, "fasta"))

        # Select the master/reference sequence
        if self.master_accession:
            master_seq = next((r for r in records if r.id == self.master_accession), None)
            if master_seq is None:
                raise ValueError(f"Master sequence with ID '{self.master_accession}' not found.")
        else:
            master_seq = records[0]
            self.master_accession = master_seq.id

        master_alignment = str(master_seq.seq)

        master_coords = []
        master_res_count = 0
        for align_pos, base in enumerate(master_alignment, start=1):
            if base != "-":
                master_res_count += 1
                master_coords.append((align_pos, master_res_count))

        results = {}
        for record in records:
            #if record.id == master_seq.id:
            #    continue  # Skip master sequence

            seq = str(record.seq)
            start_coord = None
            end_coord = None

            for (align_pos, master_coord) in master_coords:
                if seq[align_pos - 1] != "-":
                    if start_coord is None:
                        start_coord = master_coord
                    end_coord = master_coord

            results[record.id] = [
								self.master_accession,
                start_coord if start_coord is not None else "NA",
                end_coord if end_coord is not None else "NA"
            ]

        return results

