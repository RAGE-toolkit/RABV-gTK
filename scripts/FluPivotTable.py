import csv
from collections import defaultdict
from os.path import join, exists
from os import makedirs
from argparse import ArgumentParser

# Provide required segment list 
required_segments = ['1', '2', '3', '4', '5', '6', '7', '8']
print(f"Required segments: {required_segments}")

# Function to pivot the data and check if Complete genome is available for each strain
def pivot_data(input_file, output_file, required_segments):
    data_wide = defaultdict(lambda: defaultdict(list))
    segments = set()
    complete_count = 0
    incomplete_count = 0

    with open(input_file, mode='r', newline='') as csv_file:
        reader = csv.DictReader(csv_file, delimiter="\t")

        if not {'Parsed_strain', 'segment_validated', 'primary_accession', 'exclusion'}.issubset(reader.fieldnames):
            raise ValueError("Input file is missing one or more required columns: 'Parsed_strain', 'segment_validated', 'primary_accession', 'exclusion'")

        for row in reader:
            # Skip rows with exclusion
            if row['exclusion'].strip():
                continue

            parsed_strain = row['Parsed_strain']
            segment_value = row['segment_validated'].strip()
            if segment_value:
                segment = str(int(float(segment_value)))
            else:
                print(f"[WARNING] Empty segment_validated for strain: {row['primary_accession']}, {parsed_strain}")
                segment = 'Not available'

            primary_accession = row['primary_accession']

            if segment in required_segments:
                data_wide[parsed_strain][segment].append(primary_accession)
                segments.add(segment)

    segments = sorted(segments)
    # print(f"Segments found: {segments}")

    with open(output_file, mode='w', newline='') as csv_output_file:
        fieldnames = ['Parsed_strain'] + segments + ['Complete_status']
        writer = csv.DictWriter(csv_output_file, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for strain, loci_dict in data_wide.items():
            row = {'Parsed_strain': strain}
            for segment in segments:
                row[segment] = ','.join(loci_dict.get(segment, []))

            complete = all(bool(loci_dict.get(segment)) for segment in required_segments)
            row['Complete_status'] = 'Complete' if complete else 'Incomplete'

            if complete:
                complete_count += 1
            else:
                incomplete_count += 1

            writer.writerow(row)

    print(f"Number of Complete genomes: {complete_count}")
    print(f"Number of Incomplete genomes: {incomplete_count}")

if __name__ == "__main__":
    parser = ArgumentParser(description='Compile the per strain accession numbers to validate a complete genome. Requires deduplicated Genbank matrix')
    parser.add_argument('-g', '--gb_matrix', help='Genbank matrix (deduplicated) to be processed', default='./tmp/Validate-matrix/gB_matrix_validated.tsv')
    parser.add_argument('-d', '--tmp_dir', help='Path to output folder', default='tmp/Validate-matrix')
    parser.add_argument('-o', '--output_file', help='Output file name for the pivot strains matrix', default='gB_matrix_cast_strains.tsv')
    args = parser.parse_args()

    if not exists(args.tmp_dir):
        makedirs(args.tmp_dir)

    output_path = join(args.tmp_dir, args.output_file)
    pivot_data(args.gb_matrix, output_path, required_segments)

    print(f"Pivoted data has been written to {output_path}")
