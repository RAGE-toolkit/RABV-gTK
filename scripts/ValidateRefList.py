# check there are three columns, 
#1st is accession, 2nd contains master or reference, 3rd is segment (if segmented)
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Validate reference list format.")
parser.add_argument("--ref_list", required=True, help="Path to the reference list TSV file.")
parser.add_argument("--is_segmented", required=True, help="Indicates if the data is segmented ('Y' or 'N').")
params = parser.parse_args()

# read tsv
df = pd.read_csv(params.ref_list, sep="\t", header=None)

if df.shape[1] < 2:
    raise ValueError("Reference list must contain at least two columns: accession and master/reference indicator.")
if params.is_segmented == 'Y' and df.shape[1] < 3:
    raise ValueError("For segmented data, reference list must contain three columns: accession, master/reference indicator, and segment.")

# check master/reference column contains only 'master' or 'reference'
valid_indicators = {'master', 'reference'}
if not df[1].isin(valid_indicators).all():
    raise ValueError("Master/reference indicator column must contain only 'master' or 'reference' values.")

# check for each unique segment value there's a master
if params.is_segmented == 'Y':
    for segment in df[2].unique():
        segment_df = df[df[2] == segment]
        if 'master' not in segment_df[1].values:
            raise ValueError(f"No master found for segment '{segment}'. Each segment must have one master entry.")

print("Reference list validation passed.")
    