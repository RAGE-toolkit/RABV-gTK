# TING
virus genomic data integration

The script dependencies can be found at the top of the script. For any help, please run the following.
```
python <SCRIPT_NAME.py> -h
```

## GenBank XML downloader
```
$python genBank_downloader.py -b 200
```

## GenBank XML to TSV
```
$python genBank_to_tsv.py -d <DIRECTORY_OUTPUT_FROM_PREVIOUS_RUN> -o <ANY_DIRECTORY_NAME>
```
## Validate TSV
```
$python validate_tsv.py 
```
The script validates the collection date and country. The host part is still pending. The script generates missing_data.tsv file.

