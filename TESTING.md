# Pipeline Testing

This document describes the automated tests for the VGTK pipeline.

## Test Processes

The pipeline includes automated validation tests that run when `--test 1` is set in the parameters.

### TEST_SEGMENTED_OUTPUT

**When**: Runs for segmented viruses (e.g., influenza) when `is_segmented="Y"` and `test="1"`

**Validates**:
1. **Annotated BLAST file structure** - Verifies that `query_uniq_tophit_annotated.tsv` has 5 columns:
   - query accession
   - reference accession
   - alignment score
   - strand orientation
   - **segment number** (the key column for segmented viruses)

2. **Segment validation in matrix** - Checks that:
   - `segment_validated` column exists in the GenBank matrix
   - At least some records have valid segment assignments (1-8 for influenza)
   - Reports count of validated vs. total records

3. **Pivoted segments matrix** (flu only) - Verifies that:
   - Pivoted matrix has `Complete_status` column
   - Segment columns (1-8) are present
   - Reports counts of Complete vs. Incomplete genomes

**Output**: `test_segmented_results.txt` in the publish directory

### TEST_NON_SEGMENTED_OUTPUT

**When**: Runs for non-segmented viruses (e.g., RABV) when `is_segmented="N"` and `test="1"`

**Validates**:
1. **BLAST file structure** - Verifies that `query_uniq_tophits.tsv` has 4 columns only:
   - query accession
   - reference accession
   - alignment score
   - strand orientation
   - (NO segment column - this is correct for non-segmented)

2. **No segment column** - Confirms that `segment_validated` column is NOT present in the matrix (as expected)

3. **Sequence processing** - Verifies that:
   - BLAST hits were generated
   - Matrix records were created
   - Reports counts for both

**Output**: `test_non_segmented_results.txt` in the publish directory

## Running Tests

## Python Unit Tests (new)

Pytest-based unit tests now exist for the first two scripts in [scripts](scripts):
- [scripts/AddMissingData.py](scripts/AddMissingData.py)
- [scripts/BlastAlignment.py](scripts/BlastAlignment.py)

Additional unit tests now also cover:
- [scripts/CalcAlignmentCord.py](scripts/CalcAlignmentCord.py)
- [scripts/CalcGenomeCords.py](scripts/CalcGenomeCords.py)

Test files:
- [tests/unit/test_add_missing_data.py](tests/unit/test_add_missing_data.py)
- [tests/unit/test_blast_alignment.py](tests/unit/test_blast_alignment.py)
- [tests/unit/test_calc_alignment_cord.py](tests/unit/test_calc_alignment_cord.py)
- [tests/unit/test_calc_genome_cords.py](tests/unit/test_calc_genome_cords.py)

Deterministic input/expected fixtures:
- [test_data/unit/add_missing_data](test_data/unit/add_missing_data)
- [test_data/unit/blast_alignment](test_data/unit/blast_alignment)
- [test_data/unit/calc_alignment_cord](test_data/unit/calc_alignment_cord)
- [test_data/unit/calc_genome_cords](test_data/unit/calc_genome_cords)

Run locally:

```bash
pytest -q tests/unit
```

CI execution on every push/PR:
- [.github/workflows/python-unit-tests.yml](.github/workflows/python-unit-tests.yml)

### For Segmented Virus (Influenza)
```bash
nextflow run vgtk-init.nf -profile segmented_test
```

Make sure your config includes:
```groovy
params.is_segmented = "Y"
params.is_flu = "Y"
params.test = "1"
```

### For Non-Segmented Virus (RABV)
```bash
nextflow run vgtk-init.nf -profile test
```

Make sure your config includes:
```groovy
params.is_segmented = "N"
params.test = "1"
```

## Test Output

Test results are published to `${params.publish_dir}/tests/` and include:
- ✓ PASS: Test passed successfully
- ✗ FAIL: Test failed (pipeline will exit with error)
- ⚠ WARNING: Potential issue detected (pipeline continues)

## Interpreting Results

### Segmented Virus Tests

**Expected output for successful flu run**:
```
=== Testing Segmented Virus Pipeline Output ===

Test 1: Checking annotated BLAST file structure...
✓ PASS: Annotated BLAST file has 5 columns (query, reference, score, strand, segment)

Test 2: Checking segment_validated column in matrix...
✓ PASS: segment_validated column exists
  - Found 150 records with valid segments out of 200 total
✓ PASS: At least some records have segment assignments

Test 3: Checking pivoted segments matrix...
✓ PASS: Pivoted matrix has Complete_status column
  - Found 8 segment columns
  - Complete genomes: 15
  - Incomplete genomes: 10
✓ PASS: Pivoted matrix contains strain data

=== All segmented virus tests completed ===
```

### Non-Segmented Virus Tests

**Expected output for successful RABV run**:
```
=== Testing Non-Segmented Virus Pipeline Output ===

Test 1: Checking BLAST file structure...
✓ PASS: BLAST file has 4 columns (query, reference, score, strand)

Test 2: Verifying no segment column for non-segmented virus...
✓ PASS: No segment_validated column (correct for non-segmented virus)

Test 3: Checking sequence counts...
  - BLAST hits: 95
  - Matrix records: 95
✓ PASS: Pipeline processed sequences

=== All non-segmented virus tests completed ===
```

## Troubleshooting

### Common Issues

1. **"Annotated BLAST file has 4 columns, expected 5"**
   - The segmented virus pipeline isn't creating the annotated file
   - Check that `is_segmented="Y"` in params
   - Verify BlastAlignment.py has segment_file parameter

2. **"No records have valid segment assignments"**
   - BLAST couldn't match sequences to reference segments
   - Check reference list has segment information
   - Verify BLAST e-value threshold isn't too stringent

3. **"Pivoted matrix is empty"**
   - No sequences passed segment validation
   - Check exclusion criteria in ValidateSegment.py
   - Review BLAST alignment quality

4. **"segment_validated column found (unexpected for non-segmented)"**
   - Pipeline ran segmented steps for non-segmented virus
   - Verify `is_segmented="N"` in config
   - Check workflow conditional logic

## Adding New Tests

To add additional validation tests, create a new process following this template:

```groovy
process TEST_NEW_VALIDATION{
    publishDir "${params.publish_dir}/tests"
    when:
        // Your conditions here
    input:
        // Your inputs
    output:
        path "test_results.txt"
    shell:
    '''
    #!/bin/bash
    set -e
    
    echo "Test description" > test_results.txt
    
    # Your validation logic
    if [ condition ]; then
        echo "✓ PASS: Test passed" >> test_results.txt
    else
        echo "✗ FAIL: Test failed" >> test_results.txt
        exit 1
    fi
    
    cat test_results.txt
    '''
}
```

Then call it in the workflow section with appropriate inputs.
