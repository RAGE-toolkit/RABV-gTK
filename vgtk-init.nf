
//params defined in nextflow.config override with --params..

// profiles=conda,condaMamba,test,setup_rabv_full
//use conda if not running in conda env alraedy,
// use condaMamba for mamba if installed (check with mamba --version)

scripts_dir     = "${projectDir}/scripts"
// 1. List your script's explicitly defined parameters (keep this in sync!)
def scriptDefinedParams = [
    'tax_id', 'db_name', 'is_segmented', 'extra_info_fill', 'test',
    "scripts_dir", "publish_dir", "email", "ref_list", "bulk_fillup_table","is_flu",
    // Add all parameter names defined above
]

// 2. List common/allowed built-in Nextflow parameters (customize as needed)
//    You might need to add others depending on features you use (e.g., cloud options)
def knownNextflowParams = [
    'help', 'version', 'profile', 'config', 'workDir', 'resume',
    'entry', 'validateParams', 'cached', 'dumpHashes', 'listInputs',
    // Add cloud-specific ones if used: 'awsqueue', 'clusterOptions', etc.
    // Add tower related ones if used: 'tower', 'computeEnvId', etc.
    // Check Nextflow documentation for a more exhaustive list if required
]

// 3. Combine allowed parameters
def allowedParams = (scriptDefinedParams + knownNextflowParams).unique()

// 4. Get all parameters provided via command line and config files
//    The 'params' map holds the merged view of all parameters
def providedParams = params.collect { it.key }

// 5. Find parameters that were provided but are not in the allowed list
def unexpectedParams = providedParams.findAll { ! (it in allowedParams) }

// 6. Error out if any unexpected parameters were found
if (unexpectedParams) {
    // Construct a helpful error message
    def errorMsg = """
    ERROR: Unknown command line parameter(s) provided: ${unexpectedParams.join(', ')}

    Please check your command. Allowed parameters are:
    Script parameters: ${scriptDefinedParams.join(', ')}
    Nextflow options : ${knownNextflowParams.join(', ')} (subset shown)
    """.stripIndent() // Use stripIndent for cleaner multi-line string formatting

    error(errorMsg) // Stop the pipeline
}

process TEST_DEPENDENCIES{
    output:
        path "dependency_test.txt"
    shell:
    '''
    python --version > dependency_test.txt
    pip show biopython >> dependency_test.txt
    pip show pandas >> dependency_test.txt
    pip show numpy >> dependency_test.txt
    pip show openpyxl >> dependency_test.txt
    pip show requests >> dependency_test.txt
    pip show python-dateutil >> dependency_test.txt
    nextalign --version >> dependency_test.txt
    
    '''
}

process VALIDATE_REF_LIST{
    input:
        val ref_list
        val is_segmented
    output:
        path "ref_list_validated.txt"
    shell:
    '''
    # check there are three columns, 
    #1st is accession, 2nd contains master or reference, 3rd is segment (if segmented)
        python !{scripts_dir}/ValidateRefList.py -r !{ref_list} -s !{is_segmented} -o ref_list_validated.txt

    '''
}
process FETCH_GENBANK{
    publishDir "${params.publish_dir}"
    input:
        val TAX_ID
    output:
        path 'GenBank-XML', type: 'dir', emit: gen_bank_XML
    shell:
    '''
    extra=""
    if( [ "!{params.test}" -eq "1" ] )
    then
        extra="--test_run --ref_list !{params.ref_list}"
    fi
    python !{scripts_dir}/GenBankFetcher.py --taxid !{TAX_ID} -b 50 \
             ${extra} -e !{params.email} -o .
    #--update tmp/GenBank-matrix/gB_matrix_raw.tsv is gonna be problematic for this!
    #what's update doing with a tmp dir?

    '''
}

process DOWNLOAD_GFF{
    input:
        val master_acc
        path master_file_opt
    output:
         path "*.gff3"
    shell:
    '''
    MASTER_ARG="!{master_acc}"
    if [ -f "!{master_file_opt}" ]; then
        MASTER_ARG="!{master_file_opt}"
    fi
    python "!{scripts_dir}/DownloadGFF.py" --accession_ids "$MASTER_ARG" -o . -b .

    '''
}


process GENBANK_PARSER{
    publishDir "${params.publish_dir}"
    input:
        path ref_list_path
        path gen_bank_XML
    output:
        path "gB_matrix_validated.tsv" , emit: gb_matrix
        path "sequences.fa", emit: sequences_out
    shell:
    '''
        python !{scripts_dir}/GenBankParser.py -r !{ref_list_path} -d !{gen_bank_XML} -o . -b .
        python !{scripts_dir}/ValidateMatrix.py -o . -a !{projectDir}/assets -b . \
        -g gB_matrix_raw.tsv \
        -m !{projectDir}/assets/host_mapping.tsv -n !{projectDir}/assets/country_mapping.tsv  \
        -c !{projectDir}/assets/m49_country.csv
        
    '''
}
process ADD_MISSING_DATA{
    input:
        path gen_bank_table
    output:
        path "*.tsv", emit: tsvs_out
    when:
        params.extra_info_fill.toBoolean()
    shell:
    '''
        python !{scripts_dir}/AddMissingData.py -b !{gen_bank_table} \
         -t !{params.bulk_fillup_table}  -d . -f !{params.bulk_fillup_table}
    '''
}


process FILTER_AND_EXTRACT{
    input:
        path table_in
        path seqs_in
    output:
        path "query_seq.fa", emit: query_seqs_out
        path "ref_seq.fa", emit: ref_seqs_out
    shell:
    '''

        python !{scripts_dir}/FilterAndExtractSequences.py -b . -o . -r !{params.ref_list} \
         -v !{params.is_segmented} -g !{table_in} -sf !{seqs_in}
    '''
}



process BLAST_ALIGNMENT{
    input:
        path query_seqs
        path ref_seqs
        path gb_matrix
    output:
        path "query_tophits.tsv",type: "file", emit: query_tophits
        path "query_uniq_tophits.tsv", type: "file", emit: query_uniq_tophits
        path "query_uniq_tophit_annotated.tsv", type: "file", optional: true, emit: query_uniq_tophit_annotated
        path "grouped_fasta", type: 'dir', emit: grouped_fasta
        path "ref_seqs", type: 'dir', emit: ref_seqs_dir
        path "ref_seq.fa", type: 'file', emit: ref_seqs_fasta
        path "master_seq", type: 'dir', emit: master_seq_dir
        
    shell:
    '''
        if [ "!{params.is_segmented}" = "Y" ]; then
            python "!{scripts_dir}/BlastAlignment.py" -s Y -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -t . -b . -m !{params.ref_list}  -g !{gb_matrix}
        else
            python "!{scripts_dir}/BlastAlignment.py" -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -b . -t . -m !{params.ref_list} -g !{gb_matrix}
        fi
    '''
}
//workdir files:
//DB                       grouped_fasta  merged_fasta  query_tophits.tsv       ref_seq.fa  sorted_all
//gB_matrix_validated.tsv  master_seq     query_seq.fa  query_uniq_tophits.tsv  ref_seqs    sorted_fasta


//python "${scripts_dir}/NextalignAlignment.py" -m $master_acc #-gff "tmp/Gff/NC_001542.gff3"
process NEXTALIGN_ALIGNMENT{
    input:
        path genbank_matrix
        path grouped_fasta_dir
        path ref_seqs
        path ref_seqs_fasta
        path master_seq_dir
        val master_acc_str
        path master_file_opt
    output:
        path "Nextalign", type: 'dir'
    shell:
    '''
        TARGET_M="!{master_acc_str}"
        if [ -f "!{master_file_opt}" ]; then
             TARGET_M="!{master_file_opt}"
        fi

        python !{scripts_dir}/NextalignAlignment.py  -r !{ref_seqs}  \
         -q !{grouped_fasta_dir} -g !{genbank_matrix} -t . \
         -f !{ref_seqs_fasta} -m "$TARGET_M" -ms !{master_seq_dir} 
    '''
}


//"${scripts_dir}/PadAlignment-1.py" -r "/home3/sk312p/task_dir/projects/VGTK/dev_version-jun-09/TING/alUnc509RefseqsMafftHandModified.fa
process PAD_ALIGNMENT{
    publishDir "${params.publish_dir}"
    input:
        path nextalign_dir 
        val master_acc_str
        path master_file_opt
    output:
        path "*_merged_MSA.fasta", emit: merged_msa
    shell:
    '''
        TARGET_M="!{master_acc_str}"
        if [ -f "!{master_file_opt}" ]; then
             TARGET_M="!{master_file_opt}"
        fi

        python !{scripts_dir}/PadAlignment.py -nd !{nextalign_dir} \
        -m "$TARGET_M" \
        -o . -d . -i !{nextalign_dir}/query_aln --keep_intermediate_files 
    '''
}

process CALC_ALIGNMENT_CORD {
    input:
        path padded_fasta
        path gff_file
        path blast_hits
        val master_acc_str
        path master_file_opt
    output:
        path "features.tsv", emit: features
    shell:
    '''
    TARGET_M="!{master_acc_str}"
    if [ -f "!{master_file_opt}" ]; then
            TARGET_M="!{master_file_opt}"
    fi

    # CalcAlignmentCord expects a directory for -i
    mkdir padded_alignments
    cp !{padded_fasta} padded_alignments/
    
    python !{scripts_dir}/CalcAlignmentCord.py -i padded_alignments \
    -m "$TARGET_M" -g !{gff_file} -bh !{blast_hits} \
    -b . -d . -o features.tsv
    '''
}

process SOFTWARE_VERSION {
    output:
        path "Software_info/software_info.tsv", emit: software_info
    shell:
    '''
    python !{scripts_dir}/SoftwareVersion.py -d . -o Software_info \
     -f software_info.tsv
    '''
}

process GENERATE_TABLES {
    input:
        path gb_matrix
        path blast_hits
        path padded_aln
        path nextalign_dir
    output:
        path "Tables/sequence_alignment.tsv", emit: sequence_alignment
        path "Tables/insertions.tsv", emit: insertions
        path "Tables/host_taxa.tsv", emit: host_taxa
        path "Tables", emit: tables_dir
    shell:
    '''
    python !{scripts_dir}/GenerateTables.py -g !{gb_matrix} \
    -bh !{blast_hits} -p !{padded_aln} -n !{nextalign_dir} \
    -b . -o Tables -e !{params.email}
    '''
}

process CREATE_SQLITE_DB {
    publishDir "${params.publish_dir}"
    input:
        path meta_data
        path features
        path sequence_alignment
        path insertions
        path host_taxa
        path software_info
        path fasta_sequences
    output:
        path "${params.db_name}.db"
    shell:
    '''
    python !{scripts_dir}/CreateSqliteDB.py -m !{meta_data} \
    -rf !{features} -p !{sequence_alignment} \
    -i !{insertions} -ht !{host_taxa} \
    -s !{software_info} -fa !{fasta_sequences} \
    -g !{projectDir}/generic/rabv/Tables/gene_info.csv \
    -mc !{projectDir}/assets/m49_country.csv \
    -mir !{projectDir}/assets/m49_intermediate_region.csv \
    -mr !{projectDir}/assets/m49_region.csv \
    -msr !{projectDir}/assets/m49_sub_region.csv \
    -d !{params.db_name} -b . -o .
    # need to make gene info come from gff file rather than hardcoded rabv one
    '''
}

process VALIDATE_SEGMENT{

    when:
        params.is_segmented.toBoolean()
    input:
        path gb_matrix
        path blast_hits
    output:
        path "gB_matrix_validated_segment.tsv", emit: validated_matrix
    shell:
    '''
    python !{scripts_dir}/ValidateSegment.py \
        -g !{gb_matrix} \
        -s !{blast_hits} \
        -o gB_matrix_validated_segment.tsv
    '''
}

process VALIDATE_STRAIN{

    when:
        params.is_flu == "Y"
    input:
        path gb_matrix
    output:
        path "gB_matrix_validated_strain.tsv", emit: validated_matrix
    shell:
    '''
    python !{scripts_dir}/ValidateStrain.py \
        -g !{gb_matrix} \
        -o gB_matrix_validated_strain.tsv \
        -m !{projectDir}/generic/influenza/serotype_mapping.tsv
    '''
}

process PIVOT_TABLE_SEGMENTS{
    publishDir "${params.publish_dir}"
    when:
        params.is_segmented =="Y" &&
        params.is_flu == "Y"
    input:
        path gb_matrix
    output:
        path "gB_matrix_pivoted_segments.tsv", emit: pivoted_matrix
    shell:
    '''


    python !{scripts_dir}/FluPivotTable.py \
        -g !{gb_matrix} \
        -o gB_matrix_pivoted_segments.tsv
    '''
}

process TEST_SEGMENTED_OUTPUT{
    publishDir "${params.publish_dir}/tests"
    when:
        params.is_segmented == "Y" && params.test == "1"
    input:
        path annotated_blast
        path validated_matrix
        path pivoted_matrix
    output:
        path "test_segmented_results.txt"
    shell:
    '''
    #!/bin/bash
    set -e
    
    echo "=== Testing Segmented Virus Pipeline Output ===" > test_segmented_results.txt
    echo "" >> test_segmented_results.txt
    
    # Test 1: Check annotated BLAST file has 5 columns (including segment)
    echo "Test 1: Checking annotated BLAST file structure..." >> test_segmented_results.txt
    COLS=$(head -1 !{annotated_blast} | awk -F'\t' '{print NF}')
    if [ "$COLS" -eq 5 ]; then
        echo "✓ PASS: Annotated BLAST file has 5 columns (query, reference, score, strand, segment)" >> test_segmented_results.txt
    else
        echo "✗ FAIL: Annotated BLAST file has $COLS columns, expected 5" >> test_segmented_results.txt
        exit 1
    fi
    echo "" >> test_segmented_results.txt
    
    # Test 2: Check segment_validated column exists and has values
    echo "Test 2: Checking segment_validated column in matrix..." >> test_segmented_results.txt
    if head -1 !{validated_matrix} | grep -q "segment_validated"; then
        echo "✓ PASS: segment_validated column exists" >> test_segmented_results.txt
        
        # Count non-empty segment values (excluding "not found")
        SEGMENT_COUNT=$(tail -n +2 !{validated_matrix} | cut -f 67 | grep -v "^$" | grep -v "not found" | wc -l)
        TOTAL_COUNT=$(tail -n +2 !{validated_matrix} | wc -l)
        echo "  - Found $SEGMENT_COUNT records with valid segments out of $TOTAL_COUNT total" >> test_segmented_results.txt
        
        if [ "$SEGMENT_COUNT" -gt 0 ]; then
            echo "✓ PASS: At least some records have segment assignments" >> test_segmented_results.txt
        else
            echo "⚠ WARNING: No records have valid segment assignments" >> test_segmented_results.txt
        fi
    else
        echo "✗ FAIL: segment_validated column not found in matrix" >> test_segmented_results.txt
        exit 1
    fi
    echo "" >> test_segmented_results.txt
    
    # Test 3: Check pivoted matrix structure (flu only)
    if [ -f "!{pivoted_matrix}" ]; then
        echo "Test 3: Checking pivoted segments matrix..." >> test_segmented_results.txt
        HEADER=$(head -1 !{pivoted_matrix})
        
        if echo "$HEADER" | grep -q "Complete_status"; then
            echo "✓ PASS: Pivoted matrix has Complete_status column" >> test_segmented_results.txt
        else
            echo "✗ FAIL: Pivoted matrix missing Complete_status column" >> test_segmented_results.txt
            exit 1
        fi
        
        # Check for segment columns (1-8)
        SEGMENT_COLS=$(echo "$HEADER" | grep -o -E '\t[1-8]\t|\t[1-8]$' | wc -l)
        echo "  - Found $SEGMENT_COLS segment columns" >> test_segmented_results.txt
        
        # Count complete genomes
        COMPLETE=$(tail -n +2 !{pivoted_matrix} | cut -f 10 | grep -c "Complete" || true)
        INCOMPLETE=$(tail -n +2 !{pivoted_matrix} | cut -f 10 | grep -c "Incomplete" || true)
        echo "  - Complete genomes: $COMPLETE" >> test_segmented_results.txt
        echo "  - Incomplete genomes: $INCOMPLETE" >> test_segmented_results.txt
        
        if [ "$((COMPLETE + INCOMPLETE))" -gt 0 ]; then
            echo "✓ PASS: Pivoted matrix contains strain data" >> test_segmented_results.txt
        else
            echo "⚠ WARNING: Pivoted matrix is empty" >> test_segmented_results.txt
        fi
    else
        echo "Test 3: SKIPPED - Pivoted matrix not expected for this run" >> test_segmented_results.txt
    fi
    echo "" >> test_segmented_results.txt
    
    echo "=== All segmented virus tests completed ===" >> test_segmented_results.txt
    cat test_segmented_results.txt
    '''
}

process TEST_NON_SEGMENTED_OUTPUT{
    publishDir "${params.publish_dir}/tests"
    when:
        params.is_segmented == "N" && params.test == "1"
    input:
        path blast_hits
        path gb_matrix
    output:
        path "test_non_segmented_results.txt"
    shell:
    '''
    #!/bin/bash
    set -e
    
    echo "=== Testing Non-Segmented Virus Pipeline Output ===" > test_non_segmented_results.txt
    echo "" >> test_non_segmented_results.txt
    
    # Test 1: Check BLAST file has 4 columns (no segment)
    echo "Test 1: Checking BLAST file structure..." >> test_non_segmented_results.txt
    COLS=$(head -1 !{blast_hits} | awk -F'\t' '{print NF}')
    if [ "$COLS" -eq 4 ]; then
        echo "✓ PASS: BLAST file has 4 columns (query, reference, score, strand)" >> test_non_segmented_results.txt
    else
        echo "✗ FAIL: BLAST file has $COLS columns, expected 4" >> test_non_segmented_results.txt
        exit 1
    fi
    echo "" >> test_non_segmented_results.txt
    
    # Test 2: Check NO segment_validated column in matrix
    echo "Test 2: Verifying no segment column for non-segmented virus..." >> test_non_segmented_results.txt
    if head -1 !{gb_matrix} | grep -q "segment_validated"; then
        echo "⚠ WARNING: segment_validated column found (unexpected for non-segmented virus)" >> test_non_segmented_results.txt
    else
        echo "✓ PASS: No segment_validated column (correct for non-segmented virus)" >> test_non_segmented_results.txt
    fi
    echo "" >> test_non_segmented_results.txt
    
    # Test 3: Count sequences
    echo "Test 3: Checking sequence counts..." >> test_non_segmented_results.txt
    BLAST_COUNT=$(wc -l < !{blast_hits})
    MATRIX_COUNT=$(tail -n +2 !{gb_matrix} | wc -l)
    
    echo "  - BLAST hits: $BLAST_COUNT" >> test_non_segmented_results.txt
    echo "  - Matrix records: $MATRIX_COUNT" >> test_non_segmented_results.txt
    
    if [ "$BLAST_COUNT" -gt 0 ] && [ "$MATRIX_COUNT" -gt 0 ]; then
        echo "✓ PASS: Pipeline processed sequences" >> test_non_segmented_results.txt
    else
        echo "⚠ WARNING: No sequences processed" >> test_non_segmented_results.txt
    fi
    echo "" >> test_non_segmented_results.txt
    
    echo "=== All non-segmented virus tests completed ===" >> test_non_segmented_results.txt
    cat test_non_segmented_results.txt
    '''
}

workflow {

    // check some params are in right form
    // params.is_segmented should be either Y or N
    TEST_DEPENDENCIES()
    if( !(params.is_segmented in ['Y','N']) ){
        error("ERROR: params.is_segmented should be either Y or N")
    }

    VALIDATE_REF_LIST(params.ref_list, params.is_segmented)
    FETCH_GENBANK(params.tax_id)

    def ref_list_file = file(params.ref_list)

    DOWNLOAD_GFF(params.ref_list, ref_list_file)

    GENBANK_PARSER(params.ref_list, FETCH_GENBANK.out.gen_bank_XML)

    if(params.extra_info_fill){
        data=ADD_MISSING_DATA(GENBANK_PARSER.out.gb_matrix)
    }else{
        data=GENBANK_PARSER.out.gb_matrix
    }

    FILTER_AND_EXTRACT(data, 
                        GENBANK_PARSER.out.sequences_out)

    BLAST_ALIGNMENT(FILTER_AND_EXTRACT.out.query_seqs_out,
                    FILTER_AND_EXTRACT.out.ref_seqs_out,
                    data)

    // Add VALIDATE_SEGMENT here
    if (params.is_segmented == 'Y') {
        if (params.is_flu == "Y") {
            VALIDATE_STRAIN(data)
            data = VALIDATE_STRAIN.out.validated_matrix
        }
        VALIDATE_SEGMENT(data, BLAST_ALIGNMENT.out.query_uniq_tophit_annotated)
        // Update 'data' to point to the new validated matrix for downstream steps
        data = VALIDATE_SEGMENT.out.validated_matrix
        
        if (params.is_flu == "Y") {
            PIVOT_TABLE_SEGMENTS(data)
            
            // Run tests for segmented viruses
            if (params.test == "1") {
                TEST_SEGMENTED_OUTPUT(
                    BLAST_ALIGNMENT.out.query_uniq_tophit_annotated,
                    VALIDATE_SEGMENT.out.validated_matrix,
                    PIVOT_TABLE_SEGMENTS.out.pivoted_matrix
                )
            }
        }
    } else {
        // Run tests for non-segmented viruses
        if (params.test == "1") {
            TEST_NON_SEGMENTED_OUTPUT(
                BLAST_ALIGNMENT.out.query_uniq_tophits,
                data
            )
        }
    }


    NEXTALIGN_ALIGNMENT(data,
                        BLAST_ALIGNMENT.out.grouped_fasta,
                        BLAST_ALIGNMENT.out.ref_seqs_dir,
                        FILTER_AND_EXTRACT.out.ref_seqs_out,
                        BLAST_ALIGNMENT.out.master_seq_dir,
                        params.ref_list,
                        ref_list_file)
    PAD_ALIGNMENT(NEXTALIGN_ALIGNMENT.out,
                  params.ref_list,
                  ref_list_file)
    
    CALC_ALIGNMENT_CORD(PAD_ALIGNMENT.out.merged_msa, 
                        DOWNLOAD_GFF.out, 
                        BLAST_ALIGNMENT.out.query_uniq_tophits,
                        params.ref_list,
                        ref_list_file)
                        
    SOFTWARE_VERSION()
    
    GENERATE_TABLES(data, 
                    BLAST_ALIGNMENT.out.query_uniq_tophits, 
                    PAD_ALIGNMENT.out.merged_msa, 
                    NEXTALIGN_ALIGNMENT.out)


                    
    CREATE_SQLITE_DB(data, 
                     CALC_ALIGNMENT_CORD.out.features, 
                     GENERATE_TABLES.out.sequence_alignment, 
                     GENERATE_TABLES.out.insertions, 
                     GENERATE_TABLES.out.host_taxa, 
                     SOFTWARE_VERSION.out.software_info, 
                     GENBANK_PARSER.out.sequences_out
                     )
}


// if you wanted it to do an update run, would have to swap "."s for all the directories for a pre-made one

// notes, the python scripts arguments change a lot -d vs -b vs -o etc
// there's too much directory structure, I'd really strip that all out. 
// It could be base=tmp then every function has a subdir in tmp to keep things clear.
