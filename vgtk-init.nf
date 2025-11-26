
//params defined in nextflow.config override with --params..

// profiles=conda,condaMamba,test,setup_rabv_full
//use conda if not running in conda env alraedy,
// use condaMamba for mamba if installed (check with mamba --version)

scripts_dir     = "${projectDir}/scripts"
// 1. List your script's explicitly defined parameters (keep this in sync!)
def scriptDefinedParams = [
    'tax_id', 'db_name', 'master_acc', 'is_segmented', 'extra_info_fill', 'test',
    "scripts_dir", "publish_dir", "email", "ref_list", "bulk_fillup_table"
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
    output:
         path "*.gff3"
    shell:
    '''
    python "!{scripts_dir}/DownloadGFF.py" --accession_ids "!{master_acc}" -o . -b .

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
        path "grouped_fasta", type: 'dir', emit: grouped_fasta
        path "ref_seqs", type: 'dir', emit: ref_seqs_dir
        path "ref_seq.fa", type: 'file', emit: ref_seqs_fasta
        path "master_seq", type: 'dir', emit: master_seq_dir
        
    shell:
    '''
        if [ "!{params.is_segmented}" = "Y" ]; then
            python "!{scripts_dir}/BlastAlignment.py" -s Y -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -t . -b . -m !{params.master_acc}  -g !{gb_matrix}
        else
            python "!{scripts_dir}/BlastAlignment.py" -f "!{params.ref_list}" -q !{query_seqs} -r !{ref_seqs} \
             -b . -t . -m !{params.master_acc} -g !{gb_matrix}
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
    output:
        path "Nextalign", type: 'dir'
    shell:
    '''
        python !{scripts_dir}/NextalignAlignment.py  -r !{ref_seqs}  \
         -q !{grouped_fasta_dir} -g !{genbank_matrix} -t . \
         -f !{ref_seqs_fasta} -m !{params.master_acc} -ms !{master_seq_dir} 
    '''
}


//"${scripts_dir}/PadAlignment-1.py" -r "/home3/sk312p/task_dir/projects/VGTK/dev_version-jun-09/TING/alUnc509RefseqsMafftHandModified.fa
process PAD_ALIGNMENT{
    publishDir "${params.publish_dir}"
    input:
        path nextalign_dir 
    output:
        path "*_merged_MSA.fasta", emit: merged_msa
        path "*_aligned_padded.fasta", emit: padded_fastas, optional: true
    shell:
    '''
        python !{scripts_dir}/PadAlignment.py  -r !{nextalign_dir}/reference_aln/!{params.master_acc}/!{params.master_acc}.aligned.fasta \
        -o . -d . -i !{nextalign_dir}/query_aln --keep_intermediate_files 
    '''
}

process CALC_ALIGNMENT_CORD {
    input:
        path padded_fasta
        path gff_file
        path blast_hits
    output:
        path "features.tsv", emit: features
    shell:
    '''
    # CalcAlignmentCord expects a directory for -i
    mkdir padded_alignments
    cp !{padded_fasta} padded_alignments/
    
    python !{scripts_dir}/CalcAlignmentCord.py -i padded_alignments \
    -m !{params.master_acc} -g !{gff_file} -bh !{blast_hits} \
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

workflow {

    // check some params are in right form
    // params.is_segmented should be either Y or N
    TEST_DEPENDENCIES()
    if( !(params.is_segmented in ['Y','N']) ){
        error("ERROR: params.is_segmented should be either Y or N")
    }
    FETCH_GENBANK(params.tax_id)

    DOWNLOAD_GFF(params.master_acc)

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

    NEXTALIGN_ALIGNMENT(data,
                        BLAST_ALIGNMENT.out.grouped_fasta,
                        BLAST_ALIGNMENT.out.ref_seqs_dir,
                        FILTER_AND_EXTRACT.out.ref_seqs_out,
                        BLAST_ALIGNMENT.out.master_seq_dir)
    PAD_ALIGNMENT(NEXTALIGN_ALIGNMENT.out)
    
    CALC_ALIGNMENT_CORD(PAD_ALIGNMENT.out.merged_msa, 
                        DOWNLOAD_GFF.out, 
                        BLAST_ALIGNMENT.out.query_uniq_tophits)
                        
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
