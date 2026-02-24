# parameters
TAX_ID=${1:-11292} # RABV

scripts_dir="$(dirname "$0")/scripts"
generic_dir="$(dirname "$0")/generic/rabv"
db_name="rabv-jul0425"

db_file="./../../rabv-vgtk/V-gTK/tmp/SqliteDB/rabv-gDB_Dec022025.db"      # set this if updating an existing db, e.g. "rabv-gDB_Dec022025.db"
is_update=0    # 1 for update and 0 for not update

skip_fill=${2:-true}  # Use this variable to control skipping AddMissingData.py
is_segmented=${3:-N}  # segmented virus or not Y for Yes and N for Not
master_acc="NC_001542"

#generic file paths
exclusion_list="${generic_dir}/exclusion_list.txt"
ref_list="${generic_dir}/ref_list_based_on_blast.txt"

# steps to run, O to skip any step
run_genbank_fetcher=0
run_genbank_parser=0
run_download_gff=0
run_curator=0
run_validate_matrix=0
run_filter_extract=0
run_blast_alignment=0
run_nextalign=0
run_pad_alignment=0
run_calc_alignment_cord=0
run_software_version=0
run_generate_tables=0
run_host_taxa=0
run_clade_assignment=1
run_create_db=1

# ----- Run GenBankFetcher -----#
if [ "$run_genbank_fetcher" -eq 1 ]; then
  gb_args=( "--taxid" "$TAX_ID" )
  if [ "$is_update" -eq 1 ]; then
    if [ -z "$db_file" ]; then
      echo "Error: is_update=1 but db_file is empty. Set db_file to an existing .db file."
      exit 1
    fi
    gb_args+=( "--update" "--db" "$db_file" )
  fi

  python "${scripts_dir}/GenBankFetcher.py" "${gb_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: GenBankFetcher.py failed."
    exit 1
  fi

  echo "GenBankFetcher.py completed successfully."
  echo ""
else
  echo "Skipping GenBankFetcher.py"
  echo ""
fi
# ----- End of GenBankFetcher.py ----- 

# ----- Run GenBankParser -----
if [ "$run_genbank_parser" -eq 1 ]; then
  gb_parser_args=( 
    "--exclusion_list" "$exclusion_list"
    "--ref_list"       "$ref_list"
  )
  if [ "$is_update" -eq 1 ]; then
    gb_parser_args+=( "--update" )
  fi

  python "${scripts_dir}/GenBankParser.py" "${gb_parser_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: GenBankParser.py failed."
    exit 1
  fi

  echo "GenBankParser.py completed successfully."
  echo ""
else
  echo "Skipping GenBankParser.py"
  echo ""
fi
# ----- End of GenBankParser.py ----- 

# ----- python DownloadGFF.py ----- #
if [ "$run_download_gff" -eq 1 ]; then
  python "${scripts_dir}/DownloadGFF.py" --accession_ids $master_acc
  if [ $? -ne 0 ]; then
    echo "Error: DownloadGFF.py failed."
    exit 1
  fi
  echo "DownloadGFF.py completed successfully."
  echo ""
else
  echo "Skipping DownloadGFF.py"
  echo ""
fi
# ----- End of DownloadGFF.py ----- 


# ----- Curator.py  curation from MADDOG -----
if [ "$run_curator" -eq 1 ]; then

  # gb_matrix path depends on mode
  if [ "$is_update" -eq 1 ]; then
    cur_gb_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  else
    cur_gb_matrix="tmp/GenBank-matrix/gB_matrix_raw.tsv"
  fi

  cur_curated_file="generic/rabv/curation/MADDOG_curation.tsv"

  cur_args=(
    "--gb_matrix"    "$cur_gb_matrix"
    "--curated_file" "$cur_curated_file"
  )

  if [ "$is_update" -eq 1 ]; then
    cur_args+=( "--update" )
  fi

  python "${scripts_dir}/Curator.py" "${cur_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: Curator.py failed."
    exit 1
  fi
  echo "Curator.py completed successfully."
  echo ""
else
  echo "Skipping Curator.py"
  echo ""
fi
# ----- End of 21483707 Curator.py -----

# ----- Curator.py  curation from article 21483707 -----
if [ "${run_curator:-0}" -eq 1 ]; then

  # gb_matrix path depends on mode
  if [ "$is_update" -eq 1 ]; then
    cur_gb_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  else
    cur_gb_matrix="tmp/GenBank-matrix/gB_matrix_raw.tsv"
  fi

  cur_curated_file="generic/rabv/curation/curation_from_article_21483707.tsv"

  cur_args=(
    "--gb_matrix"    "$cur_gb_matrix"
    "--curated_file" "$cur_curated_file"
  )

  if [ "$is_update" -eq 1 ]; then
    cur_args+=( "--update" )
  fi

  python "${scripts_dir}/Curator.py" "${cur_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: Curator.py failed."
    exit 1
  fi
  echo "Curator.py completed successfully."
  echo ""
else
  echo "Skipping Curator.py"
  echo ""
fi
# ----- End of 21483707 Curator.py -----

# ----- Run ValidateMatrix.py -----
if [ "$run_validate_matrix" -eq 1 ]; then
  if [ "$is_update" -eq 1 ]; then
    gb_matrix_path="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  else
    gb_matrix_path="tmp/GenBank-matrix/gB_matrix_raw.tsv"
  fi

  validate_matrix_args=(
    "--gb_matrix" "$gb_matrix_path"
  )

  if [ "$is_update" -eq 1 ]; then
    validate_matrix_args+=( "--update" )
  fi

  python "${scripts_dir}/ValidateMatrix.py" "${validate_matrix_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: ValidateMatrix.py failed."
    exit 1
  fi
  echo "ValidateMatrix.py completed successfully."
  echo ""
else
  echo "Skipping ValidateMatrix.py"
  echo ""
fi
# ----- End of ValidateMatrix.py -----

# ----- Run FilterAndExtractSequences.py ----- 
if [ "$run_filter_extract" -eq 1 ]; then
  fae_args=(
    "--genbank_matrix" "tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
    "--ref_file"       "generic/rabv/ref_list_based_on_blast.txt"
  )

  # Keep your segmented-virus behavior (only add if your script supports it)
  if [ "$is_segmented" = "Y" ]; then
    fae_args+=( "--segmented" "Y" )
  else
    fae_args=(
    "--genbank_matrix" "tmp/GenBank-matrix/gB_matrix_raw.tsv"
    "--ref_file"       "$ref_list"
  )
  fi

  if [ "$is_update" -eq 1 ]; then
    fae_args+=( "--update" )
  fi

  python "${scripts_dir}/FilterAndExtractSequences.py" "${fae_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: FilterAndExtractSequences.py failed."
    exit 1
  fi
  echo "FilterAndExtractSequences.py completed successfully."
  echo ""
else
  echo "Skipping FilterAndExtractSequences.py"
  echo ""
fi
# ----- End of FilterAndExtractSequences.py ----- 

# ----- BlastAlignment.py (normal vs update mode) ----- 

if [ "$run_blast_alignment" -eq 1 ]; then
  update_gb_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  update_ref_fa="./../../rabv-vgtk/V-gTK/tmp/Sequences/ref_seq.fa"
  update_query_fa="tmp/Update/Sequences/query_seq.fa"

  if [ "$is_update" -eq 1 ]; then
    # Update mode command (as you specified)
    blast_args=(
      "--gb_matrix"   "$update_gb_matrix"
      "--ref_fa"      "$update_ref_fa"
      "--query_fa"    "$update_query_fa"
      "--master_acc"  "$master_acc"
      "--update"
    )

    # Keep segmented option if your BlastAlignment.py supports it in update mode too
    if [ "$is_segmented" = "Y" ]; then
      blast_args+=( "--segmented" "Y" )
    fi

    python "${scripts_dir}/BlastAlignment.py" "${blast_args[@]}"

  else
    # Normal mode (your existing behavior)
    if [ "$is_segmented" = "Y" ]; then
      python "${scripts_dir}/BlastAlignment.py" -s Y -f "${generic_dir}/ref_list.txt" -m "$master_acc"
    else
      python "${scripts_dir}/BlastAlignment.py" -f "${generic_dir}/ref_list.txt" -m "$master_acc"
    fi
  fi

  if [ $? -ne 0 ]; then
    echo "Error: BlastAlignment.py failed."
    exit 1
  fi
  echo "BlastAlignment.py completed successfully."
  echo ""
else
  echo "Skipping BlastAlignment.py"
  echo ""
fi
# ----- End BlastAlignment.py -----

# ----- NextalignAlignment.py (normal vs update mode) -----

if [ "$run_nextalign" -eq 1 ]; then
  # Update-mode paths (edit if yours differ)
  na_gb_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  na_query_dir="tmp/Update/Blast/grouped_fasta/"
  na_ref_dir="tmp/Update/Blast/ref_seqs/"
  na_ref_fa_file="./../../rabv-vgtk/V-gTK/tmp/Sequences/ref_seq.fa"
  na_master_seq_dir="tmp/Update/Blast/master_seq/"
  na_master_ref="$master_acc"

  if [ "$is_update" -eq 1 ]; then
    # Update mode (as you specified)
    na_args=(
      "--gB_matrix"        "$na_gb_matrix"
      "--query_dir"        "$na_query_dir"
      "--ref_dir"          "$na_ref_dir"
      "--ref_fa_file"      "$na_ref_fa_file"
      "--master_seq_dir"   "$na_master_seq_dir"
      "--master_ref"       "$na_master_ref"
      "--update"
    )

    python "${scripts_dir}/NextalignAlignment.py" "${na_args[@]}"
  else
    # Normal mode (as you specified)
    python "${scripts_dir}/NextalignAlignment.py" --master_ref "$master_acc"
  fi

  if [ $? -ne 0 ]; then
    echo "Error: NextalignAlignment.py failed."
    exit 1
  fi
  echo "NextalignAlignment.py completed successfully."
  echo ""
else
  echo "Skipping NextalignAlignment.py"
  echo ""
fi
# ----- End NextalignAlignment.py -----

# ----- PadAlignment (normal vs update mode) -----

if [ "$run_pad_alignment" -eq 1 ]; then
  pad_reference_alignment="generic/rabv/reference_alignment/alUnc509RefseqsMafftHandModified.fa"
  pad_input_dir="tmp/Update/Nextalign/query_aln/"
  pad_master_acc="$master_acc"

  if [ "$is_update" -eq 1 ]; then
    # Update mode (as you specified)
    pad_args=(
      "--reference_alignment" "$pad_reference_alignment"
      "--input_dir"           "$pad_input_dir"
      "--master_acc"          "$pad_master_acc"
      "--update"
    )

    python "${scripts_dir}/PadAlignment.py" "${pad_args[@]}"
  else
    # Normal mode (as you specified)
    python "${scripts_dir}/PadAlignment.py" --reference_alignment "${pad_reference_alignment}" --master_acc "${master_acc}"
  fi

  if [ $? -ne 0 ]; then
    echo "Error: PadAlignment.py failed."
    exit 1
  fi
  echo "PadAlignment.py for query sequence is completed successfully."
  echo ""
else
  echo "Skipping PadAlignment.py"
  echo ""
fi
# ----- End PadAlignment -----

# ----- CalcAlignmentCord (normal vs update mode) -----
if [ "$run_calc_alignment_cord" -eq 1 ]; then
    calc_paded_alignment="tmp/Update/Pad-alignment/"
    calc_blast_uniq_hits="tmp/Update/Blast/query_uniq_tophits.tsv"
    calc_master_gff="./../../rabv-vgtk/V-gTK/tmp/Gff/NC_001542.gff3"

    if [ "$is_update" -eq 1 ]; then
      # Update mode (as you specified)
      calc_args=(
        "--paded_alignment"   "$calc_paded_alignment"
        "--master_accession"  "$master_acc"
        "--blast_uniq_hits"   "$calc_blast_uniq_hits"
        "--master_gff"        "$calc_master_gff"
        "--update"
      )

      python "${scripts_dir}/CalcAlignmentCord.py" "${calc_args[@]}"
    else
      # Normal mode (your current existing command)
      python "${scripts_dir}/CalcAlignmentCord.py" -i "tmp/Pad-alignment/" -m "$master_acc" -g "tmp/Gff/NC_001542.gff3"
    fi

    if [ $? -ne 0 ]; then
      echo "Error: CalcAlignmentCord.py failed."
      exit 1
    fi
    echo "CalcAlignmentCord.py is completed successfully."
    echo ""
else
  echo "Skipping CalcAlignmentCord.py"
  echo ""
fi
# ----- End CalcAlignmentCord -----

# ----- SoftwareVersion.py (normal vs update mode) -----
if [ "$run_software_version" -eq 1 ]; then
  if [ "$is_update" -eq 1 ]; then
    python "${scripts_dir}/SoftwareVersion.py" --tmp_dir "tmp/Update"
  else
    python "${scripts_dir}/SoftwareVersion.py"
  fi

  if [ $? -ne 0 ]; then
    echo "Error: SoftwareVersion.py failed."
    exit 1
  fi
  echo "SoftwareVersion.py completed successfully."
  echo ""
else
  echo "Skipping SoftwareVersion.py"
  echo ""
fi
# ----- End SoftwareVersion.py -----

# ----- GenerateTables.py (normal vs update mode) -----
if [ "$run_generate_tables" -eq 1 ]; then
  gt_genbank_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  gt_blast_hits="tmp/Update/Blast/query_uniq_tophits.tsv"
  gt_nextalign_dir="tmp/Update/Nextalign/"
  gt_paded_aln="tmp/Update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa"

  if [ "$is_update" -eq 1 ]; then
    gt_args=(
      "--genbank_matrix" "$gt_genbank_matrix"
      "--blast_hits"     "$gt_blast_hits"
      "--nextalign_dir"  "$gt_nextalign_dir"
      "--paded_aln"      "$gt_paded_aln"
      "--update"
    )
    python "${scripts_dir}/GenerateTables.py" "${gt_args[@]}"
  else
    python "${scripts_dir}/GenerateTables.py"
  fi

  if [ $? -ne 0 ]; then
    echo "Error: GenerateTables.py failed."
    exit 1
  fi
  echo "GenerateTables.py completed successfully."
  echo ""
else
  echo "Skipping GenerateTables.py"
  echo ""
fi
# ----- End GenerateTables.py -----

# ----- HostTaxaTable.py (normal vs update mode) -----
if [ "$run_host_taxa" -eq 1 ]; then
  host_gb_matrix="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  host_names="tmp/Update/Taxa/names.dmp"
  host_nodes="tmp/Update/Taxa/nodes.dmp"

  if [ "$is_update" -eq 1 ]; then
    host_args=(
      "--gb_matrix" "$host_gb_matrix"
      "--names"     "$host_names"
      "--nodes"     "$host_nodes"
      "--update"
    )
    python "${scripts_dir}/HostTaxaTable.py" "${host_args[@]}"
  else
    python "${scripts_dir}/HostTaxaTable.py"
  fi

  if [ $? -ne 0 ]; then
    echo "Error: HostTaxaTable.py failed."
    exit 1
  fi
  echo "HostTaxaTable.py completed successfully."
  echo ""
else
  echo "Skipping HostTaxaTable.py"
  echo ""
fi
# ----- End HostTaxaTable.py -----

# ----- CladeAssignment.py (normal vs update mode) -----
if [ "$run_clade_assignment" -eq 1 ]; then
  ca_ref_aln="generic/rabv/tree/ref_plus_am3ca_am5.fa"
  ca_ref_tree="generic/rabv/tree/ref_tree_am3c_am5.treefile"
  ca_taxon_major="generic/rabv/reference_clades/ref_major_clades.tsv"
  ca_taxon_minor="generic/rabv/reference_clades/ref_minor_clades.tsv"
  ca_threads=6

  # normal vs update inputs
  if [ "$is_update" -eq 1 ]; then
    ca_query="tmp/Update/Pad-alignment/alUnc509RefseqsMafftHandModified.fa"
    ca_meta_data="tmp/Update/GenBank-matrix/gB_matrix_raw.tsv"
  else
    ca_query="tmp/Pad-alignment/alUnc509RefseqsMafftHandModified.fa"
    ca_meta_data="tmp/GenBank-matrix/gB_matrix_raw.tsv"
  fi

  ca_args=(
    "--ref-aln"     "$ca_ref_aln"
    "--ref-tree"    "$ca_ref_tree"
    "--query"       "$ca_query"
    "--taxon-major" "$ca_taxon_major"
    "--taxon-minor" "$ca_taxon_minor"
    "--threads"     "$ca_threads"
    "--meta-data"   "$ca_meta_data"
    "--steps"       "all"
    "--skip-mafft"
  )

  # update flag only in update mode
  if [ "$is_update" -eq 1 ]; then
    ca_args+=( "--update" )
  fi

  python "scripts/CladeAssignment.py" "${ca_args[@]}"
  if [ $? -ne 0 ]; then
    echo "Error: CladeAssignment.py failed."
    exit 1
  fi
  echo "CladeAssignment.py completed successfully."
  echo ""
else
  echo "Skipping CladeAssignment.py"
  echo ""
fi
# ----- End CladeAssignment.py -----