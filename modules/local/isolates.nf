process mlstSearch {
    label "mlst"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}.mlst.json")
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta
    mlst input_genome.fasta --label ${meta.alias} --json ${meta.alias}.mlst.json
    """
}


process create_sourmash_picklist {
    label "sourmash"
    memory "8 GB"
    
    input:
    path assemblies_to_exclude_file
    
    output:
    path "sourmash_picklist.csv", emit: picklist
    path "sourmash_picklist_excluded.txt", emit: excluded_list

    script:
    """
    workflow-glue sourmash_picklist \
        --fungi-db /sourmash_db/ncbi-euks-fungi-2025.01.dna.k=51.sig.zip \
        --bacteria-db /sourmash_db/gtdb-reps-rs226-k51.dna.zip \
        --exclude-file ${assemblies_to_exclude_file} \
        --output sourmash_picklist.csv
    """
}


process sourmash_species_id {
    label "sourmash"
    memory "8 GB"
    
    input:
    tuple val(meta), path(assembly)
    path picklist
    
    output:
    tuple val(meta), path("${meta.alias}_sourmash_taxonomy.csv"), emit: taxonomy
    
    script:
    """
    sourmash sketch dna --name-from-first -p k=51,scaled=10000 -o ${meta.alias}.sig ${assembly}
    sourmash gather \
        --threshold-bp 25000 \
        -o ${meta.alias}_gather.csv \
        --picklist ${picklist}:md5:md5 \
        ${meta.alias}.sig \
        /sourmash_db/ncbi-euks-fungi-2025.01.dna.k=51.sig.zip \
        /sourmash_db/gtdb-reps-rs226-k51.dna.zip

    if [ -f "${meta.alias}_gather.csv" ] && [ -s "${meta.alias}_gather.csv" ]; then
        mkdir -p ${meta.alias}_annotated_dir
        sourmash tax annotate --gather-csv ${meta.alias}_gather.csv --taxonomy /sourmash_db/lineages_combined.csv -o ${meta.alias}_annotated_dir
        mv ${meta.alias}_annotated_dir/${meta.alias}_gather.with-lineages.csv ${meta.alias}_sourmash_taxonomy.csv
    else
        echo "intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,filename,name,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp,ksize,moltype,scaled,query_n_hashes,query_abundance,query_containment_ani,match_containment_ani,average_containment_ani,max_containment_ani,potential_false_negative,n_unique_weighted_found,sum_weighted_found,total_weighted_hashes,lineage" > ${meta.alias}_sourmash_taxonomy.csv
    fi
    """
}


process getPointfinderSpecies {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("${meta.alias}.mlst.json"), path("${meta.alias}_sourmash_taxonomy.csv")
        path db_mapping_csv
    output:
        tuple val(meta), stdout
    script:
    def mapping_arg = db_mapping_csv.name != 'OPTIONAL_FILE' ? "--mapping_csv ${db_mapping_csv}" : ""
    """
    workflow-glue pointfinder_species \
        --mlst_json ${meta.alias}.mlst.json \
        --sourmash_csv ${meta.alias}_sourmash_taxonomy.csv \
        ${mapping_arg}
    """
}


process resfinder {
    label "amr"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("input_genome.fasta.gz"), val(species)
        path resfinder_db, stageAs: "resfinder_db_dir"
        path pointfinder_db, stageAs: "pointfinder_db_dir"
        val resfinder_threshold
        val resfinder_coverage
        val pointfinder_ignore_indels
        val pointfinder_ignore_stop_codons
    output:
        tuple val(meta), path("${meta.alias}_resfinder_results"), val(species)
    script:

    String ignore_indels_arg = pointfinder_ignore_indels ? "--ignore_indels" : ""
    String ignore_stop_codons_arg = pointfinder_ignore_stop_codons ? "--ignore_stop_codons" : ""

    """
    gunzip -c input_genome.fasta.gz | sed '/^>/ s/ .*//' > input_genome.fasta

    DB_RESFINDER=""
    DB_POINTFINDER=""
    # Check if databases exist and are not OPTIONAL_FILE
    if [[ -d "resfinder_db_dir" && "${resfinder_db.name}" != "OPTIONAL_FILE" ]]; then
        DB_RESFINDER="--db_path_res resfinder_db_dir"
    fi
    if [[ -d "pointfinder_db_dir" && "${pointfinder_db.name}" != "OPTIONAL_FILE" ]]; then
        DB_POINTFINDER="--db_path_point pointfinder_db_dir"
    fi

    python -m resfinder \
        -ifa input_genome.fasta \
        -o ${meta.alias}_resfinder_results \
        -j ${meta.alias}_resfinder_results/${meta.alias}_resfinder.json \
        \$DB_RESFINDER \
        \$DB_POINTFINDER \
        -s "${species}" \
        -l ${resfinder_coverage} \
        -t ${resfinder_threshold} \
        --acquired \
        --point \
        --disinfectant \
        --nanopore \
        ${ignore_indels_arg} \
        ${ignore_stop_codons_arg}

    # Extract available PointFinder genes for this species (for report)
    POINTFINDER_DIR="\${CGE_RESFINDER_RESPOINT_DB}/${species}"
    
    if [ -d "\${POINTFINDER_DIR}" ]; then
        # List all .fsa files, extract gene names, and exclude the species file itself
        ls "\${POINTFINDER_DIR}"/*.fsa 2>/dev/null | while read -r file; do
            gene=\$(basename "\$file" .fsa)
            # Skip the main species file
            if [ "\$gene" != "${species}" ]; then
                echo "\$gene"
            fi
        done | sort > "${meta.alias}_resfinder_results/pointfinder_genes_analysed.txt"
    else
        # If species directory doesn't exist, create empty file
        touch "${meta.alias}_resfinder_results/pointfinder_genes_analysed.txt"
    fi
    """
}

process serotyping {
    label "seqsero2"
    cpus 1
    memory "3 GB"
    errorStrategy 'ignore'
    input: 
        tuple val(meta), path("input_genome.fasta.gz"), val(species)
    output:
        tuple val(meta), path("${meta.alias}.serotype_results.tsv")
    script:
    """
    gunzip -c input_genome.fasta.gz > input_genome.fasta

    SeqSero2_package.py \
    -m k \
    -t '4' \
    -i input_genome.fasta \
    -p 1 \
    -b 'mem' \
    -d  output \
    -n ${meta.alias}

    cp -r output/SeqSero_result.tsv "${meta.alias}.serotype_results.tsv"
    """
}   


workflow run_isolates {
   take:
      consensus
      assemblies_to_exclude
      resfinder_threshold
      resfinder_coverage
      pointfinder_ignore_indels
      pointfinder_ignore_stop_codons
      resfinder_db
      pointfinder_db
      sourmash_pointfinder_mapping
   main:
        // MLST species ID 
        mlst_results = mlstSearch(consensus)
        // Sourmash species ID
        sourmash_picklist = create_sourmash_picklist(file(assemblies_to_exclude))
        sourmash_results = sourmash_species_id(consensus, sourmash_picklist.picklist)
        // Prep AMR input
        species_input = mlst_results.join(sourmash_results.taxonomy)
        pointfinder_species = getPointfinderSpecies(
            species_input, 
            sourmash_pointfinder_mapping.first()
        ).map{ meta, species -> [meta, species.trim()] }
        resfinder_input = consensus.join(pointfinder_species)
        // AMR 
        amr_results = resfinder(
            resfinder_input,
            resfinder_db,
            pointfinder_db,
            resfinder_threshold, 
            resfinder_coverage,
            pointfinder_ignore_indels,
            pointfinder_ignore_stop_codons
        )
        // Serotyping
        serotype = serotyping(resfinder_input
            | filter { meta, fasta, species -> species == "salmonella" }
        )
        
   emit:
      amr = amr_results.map{meta, amr, species -> [meta, amr]}
      mlst = mlst_results
      serotype = serotype
      taxonomy = sourmash_results.taxonomy
      sourmash_excluded_genomes = sourmash_picklist.excluded_list
}
