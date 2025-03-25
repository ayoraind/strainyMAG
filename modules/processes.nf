process FLYE_BUILD_ASSEMBLY {
    publishDir "${params.output_dir}/flye_out/${meta}_FLYE", mode:'copy'
    tag "flye on $meta"
    memory { 4.GB * task.attempt }

    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta")   , emit: fasta_ch
    tuple val(meta), path("*.gfa")     , emit: flye_gfa_ch
    tuple val(meta), path("*.gv")      , emit: gv_ch
    tuple val(meta), path("*.txt")     , emit: txt_ch
    tuple val(meta), path("*.log")     , emit: log_ch
    tuple val(meta), path("*.json")    , emit: json_ch
    path "versions.yml"                , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """

    flye $mode $reads --keep-haplotypes --meta --no-alt-contigs -i 0 --out-dir . --threads 1 $args

    mv assembly.fasta ${prefix}.fasta
    mv flye.log ${prefix}.flye.log
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv params.json ${prefix}.params.json
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

}

process STRAINY_SPLIT_UNITIGS {
    tag "split unitigs from $meta assemblies"
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/strainy_out/${meta}_strainy_split", mode:'copy'
    
    
    errorStrategy { task.attempt <= 5 ? "retry" : "ignore" }
    maxRetries 5
    
    input:
    tuple val(meta), path(reads), path(assembly_gfa)

    output:
    tuple val(meta), path("preprocessing_data/${meta}_long_unitigs_split.bam"),     emit: long_bam_ch
    tuple val(meta), path("preprocessing_data/${meta}_long_unitigs_split.bam.bai"), emit: long_bai_ch
    tuple val(meta), path("preprocessing_data/${meta}_long_unitigs_split.gfa"),     emit: long_unitig_ch
    tuple val(meta), path("preprocessing_data/${meta}_gfa_converted.fasta"),        emit: strainy_gfa_ch
    path "versions.yml",                                                            emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    strainy.py -g $assembly_gfa -q $reads  --unitig-split-length 50 -o . -t 1 -m nano --only_split True
    
    mv preprocessing_data/long_unitigs_split.bam preprocessing_data/${meta}_long_unitigs_split.bam
    mv preprocessing_data/long_unitigs_split.bam.bai preprocessing_data/${meta}_long_unitigs_split.bam.bai
    mv preprocessing_data/long_unitigs_split.gfa preprocessing_data/${meta}_long_unitigs_split.gfa
    mv preprocessing_data/gfa_converted.fasta preprocessing_data/${meta}_gfa_converted.fasta
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         strainy: \$( strainy --version  )
    END_VERSIONS
    """
}

process METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS {
    tag "$meta"
    
    publishDir "${params.output_dir}/bam_contig_depth", mode:'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt"), emit: depth_ch
    path "versions.yml"           , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    export OMP_NUM_THREADS=$task.cpus

    jgi_summarize_bam_contig_depths \\
        --outputDepth ${meta}.txt \\
        $args \\
        $bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}

process BINNING {
    tag "binning from $meta assembly graphs"
    
    publishDir "${params.output_dir}/bins", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(gfa_converted_fasta), path(depth)

    output:
    path("*.fa"),                   emit: bin_ch
    path "versions.yml",            emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    metabat2 $args -i $gfa_converted_fasta -a $depth -o ${meta}_bin -t $task.cpus 
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}


process MAG_QA {
    
    publishDir "${params.output_dir}", mode:'copy'
    
    
    errorStrategy { task.attempt <= 3 ? "retry" : "finish" }
    maxRetries 3
    
    input:
    path(collected_bin_fa)

    output:
    path("qa_bins/diamond_output"),        emit: diamond_dir_ch
    path("qa_bins/protein_files"),         emit: protein_dir_ch
    path("qa_bins/quality_report.tsv"),    emit: quality_report_ch
    path("qa_bins/checkm2.log"),           emit: log_ch
    path("versions.yml"),                  emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    checkm2.py predict -i $collected_bin_fa -o qa_bins --database_path ${baseDir}/CheckM2_database/uniref100.KO.1.dmnd -t $task.cpus  
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         checkm2: \$( checkm2.py 2>&1 | grep '...:::' | sed 's/.*CheckM2 v//;s/ .*//' )
    END_VERSIONS
    """
}

process UNITIG_PER_BIN {
    
    publishDir "${params.output_dir}/unitigs", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), val(meta2), path(bin_fa)

    output:
    tuple val(meta), val(meta2), path("${meta2}.lst"),  		    emit: unitig_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    cat ${meta2}.fa | grep '>' | cut -c 2- > ${meta2}.lst 
         
    """
}

process SUBSET_BAM {
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/bams", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(long_unitigs_split_bam), path(long_unitigs_split_bam_bai), val(meta2), path(unitig_lst)

    output:
    tuple val(meta), val(meta2), path("${meta2}.bam"),  		   emit: bam_ch
    path("${meta2}.bam"),						   emit: path_bam_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    split_bam.py -bam $long_unitigs_split_bam -l $unitig_lst -o ${meta2}.bam 
         
    """
}

process INDEX_BAM {
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/bams", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "ignore" }
    maxRetries 2
    
    input:
    tuple val(meta), val(meta2), path(bam)

    output:
    tuple val(meta), val(meta2), path("${meta2}.bam.bai"),  		   emit: bam_bai_ch
    path "versions.yml",                                  		   emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    samtools index ${meta2}.bam  
    
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS     
    """
}

process MAKE_SPLIT_FA {
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/split_fa", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(long_unitigs_split)

    output:
    tuple val(meta), path("${meta}.gfa_splitted.fasta"),  		   emit: fasta_split_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    awk '//^S{{print ">"\$2"\\n"\$3}}' ${meta}_long_unitigs_split.gfa | fold > ${meta}.gfa_splitted.fasta 
        
    """
}

process INDEX_FA {
    
    publishDir "${params.output_dir}/strainy_out/${meta}_strainy_split/preprocessing_data", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(gfa_converted_fasta)

    output:
    tuple val(meta), path("${meta}_gfa_converted.fasta.fai"),  		   emit: fai_ch
    path "versions.yml",                                                   emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    samtools faidx ${meta}_gfa_converted.fasta  
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS     
    """
}

process CALL_SNP {
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/SNPs", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(fai), path(gfa_converted_fasta), val(meta2), path(bam), path(bai) 

    output:
    tuple val(meta), val(meta2), path("${meta2}/merge_output.vcf.gz"),  		   emit: call_snp_ch
    path "versions.yml",                                                  	           emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    run_clair3.sh --bam_fn=${meta2}.bam --ref_fn=${meta}_gfa_converted.fasta --output=${meta2} --threads=30 --platform=ont --include_all_ctgs --model_path=\${CONDA_PREFIX}/bin/models/r941_prom_hac_g360+g422 --chunk_size=50000  --snp_min_af=0.15 --no_phasing_for_fa  
    
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh --version 2>&1) | sed 's/.*v\\([0-9.]*\\).*/\\1/')
    END_VERSIONS     
    """
}

process UNZIP_SNP {
    
    publishDir "${params.output_dir}/SNPs/${meta2}", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), val(meta2), path(zipped_vcf)

    output:
    tuple val(meta), val(meta2), path("snp.vcf"),  		   emit: unzip_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """  
    
    gunzip -c $zipped_vcf > snp.vcf  
         
    
    """
}


process CONNECT_UNITIGS {
    
    publishDir "${params.output_dir}/unitigs_connected", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(long_unitigs_split), val(meta2), path(bin_lst), path(unitigs_dir)

    output:
    tuple val(meta), val(meta2), path("${meta2}.lst"),  		   emit: connected_unitig_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """  
    
    connect_bin.py ${meta2}.lst $unitigs_dir ${meta}_long_unitigs_split.gfa ${meta2}.lst     
    
    """
}


process SUBSET_GFA {
    
    publishDir "${params.output_dir}/gfa_sub", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    tuple val(meta), path(long_unitigs_split), val(meta2), path(unitig_connected_lst)

    output:
    tuple val(meta), val(meta2), path("${meta2}.gfa"),  		   emit: subset_gfa_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """  
    
    split_gfa.py -g $long_unitigs_split -l ${meta2}.lst -outfile ${meta2}.gfa     
    
    """
}


process STRAINY {
    tag "strainy from $meta2"
    
    cache 'lenient'
    
    publishDir "${params.output_dir}/strainy_final/${meta2}/", mode:'copy'
    
    
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5
    
    input:
    tuple val(meta), path(reads), val(meta2), path(subsetted_bam), path(subsetted_bai), path(subset_gfa), path(unzipped_vcf) 

    output:
    tuple val(meta), val(meta2), path("strainy_final.gfa"),             emit: strainy_final_ch
    tuple val(meta), val(meta2), path("*"),                             emit: all_bai_ch
    path "versions.yml",                         	                emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    strainy.py -g $subset_gfa -q $reads --unitig-split-length 0 -o . -t $task.cpus -m nano -b $subsetted_bam --snp $unzipped_vcf
    
    
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         strainy: \$( echo '--version pre-release' )
    END_VERSIONS
    """
}


process BIN_TRANSFORM {
    
    publishDir "${params.output_dir}/transformed_bins", mode:'copy'
    
    
    errorStrategy { task.attempt <= 3 ? "retry" : "finish" }
    maxRetries 3
    
    input:
    tuple val(meta), val(meta2), path(strainy_final_gfa)

    output:
    tuple val(meta), val(meta2), path("${meta2}.fa"),  		   emit: transform_ch
    path("${meta2}.fa"),					   emit: path_transform_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    awk '//^S{{print ">"\$2"\\n"\$3}}' $strainy_final_gfa  > ${meta2}.fa 
        
    """
}


process FILTER_BY_MAG_QUALITY {
    
    publishDir "${params.output_dir}/", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
    input:
    path(quality_report)
    path(bam_directory)

    output:
    path("best_quality/*.bam"),  		   emit: bam_mag_ch
    path("best_quality/*.bai"),  		   emit: bai_mag_ch
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    filter_best.py -i $quality_report -o best_quality -b ${bam_directory}  
        
    """
}


process COUNT_MAG_COVERAGE {
    
    publishDir "${params.output_dir}/transformed_bins", mode:'copy'
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
  
    input:
    path(collected_bam)
    path(bam_directory)
           
    output:
    path("coverage.lst"),  		   emit: count_mag_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    calc_cov_from_bam_dir.py -i $bam_directory -o coverage.lst  
        
    """
}


process FILTER_BY_COVERAGE {
    
    publishDir "${params.output_dir}/transformed_bins", mode:'copy'
    
    
    errorStrategy { task.attempt <= 2 ? "retry" : "finish" }
    maxRetries 2
    
        
    input:
    path(coverage_list)
    path(transformed_bins_directory)
    path(collected_bin_transform)

    output:
    path("passed_coverage_threshold/*.fa"),  		   emit: fa_mag_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    cov_filtered.py -i $coverage_list -o passed_coverage_threshold -d $transformed_bins_directory  
        
    """
}


process MAG_TRANSFORMED_QA {
    
    cache 'lenient'
    
    publishDir "${params.output_dir}", mode:'copy'
    
    
    errorStrategy { task.attempt <= 3 ? "retry" : "finish" }
    maxRetries 3
    
    input:
    path(collected_bin_fa)

    output:
    path("qa_transformed_bins/diamond_output"),        emit: diamond_dir_ch
    path("qa_transformed_bins/protein_files"),         emit: protein_dir_ch
    path("qa_transformed_bins/quality_report.tsv"),    emit: quality_report_ch
    path("qa_transformed_bins/checkm2.log"),           emit: log_ch
    path("versions.yml"),                              emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    checkm2.py predict -i $collected_bin_fa -o qa_transformed_bins --database_path ${baseDir}/CheckM2_database/uniref100.KO.1.dmnd -t $task.cpus  
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         checkm2: \$( checkm2.py 2>&1 | grep '...:::' | sed 's/.*CheckM2 v//;s/ .*//' )
    END_VERSIONS
    """
}

process CUSTOM_DUMPSOFTWAREVERSIONS {
    
    publishDir "${params.output_dir}", mode:'copy'
    
    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml_ch
    path "software_versions_mqc.yml", emit: mqc_yml_ch
    path "versions.yml"             , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
