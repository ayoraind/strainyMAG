#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()

// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)

final_params = check_params(merged_params)

// starting pipeline
pipeline_start_message(version, final_params)



// include processes
include { FLYE_BUILD_ASSEMBLY; STRAINY_SPLIT_UNITIGS; METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS; BINNING; MAG_QA; UNITIG_PER_BIN; SUBSET_BAM; INDEX_BAM; MAKE_SPLIT_FA; INDEX_FA; CALL_SNP; UNZIP_SNP; CONNECT_UNITIGS; SUBSET_GFA; STRAINY; BIN_TRANSFORM; FILTER_BY_MAG_QUALITY; COUNT_MAG_COVERAGE; FILTER_BY_COVERAGE; MAG_TRANSFORMED_QA; CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/processes.nf' addParams(final_params)



workflow  {
          reads_ch = channel
                          .fromPath( final_params.reads )
                          .map { file -> tuple(file.simpleName, file) }
			  .ifEmpty { error "Cannot find any reads matching: ${final_params.reads}" }
			  
	   

         FLYE_BUILD_ASSEMBLY(reads_ch, final_params.valid_mode)
	 

         joined_ch = reads_ch.join(FLYE_BUILD_ASSEMBLY.out.flye_gfa_ch)

 
	 STRAINY_SPLIT_UNITIGS(joined_ch)
	 
	 
	 joined_strainy_bam_ch = STRAINY_SPLIT_UNITIGS.out.long_bam_ch.join(STRAINY_SPLIT_UNITIGS.out.long_bai_ch)
	 
	 
	 METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(joined_strainy_bam_ch)
	 
	 
	 joined_strainy_ch = STRAINY_SPLIT_UNITIGS.out.strainy_gfa_ch.join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth_ch)
	 
	 
	 BINNING(joined_strainy_ch)
	 
	 collected_bins_ch = BINNING.out.bin_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )
	 
	 
	 // I had to convert the unique key in the tuple to strings, remove '_bin' string from the key of the tuple, then remove the bracket [] from the resulting 'list'. Very crude, I know! But I just couldn't find any other way after struggling for hours
	// modified_bin_ch = BINNING.out.bin_ch.map { file -> tuple(file.simpleName.unique().toString().replaceAll('_bin', '')[1..-2], file) }
	 
	
	 flattened_bins_ch = collected_bins_ch
	 		                     .flatten()
	                                     .map { file -> tuple(file.simpleName.toString().replaceAll('_bin', '')[0..-1], file.baseName, file) }
					    
					     
					     
	 

        MAG_QA(collected_bins_ch)
	 
	
	UNITIG_PER_BIN(flattened_bins_ch)
	 
	// joined_unitig_bam_ch = UNITIG_PER_BIN.out.unitig_ch.cross(joined_strainy_bam_ch)
	
	joined_unitig_bam_ch = joined_strainy_bam_ch.combine(UNITIG_PER_BIN.out.unitig_ch, by: 0)
						   
	
	SUBSET_BAM(joined_unitig_bam_ch)
	 
	 
        INDEX_BAM(SUBSET_BAM.out.bam_ch)
	 
	 
        MAKE_SPLIT_FA(STRAINY_SPLIT_UNITIGS.out.long_unitig_ch)
	 
	 
	INDEX_FA(STRAINY_SPLIT_UNITIGS.out.strainy_gfa_ch)
	 
	joined_fa_and_fai_ch = INDEX_FA.out.fai_ch.join(STRAINY_SPLIT_UNITIGS.out.strainy_gfa_ch)
	joined_subset_bam_and_bai_ch = SUBSET_BAM.out.bam_ch.join(INDEX_BAM.out.bam_bai_ch, by: [0,1])
	joined_bam_fa_and_indices_ch = joined_fa_and_fai_ch.combine(joined_subset_bam_and_bai_ch, by: 0)
	
	 
	CALL_SNP(joined_bam_fa_and_indices_ch)
	 
	 
	UNZIP_SNP(CALL_SNP.out.call_snp_ch)
	 
	 
        joined_long_unitig_and_unitig_per_bin_ch = STRAINY_SPLIT_UNITIGS.out.long_unitig_ch.combine(UNITIG_PER_BIN.out.unitig_ch, by: 0)
	
	output_dir_ch = channel.fromPath( final_params.output_dir )
			  
	joined_long_unitig_and_unitig_per_bin_outputdir_ch = joined_long_unitig_and_unitig_per_bin_ch.combine(output_dir_ch)
	
	
	CONNECT_UNITIGS(joined_long_unitig_and_unitig_per_bin_outputdir_ch)
	 
	joined_long_unitig_and_connected_unitig_ch = STRAINY_SPLIT_UNITIGS.out.long_unitig_ch.combine(CONNECT_UNITIGS.out.connected_unitig_ch, by: 0)
	 
	 
	SUBSET_GFA(joined_long_unitig_and_connected_unitig_ch)
	 
	subset_bam_bai_gfa_ch = joined_subset_bam_and_bai_ch.join(SUBSET_GFA.out.subset_gfa_ch, by: [0,1])
	
	subset_bam_bai_gfa_snp_ch = subset_bam_bai_gfa_ch.join(UNZIP_SNP.out.unzip_ch, by: [0,1])
	
	all_join_ch = reads_ch.combine(subset_bam_bai_gfa_snp_ch, by: 0)
	
	 
	// join_read_and_vcf_ch = reads_ch.join(UNZIP_SNP.out.unzip_snp_ch)
	// join_bam_bai_gfa_ch = SUBSET_BAM.out.bam_ch.join(INDEX_BAM.out.bam_bai_ch.join(SUBSET_GFA.out.subset_gfa_ch))
	 
	// all_join_ch = join_read_and_vcf_ch.join(join_bam_bai_gfa_ch)
	 
        STRAINY(all_join_ch)
	
	BIN_TRANSFORM(STRAINY.out.strainy_final_ch)
	
	
       // bam_output_dir_ch = channel.fromPath( final_params.output_dir )
	bams = "bams"
	bam_output_dir_ch = channel.fromPath("${final_params.output_dir}/${bams}")
	
			  
	// quality_report_bam_ch = MAG_QA.out.quality_report_ch.combine(bam_output_dir_ch)
	
	FILTER_BY_MAG_QUALITY(MAG_QA.out.quality_report_ch, bam_output_dir_ch)
	
	// quality_report_bam_to_tuple_ch = FILTER_BY_MAG_QUALITY.out.quality_bam_ch.map { file -> tuple(file.simpleName.toString()[0..-1], file.baseName, file }
	
	
	 COUNT_MAG_COVERAGE(SUBSET_BAM.out.path_bam_ch.collect(), bam_output_dir_ch)
	transformed_bins = "transformed_bins"
	transformed_fa_ch = channel.fromPath ("${final_params.output_dir}/${transformed_bins}")
	// join_coverage_and_transformed_fa_ch = COUNT_MAG_COVERAGE.out.count_mag_ch.combine(transformed_fa_ch)
	
	 FILTER_BY_COVERAGE(COUNT_MAG_COVERAGE.out.count_mag_ch, transformed_fa_ch, BIN_TRANSFORM.out.path_transform_ch.collect())
	
	// MAG_TRANSFORMED_QA(FILTER_BY_COVERAGE.out.fa_mag_ch.collect(  sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} ))
	 
	 // Set up some basic variables
         versions_ch = Channel.empty() 
	 
	 versions_ch = versions_ch
        			.mix(FLYE_BUILD_ASSEMBLY.out.versions_ch)
        			.mix(STRAINY_SPLIT_UNITIGS.out.versions_ch)
        			.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions_ch)
				.mix(BINNING.out.versions_ch)
				.mix(MAG_QA.out.versions_ch)
				.mix(INDEX_BAM.out.versions_ch)
				.mix(CALL_SNP.out.versions_ch)
				

         CUSTOM_DUMPSOFTWAREVERSIONS (
                                versions_ch.unique().collectFile(name: 'collated_versions.yml')
    )
}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
