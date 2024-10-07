/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/
//

include { BWA_ALIGN_READS } from "../../prism-sdgmc-modules/subworkflows/align_bwa"
include { GATK_BEST_PRACTICES } from "../../prism-sdgmc-modules/subworkflows/gatk_best_practices"
include { VCF_FILTER_AND_DECOMPOSE } from "../../prism-sdgmc-modules/subworkflows/vcf_filter_and_decompose"
include { VEP_ANNOTATE } from "../../prism-sdgmc-modules/subworkflows/vep_annotation"
include { AUTOSOLVE_MULTISAMPLE } from "../../prism-sdgmc-modules/subworkflows/autosolve/autosolve_multisample"
include { BAM_QC } from "../../prism-sdgmc-modules/subworkflows/bam_qc"
include { EXOMEDEPTH_CNV_CALLING } from "../../prism-sdgmc-modules/subworkflows/exomedepth"
include { SVAFOTATE } from "../../prism-sdgmc-modules/subworkflows/svafotate"
include { EXOMEDEPTH_POSTPROCESS } from "../../prism-sdgmc-modules/subworkflows/exomedepth_postprocess"
include { GSEAPY } from "../../prism-sdgmc-modules/subworkflows/gseapy"
include { SMACA } from "../../prism-sdgmc-modules/subworkflows/smaca"
include { MITOCALLER_ANALYSIS } from "../../prism-sdgmc-modules/subworkflows/mitocaller"
include { SOMALIER } from "../../prism-sdgmc-modules/subworkflows/somalier"
include { CHECK_FILE_VALIDITY } from "../../prism-sdgmc-modules/subworkflows/file_check"
include { GENERATE_REPORT } from "../../prism-sdgmc-modules/subworkflows/generate_report"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow TARGETED_ANALYSIS {
    
    take:
    reads
    ref_genome
    ref_genome_index
    known_snps_dbsnp
    known_indels
    known_snps_dbsnp_index
    known_indels_index
    target_bed
    target_bed_covered
    vep_cache
    vep_plugins
    vcf_to_tsv_script
    mane_transcript
    autosolve_script
    panel_monoallelic
    panel_biallelic
    clingen
    mutation_spectrum
    refgene_track
    exomedepth_controls
    exomedepth_target_bed
    exomedepth_gene_bed
    chr_list
    convert_tsv_to_vcf_script_for_exomedepth
    svafotate_bed
    exomedepth_annotate_counts_script
    exomedepth_deletion_db
    exomedepth_duplication_db
    add_svaf_script
    process_script_single
    process_script_family
    panel
    decipher
    gene_sets
    gseapy_enrich_script
    header
    mitocaller_result_filter_script
    mitomap
    mitotip
    mitimpact
    modify_versions_log_script
    parameters_file
    check_file_status_script
    tabulate_samples_quality_script
    check_sample_stats_script
    rmd_template
    verifybamid_resources
    pedfile
    somalier_sites
    somalier_onekg_files
    somalier_prism_files

    ch_versions


    main:

    ch_versions = Channel.empty()
    
    BWA_ALIGN_READS(
        reads,
        ref_genome,
        ref_genome_index
    )
    ch_versions = ch_versions.mix(BWA_ALIGN_READS.out.versions)

    ch_aligned_bam = BWA_ALIGN_READS.out.aligned_bam
    GATK_BEST_PRACTICES(
        ch_aligned_bam,
        ref_genome,
        ref_genome_index,
        known_snps_dbsnp,
        known_indels,
        known_snps_dbsnp_index,
        known_indels_index,
        target_bed
    )
    ch_versions = ch_versions.mix(GATK_BEST_PRACTICES.out.versions)

    ch_raw_vcf = GATK_BEST_PRACTICES.out.raw_vcf
    VCF_FILTER_AND_DECOMPOSE(
        ch_raw_vcf,
        ref_genome,
        ref_genome_index
    )

    ch_versions = ch_versions.mix(VCF_FILTER_AND_DECOMPOSE.out.versions)

    ch_decom_norm_vcf = VCF_FILTER_AND_DECOMPOSE.out.decom_norm_vcf
    VEP_ANNOTATE(
        ch_decom_norm_vcf,
        vep_cache,
        vep_plugins,
        vcf_to_tsv_script,
        mane_transcript
    )
    ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions)

    ch_vep_tsv_filtered = VEP_ANNOTATE.out.vep_tsv_filtered.groupTuple()
    AUTOSOLVE_MULTISAMPLE(
        ch_vep_tsv_filtered,
        autosolve_script,
        clingen,
        panel_monoallelic,
        panel_biallelic,
        mutation_spectrum
    )

    ch_bqsr_bam = GATK_BEST_PRACTICES.out.bqsr_bam
    BAM_QC(
        ch_bqsr_bam,
        target_bed_covered,
        ref_genome,
        ref_genome_index,
        refgene_track,
        verifybamid_resources
    )
    ch_versions = ch_versions.mix(BAM_QC.out.versions)

    ch_bqsr_bam_collect = GATK_BEST_PRACTICES.out.bqsr_bam.collect()
    EXOMEDEPTH_CNV_CALLING(
        ch_bqsr_bam_collect,
        exomedepth_controls,
        ref_genome,
        ref_genome_index,
        exomedepth_target_bed,
        exomedepth_gene_bed,
        chr_list
    )
    ch_versions = ch_versions.mix(EXOMEDEPTH_CNV_CALLING.out.versions)

    ch_merged_tsv_flatten = EXOMEDEPTH_CNV_CALLING.out.exomedepth_merged_tsv.flatten()
    SVAFOTATE(
        ch_merged_tsv_flatten,
        convert_tsv_to_vcf_script_for_exomedepth,
        svafotate_bed
    )

    if(params.genotyping_mode == 'single'){
        ch_for_exomedepth_postprocess = VEP_ANNOTATE.out.vep_tsv_filtered
    }else if(params.genotyping_mode == 'joint'){
        ch_for_exomedepth_postprocess = VEP_ANNOTATE.out.vep_tsv_filtered_without_samplename
    }else if(params.genotyping_mode == 'family'){
        ch_for_exomedepth_postprocess = VEP_ANNOTATE.out.vep_tsv_filtered
    }

    ch_merged_tsv = EXOMEDEPTH_CNV_CALLING.out.exomedepth_merged_tsv
    EXOMEDEPTH_POSTPROCESS(
        ch_merged_tsv,
        SVAFOTATE.out.svafotate_vcf,
        exomedepth_annotate_counts_script,
        exomedepth_deletion_db,
        exomedepth_duplication_db,
        add_svaf_script,
        ch_for_exomedepth_postprocess,
        process_script_single,
        panel,
        clingen,
        mutation_spectrum,
        decipher,
        process_script_family,
        pedfile
    )

    ch_for_exomedepth_postprocess.join(EXOMEDEPTH_POSTPROCESS.exomedepth_ch).view()

    //ch_merged_filtered_tsv_for_gseapy = EXOMEDEPTH_POSTPROCESS.out.exomedepth_del_tsv_forgseapy.join(EXOMEDEPTH_POSTPROCESS.out.exomedepth_dup_tsv_forgseapy)
    ch_merged_filtered_del_tsv_for_gseapy = EXOMEDEPTH_POSTPROCESS.out.exomedepth_del_tsv_forgseapy
    ch_merged_filtered_dup_tsv_for_gseapy = EXOMEDEPTH_POSTPROCESS.out.exomedepth_dup_tsv_forgseapy
    GSEAPY(
        ch_merged_filtered_del_tsv_for_gseapy,
        ch_merged_filtered_dup_tsv_for_gseapy,
        gene_sets,
        gseapy_enrich_script)
    ch_versions = ch_versions.mix(GSEAPY.out.versions)

    ch_bqsr_bam_smaca = GATK_BEST_PRACTICES.out.bqsr_bam
    SMACA(
        ch_bqsr_bam_smaca,
        ref_genome,
        ref_genome_index
    )
    ch_versions = ch_versions.mix(SMACA.out.versions)

    ch_bqsr_bam_mito = GATK_BEST_PRACTICES.out.bqsr_bam
    MITOCALLER_ANALYSIS(
        ch_bqsr_bam_mito,
        ref_genome,
        header,
        mitocaller_result_filter_script,
        mitomap,
        mitotip,
        mitimpact
    )

    ch_bqsr_bam_somalier = GATK_BEST_PRACTICES.out.bqsr_bam
    SOMALIER(
        ch_bqsr_bam_somalier,
        ref_genome,
        ref_genome_index,
        somalier_sites,
        params.proband,
        pedfile,
        somalier_onekg_files,
        somalier_prism_files
    )

    tool_versions_ch = ch_versions.collectFile(name: 'versions.log', newLine: true, sort: false)
/*
    //CHECK_FILE_VALIDITY(tool_versions_ch, modify_versions_log_script, parameters_file, BAM_QC.out.depth_of_coverage_stats, VEP_ANNOTATE.out.vep_tsv_filtered, VCF_FILTER_AND_DECOMPOSE.out.decom_norm_vcf, check_file_status_script, tabulate_samples_quality_script, check_sample_stats_script)
    if(params.genotyping_mode == 'single'){
        ch_files_for_single_sample_check = BAM_QC.out.depth_of_coverage_stats.join(VEP_ANNOTATE.out.vep_tsv_filtered).join(VCF_FILTER_AND_DECOMPOSE.out.decom_norm_vcf).join(BAM_QC.out.edited_qualimap_output)
        ch_for_filecheck_processed = ch_files_for_single_sample_check.map { tuple ->
                                                def sampleName = tuple[0]
                                                def allFiles = tuple[1..-1].collectMany { it instanceof List ? it : [it] }
                                                [sampleName, allFiles]
                                     }
      }

    if(params.genotyping_mode == 'joint'){
        ch_for_filecheck_processed = Channel.empty()
    }

    CHECK_FILE_VALIDITY(
            tool_versions_ch, 
            modify_versions_log_script, 
            parameters_file, 
            ch_for_filecheck_processed, 
            check_file_status_script, 
            tabulate_samples_quality_script, 
            check_sample_stats_script,
            BAM_QC.out.depth_of_coverage_stats.flatten().collect(), 
            VEP_ANNOTATE.out.vep_tsv_filtered, 
            VCF_FILTER_AND_DECOMPOSE.out.decom_norm_vcf,
            //BAM_QC.out.verifybam_id_output.flatten().collect(),
            BAM_QC.out.edited_qualimap_output.collect()
    )

    if(params.genotyping_mode == 'single'){
        //ch_for_rmarkdown_single_sample = CHECK_FILE_VALIDITY.out.check_file_validity_wes_singlesample_output.join(BAM_QC.out.depth_of_coverage_stats).combine(CHECK_FILE_VALIDITY.out.version_txt)
        ch_for_rmarkdown_single_sample = CHECK_FILE_VALIDITY.out.check_file_validity_wes_output.join(BAM_QC.out.depth_of_coverage_stats).combine(CHECK_FILE_VALIDITY.out.version_txt)
        ch_for_rmarkdown_processed = ch_for_rmarkdown_single_sample.map { tuple ->
                                                  def sampleName = tuple[0]
                                                  def allFiles = tuple[1..-1].collectMany { it instanceof List ? it : [it] }
                                                  [sampleName, allFiles]
                                    }
    }

    if(params.genotyping_mode == 'joint'){
        ch_for_rmarkdown_processed = Channel.empty()
    }
    
    
    
    GENERATE_REPORT(
        ch_for_rmarkdown_processed,
        rmd_template,
        CHECK_FILE_VALIDITY.out.params_log,
        panel,
        CHECK_FILE_VALIDITY.out.version_txt,
        BAM_QC.out.depth_of_coverage_stats.flatten().collect(),
        CHECK_FILE_VALIDITY.out.check_file_validity_wes_output
    )
*/
    emit:
        GATK_BEST_PRACTICES.out.bqsr_recal_table
        GATK_BEST_PRACTICES.out.bqsr_bam
        GATK_BEST_PRACTICES.out.gvcf_file
        GATK_BEST_PRACTICES.out.gvcf_index
        GATK_BEST_PRACTICES.out.raw_vcf
        VCF_FILTER_AND_DECOMPOSE.out.filtered_vcfs
        VCF_FILTER_AND_DECOMPOSE.out.decom_norm_vcf
        VEP_ANNOTATE.out.annotated_vcf
        VEP_ANNOTATE.out.vep_tsv
        VEP_ANNOTATE.out.vep_tsv_filtered
        VEP_ANNOTATE.out.vep_tsv_filtered_highqual
        AUTOSOLVE_MULTISAMPLE.out.autosolve_tsv
        BAM_QC.out.qualimap_stats
        BAM_QC.out.depth_of_coverage_stats
        BAM_QC.out.verifybam_id_output
        EXOMEDEPTH_CNV_CALLING.out.exomedepth_tsv
        EXOMEDEPTH_CNV_CALLING.out.exomedepth_png
        EXOMEDEPTH_CNV_CALLING.out.exomedepth_rds
        EXOMEDEPTH_CNV_CALLING.out.exomedepth_merged_tsv
        SVAFOTATE.out.svafotate_vcf
        EXOMEDEPTH_POSTPROCESS.out.exomedepth_merged_filtered_tsv
        EXOMEDEPTH_POSTPROCESS.out.postprocess_result
        EXOMEDEPTH_POSTPROCESS.out.exomedepth_del_tsv_forgseapy
        EXOMEDEPTH_POSTPROCESS.out.exomedepth_dup_tsv_forgseapy
        GSEAPY.out.gseapy_output_del_tsv
        GSEAPY.out.gseapy_output_dup_tsv
        SMACA.out.smaca_tsv
        MITOCALLER_ANALYSIS.out.mitocaller_output_summary
        MITOCALLER_ANALYSIS.out.mitocaller_candidate_variants
        MITOCALLER_ANALYSIS.out.mitocaller_filtered_output
        SOMALIER.out.somalier_ancestry_output
        SOMALIER.out.somalier_relate_output
//        CHECK_FILE_VALIDITY.out.version_txt
//        CHECK_FILE_VALIDITY.out.params_log
//        CHECK_FILE_VALIDITY.out.check_file_validity_wes_output
//        GENERATE_REPORT.out.sample_report

        versions = ch_versions
}
