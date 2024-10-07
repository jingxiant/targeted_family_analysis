#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*params.reads   = "*_R{1,2}.{fq,fastq}.gz"
params.ref     = "/prism_data5/share/GATK_Bundle/hg38/newref/Homo_sapiens_assembly38.fasta"
params.ref_fai = "/prism_data5/share/GATK_Bundle/hg38/newref/Homo_sapiens_assembly38.fasta.fai"
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ref_fa = file(params.ref)
ref_fai = file(params.ref_fai)
known_snps_dbsnp = file(params.known_snps_dbsnp)
known_snps_dbsnp_index = file(params.known_snps_dbsnp + '.tbi')
known_indels_dbsnp = file(params.known_indels)
known_indels_dbsnp_index = file(params.known_indels + '.tbi')
target_bed = file(params.target_bed)
target_bed_covered = file(params.target_bed_covered)
vep_cache = file(params.vep_cache_dir)
vep_plugins = file(params.vep_plugin_dir)
vcf_to_tsv_script=file(params.vcf_to_tsv)
mane_transcript=file(params.mane_transcript)
clingen = file(params.clingen)
mutation_spectrum = file(params.mutation_spectrum)
decipher = file(params.decipher)
autosolve_script = file(params.autosolve_script)
panel_monoallelic = file(params.panel_monoallelic)
panel_biallelic = file(params.panel_biallelic)
refseq_gene_track = file(params.refseq_gene_track)
exomedepth_control = file(params.exomedepth_control)
exomedepth_target_bed = file(params.exomedepth_target_bed)
exomedepth_gene_bed = file(params.exomedepth_gene_bed)
chr = params.chr?.tokenize(',')
svafotate_bed = file(params.svafotate_bed)
add_svaf_script=file(params.add_svaf_script)
convert_tsv_to_vcf_script = file(params.convert_tsv_to_vcf_script)
exomedepth_annotate_counts_script = file(params.exomedepth_annotate_counts_script)
exomedepth_deletion_db = file(params.exomedepth_deletion_db)
exomedepth_duplication_db = file(params.exomedepth_duplication_db)
process_script_single = file(params.process_script_single)
process_script_family = file(params.process_script_family)
panel = file(params.panel)
gene_sets=file(params.gene_sets)
gseapy_enrich_script=file(params.gseapy_enrich_script)
header = file(params.header)
mitomap = file(params.mitomap)
mitotip = file(params.mitotip)
mitimpact = file(params.mitimpact)
mitocaller_result_filter_script = file(params.mitocaller_result_filter)
software_version_modify_script = file(params.software_version_modify_script)
params_file = file(params.params_file)
check_file_status_script = file(params.check_file_status_script)
tabulate_samples_quality_script = file(params.tabulate_samples_quality_script)
check_sample_stats_script = file(params.check_sample_stats_script)
rmd_template = file(params.report_template)
verifybamid_resources = file(params.verifybamid_resources_wes)
pedfile = file(params.pedfile)
somalier_sites = file(params.somalier_sites)
somalier_onekg_files = file(params.somalier_onekg_files)
somalier_prism_files = file(params.somalier_prism_files)
slivar_gnomadpath = file(params.slivar_gnomadpath)
slivar_jspath = file(params.slivar_jspath)
gff3 = file(params.gff_path)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get sample id
//
def getLibraryId( file ) {
    def filename = file.split(/\//)[-1]  // Extract filename from the full path
    def pattern = params.pattern  // Pattern to split the filename
    def join_number = params.join_number.toString().split(',').collect { it as int }
    // Split the filename by the given pattern
    def parts = filename.split(pattern)

    // Concatenate parts based on join_number indices
    def libraryId = join_number.collect { parts[it] }.join('_') 

    return libraryId
}

//def getLibraryId( file ) {
//        file.split(/\//)[-1].split(/_/)[0]
//}

include {TARGETED_ANALYSIS} from "./workflows/targeted_family_analysis"

workflow PRISM_TARGETED_ANALYSIS {

    main:

    ch_versions = Channel.empty()

    input_files = "${params.input}/${params.fastq_file_pattern}"

    Channel
        .fromFilePairs( input_files , flat: true )
//        .fromFilePairs( params.reads, flat: true )
        .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
        .groupTuple()
        .set {reads}

    TARGETED_ANALYSIS(
        reads, 
        ref_fa, 
        ref_fai,
        known_snps_dbsnp,
        known_indels_dbsnp,
        known_snps_dbsnp_index,
        known_indels_dbsnp_index,
        target_bed,
        target_bed_covered,
        vep_cache,
        vep_plugins,
        vcf_to_tsv_script,
        mane_transcript,
        autosolve_script,
        panel_monoallelic,
        panel_biallelic,
        clingen,
        mutation_spectrum,
        refseq_gene_track,
        exomedepth_control,
        exomedepth_target_bed,
        exomedepth_gene_bed,
        chr,
        convert_tsv_to_vcf_script,
        svafotate_bed,
        exomedepth_annotate_counts_script,
        exomedepth_deletion_db,
        exomedepth_duplication_db,
        add_svaf_script,
        process_script_single,
        process_script_family,
        panel,
        decipher,
        gene_sets,
        gseapy_enrich_script,
        header,
        mitocaller_result_filter_script,
        mitomap,
        mitotip,
        mitimpact,
        software_version_modify_script,
        params_file,
        check_file_status_script,
        tabulate_samples_quality_script,
        check_sample_stats_script,
        rmd_template,
        verifybamid_resources,
        pedfile,
        somalier_sites,
        somalier_onekg_files,
        somalier_prism_files,
        gff3_file,
        slivar_gnomadpath,
        slivar_jspath,
        
        ch_versions
    )
    ch_versions = ch_versions.mix(TARGETED_ANALYSIS.out.versions)
}

workflow {

    main: 
    PRISM_TARGETED_ANALYSIS()
}



    
