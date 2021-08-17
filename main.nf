#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/traits
========================================================================================
 lifebit-ai/traits Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/lifebit-ai/traits
----------------------------------------------------------------------------------------
*/

/*---------------------------------------
  Define and show help message if needed
-----------------------------------------*/

def helpMessage() {

    log.info"""
    
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --post_analysis heritability --input_gwas_statistics 
    Essential parameters:
    --post_analysis                  Type of analysis desired. Options are 'heritability' or 'genetic_correlation_h2'
    --input_gwas_statistics                     path to summary statistics from gwas. Currently compatible with SAIGE
    
    Optional parameters:
    --external_gwas_statistics                   Path to second summary statistics to be used for genetic correlation                     
    --external_gwas_cat_study_id              String containing GWAS catalogue study id
    --external_gwas_cat_study_size            Integer describing size of GWAS study
    --external_gwas_cat_ftp             Path to csv containing ftp locations of gwas catalogue files
    --hapmap3_snplist                Path to SNP list from Hapmap needed for seleting SNPs considered for analysis
    --ld_scores_tar_bz2              Path to tar.bz2 files with precomputed LD scores
    --outdir                         Path to output directory
    --output_tag                     String containing output tag
    """.stripIndent()
}

// Show help message

if (params.help) {
    helpMessage()
    exit 0
}



/*---------------------------------------------------
  Define and show header with all params information 
-----------------------------------------------------*/

// Header log info

def summary = [:]

if (workflow.revision) summary['Pipeline Release'] = workflow.revision

summary['Max Resources']                  = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output dir']                     = params.outdir
summary['Launch dir']                     = workflow.launchDir
summary['Working dir']                    = workflow.workDir
summary['Script dir']                     = workflow.projectDir
summary['User']                           = workflow.userName

summary['post_analysis']                  = params.post_analysis
summary['input_gwas_statistics']                     = params.input_gwas_statistics
summary['external_gwas_statistics']                   = params.external_gwas_statistics
summary['external_gwas_cat_study_id']              = params.external_gwas_cat_study_id
summary['external_gwas_cat_study_size']            = params.external_gwas_cat_study_size
summary['external_gwas_cat_ftp']             = params.external_gwas_cat_ftp
summary['hapmap3_snplist']                = params.hapmap3_snplist
summary['ld_scores_tar_bz2']              = params.ld_scores_tar_bz2
summary['output_tag']                     = params.output_tag
summary['outdir']                         = params.outdir

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*--------------------------------------------------
  Channel preparation
---------------------------------------------------*/

ch_ldsc_input = params.input_gwas_statistics ? Channel.value(file(params.input_gwas_statistics)) : Channel.empty()
ch_hapmap3_snplist =  params.hapmap3_snplist ? Channel.value(file(params.hapmap3_snplist)) :  "null"
ch_ld_scores_tar_bz2 =  params.ld_scores_tar_bz2 ? Channel.value(file(params.ld_scores_tar_bz2)) :  "null"
ch_gwas_summary = params.external_gwas_statistics ? Channel.value(file(params.external_gwas_statistics)) : Channel.empty()

if (params.post_analysis == 'genetic_correlation_h2'){
  if (!params.external_gwas_statistics && !params.external_gwas_cat_study_id){
    exit 1, "Second summary statistics file, or its study ID from GWAS Catalogue is required for estimating genetic correlation."
  }
}

/*--------------------------------------------------
  LDSC - Genetic correlation and heritability
---------------------------------------------------*/
if (params.post_analysis == 'heritability' || params.post_analysis == 'genetic_correlation_h2'){
  process prepare_files_ldsc {
    tag "preparation_files"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(summary_stats) from ch_ldsc_input

    output:
    file("${params.output_tag}_transformed_gwas_stats.txt") into ch_ldsc_input2

    script:

    """
    convert_output.R \
      --gwas_stats "$summary_stats" \
      --output_tag ${params.output_tag}
    """
  }
  process munge_saige_output {
    tag "munge_saige_output"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(saige_summary_stats) from ch_ldsc_input2
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${params.output_tag}_ldsc.sumstats.gz") into ch_saige_ldsc

    script:

    """
    munge_sumstats.py --sumstats $saige_summary_stats \
                      --out "${params.output_tag}_ldsc" \
                      --merge-alleles $hapmap3_snplist \
                      --a1 Allele1 \
                      --a2 Allele2 \
                      --signed-sumstats BETA,0 \
                      --p p.value \
                      --snp SNPID \
                      --info inputationInfo
    """
  }
}

if (params.post_analysis == 'heritability'){

  process heritability {
    tag "heritability"
    publishDir "${params.outdir}/heritability/", mode: 'copy'

    input:
    file(saige_output) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_h2.log") into ch_ldsc_report_input

    script:
    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
      --h2 $saige_output \
      --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
      --out ${params.output_tag}_h2
    """
  }
}

if (params.post_analysis == 'genetic_correlation_h2' && params.external_gwas_statistics){
  process prepare_gwas_summary_ldsc {
    tag "preparation_gwas_summary_ldsc"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(gwas_summary_file) from ch_gwas_summary

    output:
    file("${gwas_summary_file.simpleName}_transformed_gwas_stats.txt") into ch_gwas_summary_ldsc

    script:

    """
    convert_output.R \
      --gwas_stats "$gwas_summary_file" \
      --output_tag "${gwas_summary_file.simpleName}"
    """
  }
  //* Munge gwas stats

  process munge_gwas_summary {
    tag "munge_gwas_summary"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(summary_stats) from ch_gwas_summary_ldsc
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${summary_stats.simpleName}_gwas_summary.sumstats.gz") into ch_gwas_summary_ldsc2

    script:

    """
    munge_sumstats.py \
          --sumstats "$summary_stats" \
          --out "${summary_stats.simpleName}_gwas_summary" \
          --merge-alleles $hapmap3_snplist
    """
  }

  //* Run genetic correlation
  process genetic_correlation_h2 {
    tag "genetic_correlation_h2"
    publishDir "${params.outdir}/genetic_correlation/", mode: 'copy'

    input:
    file(gwas_summary_ldsc) from ch_gwas_summary_ldsc2
    file(saige_ldsc) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_genetic_correlation.log") into ch_ldsc_report_input

    script:

    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
          --rg $saige_ldsc,$gwas_summary_ldsc \
          --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --out ${params.output_tag}_genetic_correlation \
          --no-intercept
    """
  }

}

 //* gwas catalogue

if (params.post_analysis == 'genetic_correlation_h2' && params.external_gwas_cat_study_id){

  external_gwas_cat_ftp_ch = Channel.fromPath(params.external_gwas_cat_ftp, checkIfExists: true)
    .ifEmpty { exit 1, "GWAS catalogue ftp locations not found: ${params.external_gwas_cat_ftp}" }
    .splitCsv(header: true)
    .map { row -> tuple(row.study_accession, row.ftp_link_harmonised_summary_stats) }
    .filter{ it[0] == params.external_gwas_cat_study_id}
    .ifEmpty { exit 1, "The GWAS study accession number you provided does not come as a harmonized dataset that can be used as a base cohort "}
    .flatten()
    .last()

  process download_gwas_catalogue {
    label "high_memory"
    publishDir "${params.outdir}/GWAS_cat", mode: "copy"
    
    input:
    val(ftp_link) from external_gwas_cat_ftp_ch
    
    output:
    file("*.h.tsv*") into downloaded_gwas_catalogue_ch
    
    script:
    """
    wget ${ftp_link}
    """
  }

  process transform_gwas_catalogue {
    label "high_memory"
    publishDir "${params.outdir}/GWAS_cat", mode: "copy"
    
    input:
    file gwas_catalogue_base from downloaded_gwas_catalogue_ch
    
    output:
    file("${params.external_gwas_cat_study_id}.data") into transformed_base_ch
    
    script:
    """
    transform_gwas_catalogue.R --input_gwas_cat "${gwas_catalogue_base}" \
                               --outprefix "${params.external_gwas_cat_study_id}"
    """
    }

  
  //* Munge gwas cat stats

  process munge_gwas_cat_summary {
    tag "munge_gwas_summary"
    publishDir "${params.outdir}/ldsc_inputs/", mode: 'copy'

    input:
    file(summary_stats) from transformed_base_ch
    file(hapmap3_snplist) from ch_hapmap3_snplist

    output:
    file("${params.external_gwas_cat_study_id}_gwas_summary.sumstats.gz") into ch_gwas_summary_ldsc2

    script:

    """
    munge_sumstats.py \
          --sumstats "$summary_stats" \
          --out "${params.external_gwas_cat_study_id}_gwas_summary" \
          --merge-alleles $hapmap3_snplist \
          --signed-sumstats BETA,0 \
          --N ${params.external_gwas_cat_study_size}
    """
  }

  //* Run genetic correlation
  process genetic_correlation_h2_gwas_cat {
    tag "genetic_correlation_h2"
    publishDir "${params.outdir}/genetic_correlation/", mode: 'copy'

    input:
    file(gwas_summary_ldsc) from ch_gwas_summary_ldsc2
    file(saige_ldsc) from ch_saige_ldsc
    file(ld_scores_tar_bz2) from ch_ld_scores_tar_bz2

    output:
    file("${params.output_tag}_genetic_correlation.log") into ch_ldsc_report_input

    script:
    """
    tar -xvjf ${ld_scores_tar_bz2}

    ldsc.py \
          --rg $saige_ldsc,$gwas_summary_ldsc \
          --ref-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --w-ld-chr ${ld_scores_tar_bz2.simpleName}/ \
          --out ${params.output_tag}_genetic_correlation \
          --no-intercept
    """
  }

}



def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  lifebit-ai/traits v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
