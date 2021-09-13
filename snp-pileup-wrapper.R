#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(argparse)
})

args = commandArgs(TRUE)
if (length(args) == 0) {
    message('Run snp-pileup-wrapper.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Generate SNP read counts from matched tumor-normal BAM files.')

parser$add_argument('-v', '--verbose', action = "store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-co', '--count-orphans', help= "count anamolous read pairs", required = FALSE)
parser$add_argument('-sp', '--snp-pileup-path', required = FALSE,
                    help = 'Path to snp-pileup executable [default environment variable $SNP_PILEUP]')
parser$add_argument('-vcf', '--vcf-file', required = TRUE,
                    help = 'Path to VCF file containing SNP positions')
parser$add_argument('-n', '--normal-bam', required = FALSE,
                    help = 'Path to normal sample BAM file')
parser$add_argument('-t', '--tumor-bam', required = FALSE,
                    help = 'Path to tumor sample BAM file')
parser$add_argument('-o', '--output-prefix', required = TRUE,
                    help = 'Name prefix for output file')
parser$add_argument('-p', '--pseudo-snps', required = FALSE, default = 50,
                    help = 'Do pileup at every p:th position [default %(default)s]')
parser$add_argument('-d', '--max-depth', required = FALSE, default = 20000,
                    help = 'Maximum read depth [default %(default)s]')
parser$add_argument('-q', '--min-map-quality', required = FALSE, default = 0,
                    help = 'Sets the minimum threshold for mapping quality')
parser$add_argument('-Q', '--min-base-quality', required = FALSE, default = 0,
                    help = 'Sets the minimum threshold for base quality.')
parser$add_argument('-un', '--unmatched-normal-BAMS', required = FALSE, default = FALSE,
                    help = 'full path(s) as quoted string to unmatched normal BAM(s) to use for log ratio normalization, e.g. "<some/path>/*-NS_*bam"')

args = parser$parse_args()

# Prepare output --------------------------------------------------------------------------------------------------

snp_pileup_env = Sys.which("snp-pileup")

if (!is.null(args$snp_pileup_path)) {
    if (snp_pileup_env == '') {
        stop(paste('No snp-pileup path provided or in user environment.'), call. = F)
    } else {
        snp_pileup_path = snp_pileup_env
    }
}

output_file = paste0(args$output_prefix, '.snp_pileup.gz')

if (file.exists(output_file)) {
    stop(paste(output_file, 'already exists. Remove before running.'), call. = F)
}
default_args = c('--gzip')

#enforce minmum read count of 10 for all normals analyzed
if (args$unmatched_normal_BAMS!=FALSE){
  unmatched_bam_count = length(system(paste("ls ", args$unmatched_normal_BAMS), intern=TRUE))
  if (is.null(args$normal_bam) & is.null(args$tumor_bam)) {
    message("generating counts file for ", unmatched_bam_count, " BAMs")
    min_read_counts =gsub(", ", ",", toString(c(" --min-read-counts 1", rep("1", unmatched_bam_count-1))))
  }else{
    min_read_counts = gsub(", ", ",", toString(c(" --min-read-counts 10,0", rep("0", unmatched_bam_count))))
    message("incorporating ", unmatched_bam_count, " unmatched bams into analysis")
  }
}else{
  min_read_counts = " --min-read-counts 10,0"
}

pileup_cmd = paste(
    snp_pileup_path,
    default_args,
    '-P', args$pseudo_snps,
    '-d', args$max_depth,
    '-q', args$min_map_quality,
    '-Q', args$min_base_quality,
    args$count_orphans,
    min_read_counts,
    args$vcf_file,
    output_file,
    args$normal_bam,
    args$tumor_bam,
    args$unmatched_normal_BAMS

)

message("pileup cmd: ", pileup_cmd)
system(pileup_cmd)
