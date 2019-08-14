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

parser$add_argument('-v', '--verbose', action="store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-vcf', '--vcf-file', required = T,
                    help = 'Path to VCF file containing SNP positions')
parser$add_argument('-n', '--normal-bam', required = T,
                    help = 'Path to normal sample BAM file')
parser$add_argument('-t', '--tumor-bam', required = T,
                    help = 'Path to tumor sample BAM file')
parser$add_argument('-nn', '--normal-name', required = F,
                    help = 'Name of normal sample')
parser$add_argument('-tn', '--tumor-name', required = F,
                    help = 'Name of tumor sample')
parser$add_argument('-o', '--output-file', required = F,
                    help = 'Name of output file [default countsMerged_tumor_normal.gz]')
parser$add_argument('-p', '--pseudo-snps', required = F, default = 50,
                    help = 'Do pileup at every p:th position [default %(default)s]')
parser$add_argument('-d', '--max-depth', required = F, default = 4000,
                    help = 'Maximum read depth [default %(default)s]')

args = parser$parse_args()

# Prepare output --------------------------------------------------------------------------------------------------
snp_pileup_path = '$HOME/git/facets/inst/extcode/snp-pileup'

if (is.null(args$normal_name)) args$normal_name = gsub('.bam$', '', basename(args$normal_bam))
if (is.null(args$tumor_name)) args$tumor_name = gsub('.bam$', '', basename(args$tumor_bam))
if (is.null(args$output_file)) args$output_file = paste0(args$tumor_name, '_', args$normal_name, 'snp_pileup.gz')
if (file.exists(args$output_file)) stop(paste(args$output_file, 'already exists. Remove before running.'), call. = F)

default_args = c('--count-orphans', '--gzip')

pilup_cmd = paste(
    snp_pileup_path,
    default_args,
    '-P', args$pseudo_snps,
    '-d', args$max_depth,
    args$vcf_file,
    args$output_file,
    args$normal_bam,
    args$tumor_bam
    )

system(pilup_cmd)



