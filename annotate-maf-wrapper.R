#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(parallel)
    library(facetsSuite)
})

args = commandArgs(TRUE)
if (length(args) == 0) {
    message('Run annotate-maf-wrapper.R --help for list of input arguments.')
    quit()
}

parser = ArgumentParser(description = 'Annotate MAF with local copy number and CCF estimates.')
parser$add_argument('-m', '--maf-file', required = T,
                    help = 'Mutations in MAF format')
parser$add_argument('-s', '--sample-mapping', required = T,
                    help = 'Tab-separated file with header where first column contains sample ID and second column path to FACETS .rds output files.')
parser$add_argument('-a', '--facets-algorithm', required = F,
                    default = 'em', help = 'Which FACETS algorithm to use [default %(default)s]')
parser$add_argument('-o', '--output', required = F,
                    help = 'Output directory [default input.ccf.maf')
parser$add_argument('-p', '--parallel', required = F,
                    default = FALSE,  help = 'Parallelize [default FALSE]')
args = parser$parse_args()

output = args$output %||% paste0(gsub('\\.[a-z]+$', '', args$maf_file), '.ccf.maf')

# Match sample IDs in input ---------------------------------------------------------------------------------------
sample_map = fread(args$sample_mapping, header = TRUE, col.names = c('sample', 'file'))
maf = fread(args$maf_file, header = T)

if (any(duplicated(sample_map$sample))) stop('Sample map contains duplicate sample names', call. = FALSE)

maf_unique_samples = length(setdiff(unique(maf$Tumor_Sample_Barcode), sample_map$sample))
sample_map_unique_samples = length(setdiff(sample_map$sample, unique(maf$Tumor_Sample_Barcode)))
common_samples = intersect(unique(maf$Tumor_Sample_Barcode), sample_map$sample)
if (maf_unique_samples > 0) message(paste(maf_unique_samples, 'samples in input MAF not in sample map.'))
if (sample_map_unique_samples > 0) message(paste(sample_map_unique_samples, 'samples in sample map not in input MAF.'))

# Annotate --------------------------------------------------------------------------------------------------------
message(paste('Annotating', length(common_samples), 'samples'))

annotate_sample = function(sample_id) {
    sample_maf = maf[maf$Tumor_Sample_Barcode == sample_id,]
    sample_facets = readRDS(sample_map$file[which(sample_map$sample == sample_id)])
    annotate_maf(sample_maf, sample_facets$segs, sample_facets$purity)
}

if (args$parallel == TRUE) {
    output_maf = mclapply(common_samples, annotate_sample, mc.cores = detectCores())
} else {
    output_maf = lapply(common_samples, annotate_sample)
}
output_maf = rbindlist(output_maf)

# Add back samples that were missing in sample map
output_maf = rbind.fill(output_maf,
                        maf[!which(maf$Tumor_Sample_Barcode %in% sample_map$sample), ])

write.table(output_maf, output, quote = F, sep = '\t', col.names = T, row.names = F)
 