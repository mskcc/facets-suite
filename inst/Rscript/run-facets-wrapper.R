#!/usr/bin/env Rscript --vanilla
suppressPackageStartupMessages({
    library(facetsSuite)
    library(argparse)
    library(dplyr)
    library(facets)
    library(ggplot2)
    library(Cairo)
    library(egg)
})

parser = ArgumentParser(description = 'Execute FACETS and generate various output, input SNP read counts from snp-pileup.')

parser$add_argument('-v', '--verbose', action="store_true", default = TRUE,
                    help = 'Print run info')
parser$add_argument('-f', '--counts-file', required = T,
                    help = 'Merged, gzipped tumor-normal output from snp-pileup')
parser$add_argument('-s', '--sample-id', required = F,
                    help = 'Sample ID, preferrable Tumor_Normal to keep track of the normal used')
parser$add_argument('-D', '--directory', required = F,
                    help = 'Output directory [default sample_id')
parser$add_argument('-g', '--genome', required = F,
                    choices = c('hg18', 'hg19', 'hg38'),
                    default = 'hg19', help = 'Reference genome [default %(default)s]')
parser$add_argument('-c', '--cval', required = F, type = 'integer',
                    default = 100, help = 'Segmentation parameter (cval) [default %(default)s]')
parser$add_argument('-pc', '--purity-cval', required = F, type = 'integer',
                    default = NULL, help = 'If two-pass, purity segmentation parameter (cval)')
parser$add_argument('-m', '--min-nhet', required = F,
                    default = 15, help = 'Min. number of heterozygous SNPs required for clustering [default %(default)s]')
parser$add_argument('-pm', '--purity-min-nhet', required = F,
                    default = 15, help = 'If two-pass, purity min. number of heterozygous SNPs (cval) [default %(default)s]')
parser$add_argument('-n', '--snp-window-size', required = F,
                    default = 250, help = 'Window size for heterozygous SNPs [default %(default)s]')
parser$add_argument('-nd', '--normal-depth', required = F,
                    default = 35, help = 'Min. depth in normal to keep SNPs [default %(default)s]')
parser$add_argument('-d', '--diplogr', required = F,
                    default = NULL, help = 'Manual diplogr')
parser$add_argument('-S', '--seed', required = F,
                    default = 100, help = 'Manual seed value [default %(default)s]')
args = parser$parse_args()

# Helper functions ------------------------------------------------------------------------------------------------
# Print run details
print_run_details = function(outfile,
                             cval,
                             min_nhet,
                             purity,
                             ploidy,
                             diplogr,
                             flags = NULL) {
    
    run_details = c(
        '# Input parameters',
        paste0('TAG = ', args$sample_id),
        paste0('Facets version = ', packageVersion('facets')),
        paste0('Input = ', basename(args$counts_file)),
        paste0('Seed = ', args$seed),
        paste0('snp_nbhd = ', args$snp_window_size),
        paste0('ndepth = ', args$normal_depth),
        paste0('cval = ', cval),
        paste0('min.nhet = ', min_nhet),
        paste0('genome = ', args$genome),
        '\n# FACETS output',
        paste0('Purity = ', purity),
        paste0('Ploidy = ', ploidy),
        paste0('dipLogR = ', diplogr),
        paste0('output flags = ', flags)
    )
    
    invisible(
        Map(f = function(x) {
            cat(x, file = outfile, sep = '\n', append = T)
        },
        run_details)
    )
}

# Default set of output plots
print_plots = function(outfile,
                       facets_output,
                       cval) {
    
    plot_title = paste0(args$sample_id,
                        ' | cval=', cval,
                        ' | purity=', round(facets_output$purity, 2),
                        ' | ploidy=', round(facets_output$ploidy, 2),
                        ' | dipLogR=', round(facets_output$diplogr, 2))
    
    CairoPNG(file = outfile, width = 850, height = 999, units = 'px')
    suppressWarnings(
        egg::ggarrange(
            plots = list(
                cnlr_plot(facets_output),
                valor_plot(facets_output),
                icn_plot(facets_output, method = 'em'),
                cf_plot(facets_output, method = 'em'),
                icn_plot(facets_output, method = 'cncf'),
                cf_plot(facets_output, method = 'cncf')
            ),
            ncol = 1,
            nrow = 6,
            heights = c(1, 1, 1, .25, 1, .25),
            top = plot_title)
    )
    dev.off()
}

# Print segmentation
print_segments = function(outfile,
                          facets_output) {
    write.table(facets_output$segs, file = outfile, sep = '\t', quote = F, col.names = T, row.names = F)
}

# Print IGV-style .seg file
print_igv = function(outfile,
                     facets_output) {
    
    ii = format_igv_seg(facets_output = facets_output,
                        sample_id = args$sample_id,
                        normalize = T)
    
    write.table(ii, file = outfile, sep = '\t', quote = F, col.names = T, row.names = F)
}

# Define facets iteration
# Given a set of parameters, do:
# 1. Run facets
# 2. Generate and save plots
# 3. Print run iformation, IGV-style seg file, segmentation data
facets_iteration = function(name_prefix, ...) {
    params = list(...)
    
    output = run_facets(read_counts = read_counts,
                        cval = params$cval,
                        diplogr = params$diplogr,
                        ndepth = params$ndepth,
                        snp_nbhd = params$snp_nbhd,
                        min_nhet = params$min_nhet,
                        genome = params$genome,
                        seed = params$seed)
    
    print_run_details(outfile = paste0(name_prefix, '.out'),
                      cval = params$cval,
                      min_nhet = params$min_nhet,
                      purity = output$purity,
                      ploidy = output$ploidy,
                      diplogr = output$diplogr,
                      flags = output$flags)
    
    print_segments(outfile = paste0(name_prefix, '.cncf.txt'),
                   facets_output = output)
    
    print_igv(outfile = paste0(name_prefix, '.seg'),
              facets_output = output)
    
    print_plots(outfile = paste0(name_prefix, '.png'),
                facets_output = output,
                cval = params$cval)
    
    output
}

# Run -------------------------------------------------------------------------------------------------------------

# Name files and create output directory
if (is.null(args$sample_id)) args$sample_id = gsub('(.dat.gz$|.gz$)', '', basename(args$counts_file))
if (is.null(args$directory)) {
    args$directory = paste0(getwd(), '/', args$sample_id)
} else {
    args$directory = paste0(gsub('[\\/]$', '', args$directory), '/', args$sample_id)
}

if (dir.exists(args$directory)) {
    stop('Output directory already exists, specify different one.')
} else {
    system(paste('mkdir', args$directory))
}

message(paste('Reading', args$counts_file))
message(paste('Writing to', args$directory))

# Read SNP counts file
read_counts = read_snp_matrix(args$counts_file)

# Determine if running two-pass
if (!is.null(args$purity_cval)) {
    purity_prefix = paste0(args$directory, '/', args$sample_id, '_purity')
    purity_output = facets_iteration(name_prefix = purity_prefix, 
                                     sample_id = args$sample_id,
                                     diplogr = args$diplogr,
                                     cval = args$purity_cval,
                                     ndepth = args$normal_depth,
                                     snp_nbhd = args$snp_window_size,
                                     min_nhet = args$purity_min_nhet,
                                     genome = args$genome,
                                     seed = args$seed)
    
    hisens_prefix = paste0(args$directory, '/', args$sample_id, '_hisens')
    hisens_output = facets_iteration(name_prefix = hisens_prefix, 
                                     sample_id = args$sample_id,
                                     diplogr = purity_output$diplogr,
                                     cval = args$cval,
                                     ndepth = args$normal_depth,
                                     snp_nbhd = args$snp_window_size,
                                     min_nhet = args$purity_min_nhet,
                                     genome = args$genome,
                                     seed = args$seed)
    
    saveRDS(purity_output, paste0(args$directory, '/', args$sample_id, '_purity', '.rds'))
    saveRDS(hisens_output, paste0(args$directory, '/', args$sample_id, '_hisens', '.rds'))
    
} else {
    output = facets_iteration(name_prefix = paste0(args$directory, '/', args$sample_id), 
                              sample_id = args$sample_id,
                              cval = args$cval,
                              ndepth = args$normal_depth,
                              snp_nbhd = args$snp_window_size,
                              min_nhet = args$min_nhet,
                              genome = args$genome,
                              seed = args$seed)
    saveRDS(output, paste0(args$directory, '/', args$sample_id, '.rds'))
}
