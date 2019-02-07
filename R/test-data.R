#' Test data
#'
#' Data for user and development testing, based on TCGA sample TCGA-06-0154.
#' All data from GDC's data portal, processed as described below.
#'
#' @source \url{portal.gdc.cancer.gov}
#' @name plot_facets

# Output from snp-pileup, downsampled to reduce file size
# VCF: ftp://ftp.ncbi.nih.gov/snp/.redesign/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
# Normal BAM file: TCGA-06-0154-01A-03D-1491-08
# Tumor BAM file: TCGA-06-0154-10A-01D-1491-08
# snp-pileup -A -g -r 10,0 --pseudo-snps=50 dbsnp.vcf normal.bam tumor.bam countsMerged.gz
# test_read_counts = read_snp_matrix('countsMerged.gz') %>% 
#     group_by(Chromsome) %>% 
#     sample_frac(size = .75) %>%
#     arrange(Chromosome, Position)
'test_read_counts'

# Output from run-facets.R
# test_facets_output = run_facets(test_read_counts, cval = 500, genome = 'hg38')
'test_facets_output'

# MAF file downloaded from GDC's data portal
# test_maf = mutate(test_maf, Chromosome = str_replace(Chromosome, 'chr', '')) %>% 
#     filter(GDC_FILTER == '')
'test_maf'

# usethis::use_data(test_read_counts, test_facets_output, test_maf)





