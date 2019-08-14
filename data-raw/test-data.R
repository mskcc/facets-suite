
# Generate test data ----------------------------------------------------------------------------------------------
library(dplyr)
library(facetsSuite)
library(usethis)

# test_read_counts
# Output from snp-pileup, downsampled to reduce file size
# VCF: ftp://ftp.ncbi.nih.gov/snp/.redesign/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
# Normal BAM file: TCGA-06-0154-01A-03D-1491-08
# Tumor BAM file: TCGA-06-0154-10A-01D-1491-08
# snp-pileup -A -g -r 10,0 --pseudo-snps=50 dbsnp.vcf normal.bam tumor.bam countsMerged.gz
# zcat countsMerged.gz | head -100 > countsMerged_mini.gz
test_read_counts = read_snp_matrix('data-raw/countsMerged.gz') %>%
    group_by(Chromosome) %>%
    sample_frac(size = .75) %>%
    arrange(Chromosome, Position)

# Output from run-facets.R
read_counts = read_snp_matrix('data-raw/countsMerged.gz')
test_facets_output = run_facets(read_counts, cval = 500, genome = 'hg38')

# MAF file downloaded from GDC's data portal
test_maf = fread('data-raw/TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf.gz') %>% 
    mutate( Chromosome = str_replace(Chromosome, 'chr', '')) %>%
    filter(GDC_FILTER == '')

use_data(test_read_counts, test_facets_output, test_maf, compress = 'bzip2', overwrite = T)



