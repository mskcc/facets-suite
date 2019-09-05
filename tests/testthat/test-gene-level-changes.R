
test_that("gene-level function works", {
  
    test_gene_level = gene_level_changes(test_facets_output, genome = 'hg38')
    
    # check output is data frame
    expect_is(test_gene_level, 'data.frame')
    
    # check output columns
    
    test_names = c('gene', 'chrom', 'gene_start', 'gene_end', 'seg', 'median_cnlr_seg', 'segclust', 'seg_start',
                   'seg_end', 'cf', 'tcn', 'lcn', 'seg_length', 'mcn', 'genes_on_seg', 'gene_snps', 'gene_het_snps',
                   'spans_segs', 'cn_state')
    expect_true(all(test_names %in% names(test_gene_level)))
})
