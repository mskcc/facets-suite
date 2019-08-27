

test_that('function returns output', {
    
    test_annotation = arm_level_changes(test_facets_output$segs, test_facets_out$ploidy)
    
    # check output contains all items
    test_names = c('genome_doubled', 'fraction_cna', 'weighted_fraction_cna', 'aneuploidy_score', 'full_output')
    expect_true(all(test_names %in% names(test_annotation)))
})
