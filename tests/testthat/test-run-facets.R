
test_that('function returns output', {
    
    test_facets_run = run_facets(test_read_counts)
    
    # check that function returns a list
    expect_is(test_facets_run, 'list')
    
    # check that the list contains the following items
    test_output_names = c('snps', 'segs', 'purity', 'ploidy', 'diplogr')
    expect_true(all(test_output_names %in% names(test_facets_run)))
    
    # check that error is returned if bad input is passed
    expect_error(run_facets(read_counts = mtcars))
})
