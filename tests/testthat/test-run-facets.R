
library(facetsSuite)
library(pctGCdata)

test_that('function returns output', {
    
    test_facets_run = run_facets(test_read_counts)
    
    test_output_names = c('snps', 'segs', 'purity', 'ploidy', 'diplogr')
    expect_is(test_facets_run, 'list')
    expect_true(all(test_output_names %in% names(test_facets_run)))
})
