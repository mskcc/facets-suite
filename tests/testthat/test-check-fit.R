
test_check_fit = check_fit(test_facets_output, maf = test_maf)

test_that('function works', {
    
    expect_is(test_check_fit, 'list')
    
    test_names = c('diplogr_flag', 'n_alternative_diplogr', 'n_homdel_muts', 'median_vaf_homdel_muts')
    expect_true(all(test_names %in% names(test_check_fit)))
})
