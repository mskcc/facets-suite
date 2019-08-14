
test_annotation = annotate_maf(test_maf, test_facets_output$segs, test_facets_output$purity)

test_that('output is data frame', {
    expect_is(test_annotation, 'data.frame')
})

test_that('output has annotation columns', {
    test_columns = c('tcn', 'lcn', 'ccf_expected_copies', 'ccf_expected_copies_lower', 'ccf_expected_copies_upper')
    expect_true(all(test_columns %in% names(test_annotation)))
})
