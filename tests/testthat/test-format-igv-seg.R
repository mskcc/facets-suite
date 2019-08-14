
test_output = format_igv_seg(test_facets_output, 'test')

test_that('output is data frame', {
    expect_is(test_output, 'data.frame')
    
    test_igv_columns = c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')
    expect_true(all(test_igv_columns %in% names(test_output)))
})

