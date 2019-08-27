
test_that('output is data frame', {
    
    # check that data.frame is returned
    test_output = format_igv_seg(test_facets_output, 'test')
    expect_is(test_output, 'data.frame')
    
    # check the data.frame contains the right columns
    test_igv_columns = c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')
    expect_true(all(test_igv_columns %in% names(test_output)))
    
    # check error if bad input
    expect_error(format_igv_seg(test_maf))
})

