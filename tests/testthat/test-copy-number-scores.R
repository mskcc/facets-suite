
test_that('fcna output works', {
    
    test_fcna = calculate_fraction_cna(test_facets_output$segs, test_facets_output$ploidy)
    
    test_fcna_columns = c('genome_doubled', 'fraction_cna')
    expect_true(all(test_fcna_columns %in% names(test_fcna)))
})

test_that('lst output works', {
    
    test_lst = calculate_lst(test_facets_output$segs, test_facets_output$ploidy)
    
    test_lst_columns = c('lst')
    expect_true(all(test_lst_columns %in% names(test_lst)))
})

test_that('hrd-loh output works', {
    
    test_hrdloh = calculate_hrdloh(test_facets_output$segs, test_facets_output$ploidy)
    
    test_hrdloh_columns = c('hrd_loh')
    expect_true(all(test_hrdloh_columns %in% names(test_hrdloh)))
})

test_that('ntai output works', {

    test_ntai = calculate_ntai(test_facets_output$segs, test_facets_output$ploidy)

    test_ntai_columns = c('ntelomeric_ai')
    expect_true(all(test_ntai_columns %in% names(test_ntai)))
})
