
test_that('cnlr plot works', {
    
    # check that function returns ggplot object
    test_cnlr = cnlr_plot(test_facets_output, return_object = TRUE)
    expect_is(test_cnlr, 'ggplot')
})

test_that('valor plot works', {
    
    # check that function returns ggplot object
    test_valor = valor_plot(test_facets_output, return_object = TRUE)
    expect_is(test_valor, 'ggplot')
})

test_that('cf plot works', {
    
    # check that function returns ggplot object
    test_cf = cf_plot(test_facets_output, return_object = TRUE)
    expect_is(test_cf, 'ggplot')
})

test_that('icn plot works', {
    
    # check that function returns ggplot object
    test_icn = icn_plot(test_facets_output, return_object = TRUE)
    expect_is(test_icn, 'ggplot')
})

test_that('closeup plot works', {
    
    # check that function returns ggplot objects in a list
    test_closeup = closeup_plot(test_facets_output, highlight_gene = 'TP53', return_object = TRUE)
    expect_is(test_closeup, 'list')
    expect_is(test_closeup[[1]], 'ggplot')
})

test_that('help functions work', {
    
    # check expected output from get_cum_chr_maploc
    test_cumloc = get_cum_chr_maploc(test_facets_output$snps)
    expect_true(all(c('snps', 'mid', 'centromeres') %in% names(test_cumloc)))
    
    # check expected output from gene pos
    test_one_gene = get_gene_position('BRCA1')
    expect_true(nrow(test_one_gene) == 1)
    expect_true(all(test_one_gene$gene == 'BRCA1'))
    
    test_two_genes = get_gene_position(c('BRCA1', 'BRCA2'))
    expect_true(nrow(test_two_genes) == 2)
    expect_true(all(c('BRCA1', 'BRCA2') %in% test_two_genes$gene))
    
    expect_error(get_gene_position('no_gene'))
})