
test_cnlr = cnlr_plot(test_facets_output, return_object = TRUE)
test_valor = valor_plot(test_facets_output, return_object = TRUE)
test_cf = cf_plot(test_facets_output, return_object = TRUE)
test_icn = icn_plot(test_facets_output, return_object = TRUE)
test_closeup = closeup_plot(test_facets_output, highlight_gene = 'TP53', return_object = TRUE)

test_that('cnlr plot works', {
    expect_is(test_cnlr, 'ggplot')
})

test_that('valor plot works', {
    expect_is(test_valor, 'ggplot')
})

test_that('cf plot works', {
    expect_is(test_cf, 'ggplot')
})

test_that('icn plot works', {
    expect_is(test_icn, 'ggplot')
})

test_that('closeup plot works', {
    expect_is(test_closeup, 'list')
    expect_is(test_closeup[[1]], 'ggplot')
})
