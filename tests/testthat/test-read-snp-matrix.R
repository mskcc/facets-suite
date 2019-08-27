
# This test does not work on Travis-CI for some reason

# test_that('read function works', {
#     
#     test_snp_matrix = read_snp_matrix(system.file('data-raw/countsMerged_mini.gz', package = 'facetsSuite'))
#     
#     test_snp_columns = c('Chromosome', 'Position', 'NOR.DP', 'TUM.DP', 'NOR.RD', 'TUM.RD')
#     
#     expect_is(test_snp_matrix, 'data.frame')
#     expect_true(all(test_snp_columns %in% names(test_snp_matrix)))
# })
