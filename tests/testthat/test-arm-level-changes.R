
test_annotation = arm_level_changes(test_facets_output$segs, test_facets_out$ploidy)

test_that('function returns output', {
    test_names = c('genome_doubled', 'fcna', 'weighted_fcna', 'aneuploidy_score', 'altered_arms')
    expect_true(all(test_names %in% names(test_annotation)))
})
