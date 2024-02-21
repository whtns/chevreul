test_that("plot gets made", {

  expect_error(
    plot_readcount(human_gene_transcript_sce),
    NA)
})
