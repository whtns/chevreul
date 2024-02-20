test_that("plot gets made", {

  expect_error(
    plot_readcount(genes_to_transcripts),
    NA)
})
