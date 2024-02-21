test_that("percentage mito added", {
  expect_error(
    add_percent_mito(human_gene_transcript_sce),
    NA)
})
