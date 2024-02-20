test_that("plotting works", {
  expect_error(
    plot_var(human_gene_transcript_sce, "batch"),
    NA)
})
