test_that("plotting works", {
  expect_error(
    plot_violin(human_gene_transcript_sce, plot_var = "batch", features = "NRL"),
    NA)
})
