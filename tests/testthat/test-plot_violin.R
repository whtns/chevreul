test_that("plotting works", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  expect_error(
    plot_violin(chevreul_sce, plot_var = "batch", features = "NRL"),
    NA)
})
