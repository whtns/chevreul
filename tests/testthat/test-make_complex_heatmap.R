test_that("plotting of heatmap works", {
  expect_error(
    make_complex_heatmap(human_gene_transcript_sce, features = "NRL"),
    NA)
})
