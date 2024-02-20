test_that("plot_transcript_composition works", {
  expect_error(
    plot_transcript_composition(human_gene_transcript_sce, "NRL"),
    NA)
})
