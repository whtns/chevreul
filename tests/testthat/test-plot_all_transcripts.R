test_that("multiplication works", {
  expect_error(
    plot_all_transcripts(human_gene_transcript_sce, "NRL"),
    NA)
})
