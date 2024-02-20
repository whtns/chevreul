test_that("returns a character vector", {
  expect_contains(
    get_transcripts_from_object(human_gene_transcript_sce, "NRL"),
    c("ENST00000397002", "ENST00000561028"))
})
