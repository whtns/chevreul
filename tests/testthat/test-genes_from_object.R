test_that("returns a character vector", {
  expect_type(
    genes_from_object(human_gene_transcript_seu),
    "character")
})
