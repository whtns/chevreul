test_that("returns a character vector", {
  expect_type(
    metadata_from_object(human_gene_transcript_sce),
    "character"
  )
})
