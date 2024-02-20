test_that("multiplication works", {
  expect_error(
    unite_metadata(human_gene_transcript_sce, "nFeature_gene"),
    NA)
})
