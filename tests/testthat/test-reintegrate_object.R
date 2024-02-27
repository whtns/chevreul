test_that("Object reintegrated", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  expect_error(
    reintegrate_object(chevreul_sce),
    NA)
})
