test_that("Feature type retrived", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce

  expect_error(get_feature_types(chevreul_sce), NA)
})
