test_that("Varisble features retrieved", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  expect_error(get_variable_features(chevreul_sce), NA)
})
