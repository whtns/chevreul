test_that("object metadata retrieved", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  expect_type(get_object_metadata(chevreul_sce), "list")

})
