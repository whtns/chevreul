test_that("Feature type set", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce
  new_object<-set_feature_types(chevreul_sce, "transcript")
  expect_equal(mainExpName(new_object), "transcript")
})
