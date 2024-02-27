test_that("original experiment preserved as alt experiment", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  regressed_object <- regress_cell_cycle(chevreul_sce)
  expect_contains(mainExpName(regressed_object), "gene_regressed")
})
