test_that("original experiment preserved as alt experiment", {
  regressed_object <- regress_cell_cycle(human_gene_transcript_sce)
  expect_contains(altExpNames(regressed_object, "original"))
})
