test_that("integrated experiment created", {
  batches <- splitByCol(human_gene_transcript_sce, "batch")
  integrated_object <- integration_workflow(batches)

  expect_match(mainExpName(integrated_object), "integrated")
})
