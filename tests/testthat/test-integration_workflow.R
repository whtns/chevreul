test_that("integrated experiment created", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    batches <- splitByCol(chevreul_sce, "batch")
    integrated_object <- integration_workflow(batches)

    expect_match(mainExpName(integrated_object), "integrated")
})
