test_that("Integration pipeline successful", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    batches <- splitByCol(chevreul_sce, "batch")
    integrated_object <- object_integration_pipeline(batches)
    expect_match(mainExpName(integrated_object), "integrated")
})
