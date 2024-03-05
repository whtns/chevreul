test_that("object integrated", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    batches <- splitByCol(chevreul_sce, "batch")
    integrated_object <- object_integrate(batches)
    expect_match(mainExpName(integrated_object), "integrated")
})
