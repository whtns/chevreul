test_that("pipeline runs", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        object_pipeline(chevreul_sce),
        NA
    )
})
