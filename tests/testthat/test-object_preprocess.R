test_that("preprocessing works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        object_preprocess(chevreul_sce),
        NA
    )
})
