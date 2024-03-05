test_that("metadata united", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        unite_metadata(chevreul_sce, "nFeature_gene"),
        NA
    )
})
