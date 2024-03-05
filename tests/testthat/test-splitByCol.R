test_that("SCE split", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_contains(
        obj <- splitByCol(chevreul_sce, "batch"),
        list()
    )
})
