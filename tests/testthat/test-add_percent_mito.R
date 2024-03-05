test_that("percentage mito added", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        add_percent_mito(chevreul_sce),
        NA
    )
})
