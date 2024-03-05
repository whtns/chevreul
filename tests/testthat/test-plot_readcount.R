test_that("plot gets made", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()

    expect_error(
        plot_readcount(chevreul_sce),
        NA
    )
})
