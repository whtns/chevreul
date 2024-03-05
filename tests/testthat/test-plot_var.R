test_that("plotting works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        plot_var(chevreul_sce, "batch"),
        NA
    )
})
