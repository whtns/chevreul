test_that("plot_transcript_composition works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        plot_transcript_composition(chevreul_sce, "NRL"),
        NA
    )
})
