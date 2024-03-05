test_that("plotting of heatmap works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        make_complex_heatmap(chevreul_sce, features = "NRL"),
        NA
    )
})
