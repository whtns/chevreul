test_that("plotting works for sce", {
   chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        plot_feature(chevreul_sce, embedding = "UMAP", features = "NRL"),
        NA
    )
})
