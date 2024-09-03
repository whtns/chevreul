test_that("plotting works for sce", {
   chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        plot_markers(chevreul_sce, group_by = "gene_snn_res.0.2"),
        NA
    )
})
