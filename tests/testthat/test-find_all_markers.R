test_that("Markers found", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        find_all_markers(chevreul_sce),
        NA
    )
})
