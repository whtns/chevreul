test_that("Feature names obained", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    names <- get_features(chevreul_sce)
    expect_equal(names, rownames(chevreul_sce))
})
