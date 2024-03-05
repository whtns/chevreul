test_that("SCE objects merged", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_contains(obj <- merge_small_objects(
        "small_batch1" = chevreul_sce[, seq(100)],
        "small_batch2" = chevreul_sce[, seq(101, 300)]
    ), list())
})
