test_that("clustering done", {
	chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(object_cluster(chevreul_sce), NA)
})
