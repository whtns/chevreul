test_that("multiplication works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_equal(query_experiment(chevreul_sce, "gene"), TRUE)
})
