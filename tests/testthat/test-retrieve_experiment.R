test_that("Assay retrieved", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    retrieve_data <- retrieve_experiment(chevreul_sce, experiment = "gene")
    expect_equal(mainExpName(retrieve_data), "gene")
})
