test_that("Read count metrics calculated", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(object_calcn(chevreul_sce), NA)
})
