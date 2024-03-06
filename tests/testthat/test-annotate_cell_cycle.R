test_that("Annotation of readcount works", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    data(cc.genes.cyclone)
    expect_error(
        annotate_cell_cycle(chevreul_sce),
        NA
    )
})
