test_that("returns a character vector", {
    chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_type(
        metadata_from_object(chevreul_sce),
        "character"
    )
})
