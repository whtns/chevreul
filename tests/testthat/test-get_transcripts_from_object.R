test_that("returns a character vector", {
	chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_contains(
        get_transcripts_from_object(chevreul_sce, "NRL"),
        c("ENST00000397002", "ENST00000561028")
    )
})
