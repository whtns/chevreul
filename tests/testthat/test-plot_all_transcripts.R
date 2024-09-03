test_that("multiplication works", {
	chevreul_sce <- chevreuldata::human_gene_transcript_sce()
    expect_error(
        plot_all_transcripts(chevreul_sce, "NRL"),
        NA
    )
})
