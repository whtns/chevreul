test_that("plot_transcript_composition works", {
    expect_error(
        plot_transcript_composition(tiny_sce, "NRL"),
        NA
    )
})
