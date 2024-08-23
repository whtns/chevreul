test_that("plotting of heatmap works", {
    
    expect_error(
        make_complex_heatmap(chevreul_sce, features = "NRL"),
        NA
    )
})
