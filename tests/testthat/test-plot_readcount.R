test_that("plot gets made", {
		small_example_dataset <- object_calcn(small_example_dataset)
    expect_error(
        plot_readcount((small_example_dataset), return_plotly = TRUE),
        NA
    )
})
