test_that("plotting works", {
    expect_error(
    	plot_var(small_example_dataset, "Mutation_Status", return_plotly = FALSE),
        NA
    )
})
