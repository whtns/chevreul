test_that("maintains object type", {
    expect_type(record_experiment_data(small_example_dataset), typeof(small_example_dataset))
})

test_that("experiment metadata is a list", {
   
    expect_type(
        get_object_metadata(record_experiment_data(small_example_dataset))[["experiment"]],
        "list"
    )
})

test_that("names are correct", {
   
    expect_named(
        get_object_metadata(record_experiment_data(small_example_dataset))[["experiment"]],
        c(
            "experiment_name", "organism", "date_of_export", "date_of_analysis",
            "parameters", "filtering", "sessionInfo", "SingleCellExperiment_version",
            "chevreul_version"
        )
    )
})
