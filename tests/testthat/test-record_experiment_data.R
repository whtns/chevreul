test_that("maintains object type", {
   chevreul_sce <- small_example_dataset
    expect_type(record_experiment_data(chevreul_sce), typeof(chevreul_sce))
})

test_that("experiment metadata is a list", {
   chevreul_sce <- small_example_dataset
    expect_type(
        get_object_metadata(record_experiment_data(chevreul_sce))[["experiment"]],
        "list"
    )
})

test_that("names are correct", {
   chevreul_sce <- small_example_dataset
    expect_named(
        get_object_metadata(record_experiment_data(chevreul_sce))[["experiment"]],
        c(
            "experiment_name", "organism", "date_of_export", "date_of_analysis",
            "parameters", "filtering", "sessionInfo", "SingleCellExperiment_version",
            "chevreul_version"
        )
    )
})
