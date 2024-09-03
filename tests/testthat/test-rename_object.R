test_that("Project renamed", {
    obj <- rename_object(small_example_dataset, "new_name")
    expect_contains(metadata(obj)["project.name"], "new_name")
})
