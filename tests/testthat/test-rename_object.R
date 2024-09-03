test_that("Project renamed", {
   chevreul_sce <- small_example_dataset
    obj <- rename_object(chevreul_sce, "new_name")
    expect_contains(metadata(obj)["project.name"], "new_name")
})
