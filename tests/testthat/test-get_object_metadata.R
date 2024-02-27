test_that("object metadata retrieved", {
  expect_equal(get_object_metadata(chevreul_sce), metadata(chevreul_sce))
})
