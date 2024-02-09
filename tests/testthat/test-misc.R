test_that("x notation renaming works", {
    # Test case 1: Basic test
    cell_ids <- c("X1", "X22", "X333")
    batch_id <- "B"
    expected_result <- c("B001", "B022", "B333")
    result <- rename_from_x_notation(cell_ids, batch_id)
    expect_equal(result, expected_result)

    # Test case 2: Empty input
    cell_ids <- character(0)
    batch_id <- "B"
    expected_result <- "B"
    result <- rename_from_x_notation(cell_ids, batch_id)
    expect_equal(result, expected_result)

    # Test case 3: Single cell ID
    cell_ids <- "X9999"
    batch_id <- "B"
    expected_result <- "B9999"
    result <- rename_from_x_notation(cell_ids, batch_id)
    expect_equal(result, expected_result)

    # Test case 4: Batch ID with leading zeros
    cell_ids <- c("X1", "X22", "X333")
    batch_id <- "B001"
    expected_result <- c("B001001", "B001022", "B001333")
    result <- rename_from_x_notation(cell_ids, batch_id)
    expect_equal(result, expected_result)
})
