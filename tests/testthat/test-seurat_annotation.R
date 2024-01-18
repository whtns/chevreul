test_that("cell cycle annotated", {
    cell_cycle_annotated_seu <- annotate_cell_cycle(panc8)

    cell_cycle_cols <- c("S.Score", "G2M.Score", "Phase")

    expect_equal(c(colnames(panc8@meta.data), cell_cycle_cols), colnames(cell_cycle_annotated_seu@meta.data))
})
