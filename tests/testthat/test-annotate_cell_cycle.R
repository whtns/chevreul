test_that("Annotation of readcount works", {
  expect_error(
    annotate_cell_cycle(human_gene_transcript_sce),
    NA)
})
