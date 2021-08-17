context("species comparison")

test_that("gene symbols converted between species", {
  src_genes <- rownames(baron2016singlecell$gene)

  human_conversion <- convert_symbols_by_species(src_genes, "mouse")

  mouse_conversion <- convert_symbols_by_species(human_conversion, "human")

  expect_equal(src_genes, src_genes)
})
