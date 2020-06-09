# batches <- seurat_pancreas_reduced %>%
#   purrr::map(Seurat::SplitObject, split.by = "dataset") %>%
#   purrr::transpose() %>%
#   identity()
#
# test1 <- integration_workflow(batches)
