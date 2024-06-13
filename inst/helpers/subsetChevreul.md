# Subset a SingleCellExperiment Object

Often we want to subset a dataset to remove outlier cells. Such cells may be present for technical or biological reasons. If excluded for principled reasons, such filtering can allow meaningful analysis. 

Dataset can be subset either by:
1. Selection on a dimensionally reduced plot
2. Selection via a correctly formatted csv file 

Selection on a dimensionally-reduced plot can be accomplished by first selecting cells using the lasso tool then selecting rows of cells in the associated table and clicking 'Subset by Selection'

Selection via an uploaded csv can be accomplished by formatting a csv file with cells to retain listed in the first column then uploading and selecting 'Subset by CSV'

If the dataset to be subset consists of a single batch (not integrated) subsetting will perform clustering, marker gene, and pathway enrichment as defined in `seurat_pipeline`

If the dataset to be subset is already integrated, subsetting will result in integration based on the 'batch' variable then proceed with clustering, marker gene, and pathway enrichment as defined in `seurat_integration_pipeline`

