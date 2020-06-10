# Regress Features

Often we would like to correct for expression variation in our analysis due to phenomena unrelated to our focus of study. Such unwanted effects can include cell-cycle variability or differences in mitochondrial expression. 

Rather than simply excluding count values attributable to a given process, it is better to adjust the expression of all remaining genes or transcripts in each cell based on the sum score of relevant genes, to _regress_ out cell-cycle effects, for [example](https://satijalab.org/seurat/v3.0/cell_cycle_vignette.html)

Such regression can be performed for any arbitrary gene set. Some preset gene lists are included in `seuratTools` including cell-cycle, mitochondrial, and apoptotsis related gene sets. 
