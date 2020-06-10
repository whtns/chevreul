# Pseudotime Analysis by [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/)

Pseudotime analysis is broken into several steps, reflecting a typical command-line workflow. Starting with a dimensionally reduced dataset: 

1. (Optionally) subset seurat embedding
2. Calculate trajectory graph using Monocle3
3. Order cells in pseudotime after specifying 'root' cells
4. Identify genes which vary signicantly over pseduotime based on autocorrelation over pseodtime. This analysis is performed assuming a negative binomial model of gene expression. 
5. Identify modules of genes which vary similarly over pseudotime

Module expression is displayed at cluster or cell level resolution in an associated heatmap. If cell-level expression, cells are ordered based on pseudotime values.  
