# Differential expression

Differential expression can be used to identify genes which vary significantly between two populations of cells. Testing is executed via `Seurat::FindMarkers` and can be specified to use several methods including:
1. Student's t-test
2. Wilcoxon rank-sum test
3. Likelihood ratio test assuming an underlying negative binomial distribution
4. [MAST](https://github.com/RGLab/MAST): GLM-framework that treates cellular detection rate as a covariate ([Finak et al, Genome Biology, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4676162/)) ([Installation instructions](https://github.com/RGLab/MAST))

Comparison groups can be selected based on pre-computed cluster information or can be specified based on lasso-selection from a dimensionally reduced plot. Differential expression results for each method can be displayed in an associated table. 
