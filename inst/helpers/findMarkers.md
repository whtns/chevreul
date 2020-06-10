# Find Marker Genes

Marker genes can be plotted for each cluster at any pre-computed cluster resolution. Marker genes are determined by Wilcoxon rank sum test using [presto](https://www.biorxiv.org/content/10.1101/653253v1) or by specificity using [genesorteR](https://github.com/mahmoudibrahim/genesorteR). Any number of marker genes per cluster can be specified. Available markers for each cluster are sorted based on adjusted p value and thresholded at 5E-2, then the top n genes according to user input are chosed based on adjusted p value. 
