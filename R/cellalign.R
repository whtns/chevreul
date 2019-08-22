#' Align two trajectories via cellalign
#'
#' @param exp_ref
#' @param exp_query
#' @param traj_ref
#' @param traj_query
#'
#' @return
#' @export
#'
#'
#' @examples
local_align <- function(exp_ref, exp_query, traj_ref, traj_query) {
  #interpolating and scaling data for local alignment
  interLocalRef = cellAlign::interWeights(expDataBatch = exp_ref, trajCond = traj_ref, winSz = 0.1, numPts = numPts)
  interLocalQuery = cellAlign::interWeights(expDataBatch = exp_query, trajCond = traj_query, winSz = 0.1, numPts = numPts)
  interScaledLocalRef = cellAlign::scaleInterpolate(interLocalRef)
  interScaledLocalQuery = cellAlign::scaleInterpolate(interLocalQuery)

  #calculate dissimilarity matrix gated with threshold and align at local minima
  A=calcDistMat(interScaledLocalQuery$scaledData,interScaledLocalRef$scaledData, dist.method = 'Euclidean')
  A[A > 10*Thresh] <- max(A)
  alignment = localAlign(interScaledLocalQuery$scaledData,interScaledLocalRef$scaledData,threshPercent = Thresh)

  costMat = t(apply(A,1,function(x){return(as.numeric(x))}))
  linearInd = cellAlign::sub2ind(nrow(A), alignment$align[[1]]$index1, alignment$align[[1]]$index2)
  costMat[linearInd] = NA
  costMat = data.frame(costMat, row.names=1:numPts)
  colnames(costMat) = 1:numPts
  # pheatmap(costMat, cluster_cols = F, cluster_rows=F, border_color = NA,
  #          main = 'gated search region',
  #          show_rownames = F, show_colnames = F)

  alignment_plot <- plotAlign(alignment)

  #plots regions in pseudotime that are conserved
  BRef=colMeans(interScaledLocalRef$scaledData)
  BQuery=colMeans(interScaledLocalQuery$scaledData)
  p=unique(alignment$align[[1]]$index1)
  q=unique(alignment$align[[1]]$index2)
  plot(1:200,BQuery,xlab = "pseudotime",ylab = "mean interpolated expression", main = "unaligned mean expression",ylim = c(0,1.1))
  points(p,BQuery[p],col="red")
  points(1:200,BRef,col="grey60")
  points(q,BRef[q],col="red")
  text(90,1,"Ref")
  text(150,1,"Query")
  text(125,.3,"red points are conserved")

  return(alignment, alignment_plot)
}

#' Global Alignment of two trajectories with cellAlign
#'
#' @param exp_ref
#' @param exp_query
#' @param traj_ref
#' @param traj_query
#'
#' @return
#' @export
#'
#' @examples
global_align <- function(exp_ref, exp_query, traj_ref, traj_query){
  interRef = cellAlign::interWeights(expDataBatch = exp_ref, trajCond = traj_ref,
                                     winSz = 0.1, numPts = numPts)
  interQuery = cellAlign::interWeights(expDataBatch = exp_query, trajCond = traj_query,
                                       winSz = 0.1, numPts = numPts)

  #scale the interpolated data (Recommended):
  interScaledGlobalRef = cellAlign::scaleInterpolate(interRef)
  interScaledGlobalQuery = cellAlign::scaleInterpolate(interQuery)

  #test small DistMat
  A=calcDistMat(interScaledGlobalQuery$scaledData[,1:10],interScaledGlobalRef$scaledData[,1:10], dist.method = 'Euclidean')
  pheatmap(A, cluster_cols = F, cluster_rows=F, main = "Ref vs Query distances, 1st 10 points",
           show_rownames = F, show_colnames = F,display_numbers = TRUE)


  #perform global alignment of all genes
  alignment = globalAlign(interScaledGlobalQuery$scaledData, interScaledGlobalRef$scaledData,
                          scores = list(query = interScaledGlobalQuery$traj,
                                        ref = interScaledGlobalRef$traj),
                          sigCalc = F, numPerm = 20)
  p <- cellalign::plotAlign(alignment)
  #map interpolation to real data
  mapping = mapRealDataGlobal(alignment,intTrajQuery = interScaledGlobalQuery$traj, realTrajQuery = traj_query,
                              intTrajRef = interScaledGlobalRef$traj, realTrajRef = traj_ref)
  mapping_plot <- cellalign::plotMapping(mapping)
  print(p)
  print(mapping_plot)
  return(list(alignment, mapping))
}

#' Plot Gene expression over pseudotime on reference and query trajectories from cellAlign
#'
#' @param expGlobalRef
#' @param expGlobalQuery
#' @param trajRef
#' @param trajQuery
#'
#' @return
#' @export
#'
#' @examples
gene_test_plot <- function(expGlobalRef, expGlobalQuery, trajRef, trajQuery) {
  numPts = 200
  interGlobalRef = cellAlign::interWeights(expDataBatch = expGlobalRef, trajCond = trajRef,
                                           winSz = 0.1, numPts = numPts)
  interGlobalQuery = cellAlign::interWeights(expDataBatch = expGlobalQuery, trajCond = trajQuery,
                                             winSz = 0.1, numPts = numPts)

  # require(ggplot2)
  # require(reshape2)
  # require(pheatmap)
  sharedMarkers = intersect(rownames(expGlobalRef), rownames(expGlobalQuery))
  #whichgene="NRL"
  whichgene=sharedMarkers[1]
  selectedRef<-interGlobalRef$interpolatedVals[whichgene,]
  selectedQuery<-interGlobalQuery$interpolatedVals[whichgene,]

  dfRefi = data.frame(traj = interGlobalRef$traj, value=(selectedRef), error=interGlobalRef$error[whichgene,])
  dfRef = data.frame(traj = trajRef, t(expGlobalRef[whichgene,]))
  dfQueryi = data.frame(traj = interGlobalQuery$traj, value=(selectedQuery), error=interGlobalQuery$error[whichgene,])
  dfQuery = data.frame(traj = trajQuery, t(expGlobalQuery[whichgene,]))
  dfRefM = melt(dfRef, id.vars = 'traj')
  dfQueryM = melt(dfQuery, id.vars = 'traj')
  #plot of an example gene and its interpolation with error bars
  p <- ggplot(dfRefi, aes(x=traj,y=value)) +  geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + geom_line(size = 2) + geom_point(data=dfRefM, aes(x=traj,y=value)) + ggtitle(whichgene)
  print(p)
}
