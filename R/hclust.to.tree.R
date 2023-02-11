#' Classification by hclust package
#'
#' Computes the partition matrix of an object in the hclust package.
#'
#' @param hclust object of the hclust package
#'
#' @return parition matrix
#' @export
#'
#' @examples
hclust.to.tree = function(hclust){
  n = nrow(hclust$merge)+1
  P = matrix(0,n,n)
  P[1,] = 1:n
  for (i in 1:(n-1)){
    clusters = rep(0,2)
    for (j in 1:2){
      if (hclust$merge[i,j] > 0){
        clusters[j] = P[i,which.max(P[hclust$merge[i,j],]-P[hclust$merge[i,j]+1,])]
      }
      else{
        clusters[j] = P[i,-hclust$merge[i,j]]
      }
    }
    imin = min(clusters)
    imax = max(clusters)
    P[i+1,] = P[i,]
    P[i+1,P[i+1,]==imax] = imin
    P[i+1,P[i+1,]>imax] = P[i+1,P[i+1,]>imax]-1
  }
  return(P)
}
