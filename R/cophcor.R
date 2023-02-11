#' Cophenetic correlation
#'
#' Computes the cophenetic correlation between two dendograms
#'
#' @param hclust1 first list having as elements a matrix P of partitions and a vector coph of kinetic distances
#' @param hclust2 second list having as elements a matrix P of partitions and a vector coph of kinetic distances
#'
#' @return cophenetic correlation between the dendrograms of hclust1 and hclust2
#' @export
#'
#' @importFrom stats cor
#'
#' @examples
cophcor = function(hclust1,hclust2){
  coph1 = hclust1$coph
  coph2 = hclust2$coph
  P1 = hclust1$P
  P2 = hclust2$P
  n = ncol(P1)
  t1 = matrix(0,n,n)
  t2 = matrix(0,n,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      k = 2
      while (P1[k,i] != P1[k,j]){
        k = k+1
      }
      t1[i,j] = coph1[k-1]
      k = 2
      while (P2[k,i] != P2[k,j]){
        k = k+1
      }
      t2[i,j] = coph2[k-1]
    }
  }
  c = cor(t1[upper.tri(t1)]-mean(t1[upper.tri(t1)]),t2[upper.tri(t2)]-mean(t2[upper.tri(t2)]))
  return(c)
}
