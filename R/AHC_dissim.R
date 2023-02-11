
#' Classification by dissimilarities matrix
#'
#' Performs a ascending hierarchical classification by classical algorithm from the dissimilarities matrix
#'
#' @param D dissimilarities matrix if already computed
#' @param X Data matrix (only if the similarity matrix S is not given)
#' @param method method for calculating dissimilarities between clusters (e.g. "single", "complete" ..), default = "average"
#'
#' @return list containing partition matrix and vector of cophenetic distances
#' @export
#'
#' @importFrom stats dist
#'
#' @examples
#' X = as.matrix(iris[-5])
#' method = "single"
#' P = AHC_dissim(X=X, method = "median")
#'
AHC_dissim = function(D=0,X=0,method="average"){
  if (!is.matrix(D)){
    D = as.matrix(dist(X))
  }
  diag(D) = Inf
  n0 = nrow(D)
  rownames(D) = 1:n0
  colnames(D) = 1:n0
  n0 = dim(D)[1]
  # the clustering parameters are initialized
  alpha.i = 0.5
  alpha.j = 0.5
  beta = 0
  gamma = 0
  if (method == "single"){
    gamma = -0.5
  }
  else if (method == "complete"){
    gamma = 0.5
  }
  else if (method == "median"){
    beta = -0.25
  }
  P = matrix(0,n0,n0)
  P[1,] = 1:n0
  t = rep(0,n0-1)
  n = n0
  for (k in 2:n0){
    # We look for the two closest classes
    arg.min = which.min(D)
    j = (arg.min-1)%/%n+1
    i = arg.min-(j-1)*n
    # Create the new partition
    P[k,] = P[k-1,]
    P[k,P[k,]==i] = j
    P[k,P[k,]>i] = P[k,P[k,]>i]-1
    # Update alpha and beta
    if (method %in% c("average","centroid","Ward")){
      card.i = sum(P[k-1,]==i)
      card.j = sum(P[k-1,]==j)
      if (method %in% c("average","centroid")){
        alpha.i = card.i/(card.i+card.j)
        alpha.j = card.j/(card.i+card.j)
        if (method == "centroid"){
          beta = -alpha.i*alpha.j
        }
      }
      else if (method == "Ward"){
        alpha.i = rep(card.i,n)
        alpha.j = rep(card.j,n)
        beta = rep(0,n)
        for (m in 1:n){
          card.m = sum(P[k-1,]==m)
          S = card.i+card.j+card.m
          alpha.i[m] = (alpha.i[m]+card.m)/S
          alpha.j[m] = (alpha.j[m]+card.m)/S
          beta[m] = -card.m/S
        }
      }
    }
    t[k-1] = D[i,j]
    # Update D
    D[j,] = alpha.i*D[i,]+alpha.j*D[j,]+beta*D[i,j]+gamma*abs(D[i,]-D[j,])
    D[j,j] = Inf
    D[,j] = D[j,]
    D = D[-i,-i]
    n = n-1
  }
  return(list(P=P,coph=t))
}
