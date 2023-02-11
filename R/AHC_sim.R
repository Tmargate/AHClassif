
#' Classification by similarities
#'
#' Performs a ascending hierarchical classification by the classical algorithm from the similarity matrix
#'
#' @param S Similarity matrix if already computed
#' @param X Data matrix (only if the similarity matrix S is not given)
#' @param sparsity threshold below which similarities are zero
#' @param method method for calculating dissimilarities between clusters (e.g. "single", "complete" ..), default = "average"
#' @param Kernel.FUN kernel function to be used if we want to calculate S from X (default :"linear")
#' @param ... other parameters that will be passed to the Kernel procedure (see Kernel procedure for details)
#'
#' @return P : list containing parition matrix and vector of kophenetic distances
#'
#' @export
#'
#' @examples
#' X = as.matrix(iris[-5])
#' method = "single"
#' P = AHC_sim(X=X, method = "median", Kernel.FUN = "sigmoid")
#'
AHC_sim = function(S=0,X=0,sparsity=0,method="average",Kernel.FUN="linear",...){
  if (!is.matrix(S)){
    n0 = dim(X)[1]
    S = diag(1,n0)
    for (i in 1:(n0-1)){
      for (j in (i+1):n0){
        S[i,j] = Kernel(X[i,],X[j,],Kernel.FUN,...)/2+0.5
        S[j,i] = S[i,j]
      }
    }
  }
  S[S<sparsity] = 0
  n0 = nrow(S)
  rownames(S) = 1:n0
  colnames(S) = 1:n0
  # On initialise les paramètres du clustering
  alpha.i = 0.5
  alpha.j = 0.5
  beta = 0
  gamma = 0
  delta.i = 0.5
  delta.j = 0.5
  if (method == "single"){
    gamma = -0.5
  }
  else if (method == "complete"){
    gamma = 0.5
  }
  else if (method == "median"){
    beta = -0.25
    delta.i = 0.25
    delta.j = 0.25
  }
  P = matrix(0,n0,n0)
  P[1,] = 1:n0
  t = rep(0,n0-1)
  n = n0
  for (k in 2:n0){
    # On cherche les 2 classes les plus proches
    zeros = S == 0
    S[zeros] = NA
    diag.S = diag(S)
    diag(S) = NA
    arg.max = which(t(t(S)-0.5*diag.S)-0.5*diag.S == max(t(t(S)-0.5*diag.S)-0.5*diag.S,na.rm=T),arr.ind=T)
    S[is.na(S)] = 0
    diag(S) = diag.S
    i = arg.max[1,1]
    j = arg.max[1,2]
    # On crée la nouvelle partition
    P[k,] = P[k-1,]
    P[k,P[k,]==i] = j
    P[k,P[k,]>i] = P[k,P[k,]>i]-1
    # On met à jour alpha, beta et delta
    if (method %in% c("average","centroid","Ward")){
      card.i = sum(P[k-1,]==i)
      card.j = sum(P[k-1,]==j)
      if (method %in% c("average","centroid")){
        alpha.i = card.i/(card.i+card.j)
        alpha.j = card.j/(card.i+card.j)
        if (method == "centroid"){
          beta = -alpha.i*alpha.j
          delta.i = alpha.i**2
          delta.j = alpha.j**2
        }
      }
      else if (method == "Ward"){
        alpha.i = rep(card.i,n)
        alpha.j = rep(card.j,n)
        beta = rep(0,n)
        for (m in 1:n){
          card.m = sum(P[k-1,]==m)
          sum = card.i+card.j+card.m
          alpha.i[m] = (alpha.i[m]+card.m)/sum
          alpha.j[m] = (alpha.j[m]+card.m)/sum
          beta[m] = -card.m/sum
        }
      }
    }
    t[k-1] = S[i,i]+S[j,j]-2*S[i,j]
    # On met à jour S
    S.jj = S[j,j]
    S[j,] = alpha.i*S[i,]+alpha.j*S[j,]+beta*S[i,j]-gamma*abs(S[i,]-S[j,])
    S[j,j] = delta.i*S[i,i]+delta.j*S.jj
    S[,j] = S[j,]
    S = S[-i,-i]
    n = n-1
  }
  return(list(P=P,coph=t))
}
