#' Nearest Neighbors Chain algorithm (similarity)
#'
#' Performs a ascending hierarchical classification using the Nearest-Neighbor Chain algorithm from the similarity matrix.
#'
#' @param S Similarity matrix if already computed
#' @param X Data matrix (only if the similarity matrix S is not given)
#' @param sparsity threshold below which similarities are zero
#' @param method method for calculating dissimilarities between clusters (e.g. "single", "complete" ..), default = "average"
#' @param Kernel.FUN kernel function to be used if we want to calculate S from X (default :"linear")
#' @param ... other parameters that will be passed to the Kernel procedure (see Kernel procedure for details)
#'
#' @return list containing partition matrix and vector of cophenetic distances
#' @export
#'
#' @examples
#' X = as.matrix(iris[-5])
#' method = "single"
#' P = NNC_sim(X=X, method = "median", Kernel.FUN = "sigmoid")
#'
NNC_sim = function(S=0,X=0,sparsity=0,method="average",Kernel.FUN="linear",...){
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
  colnames(S) = 1:n0
  P = matrix(0,n0,n0)
  P[1,] = 1:n0
  t = rep(0,n0-1)
  card = rep(1,n0)
  stack = c()
  l.stack = 0
  n = n0
  #parameters of the LW formula
  alpha1 = 0.5
  alpha2 = 0.5
  beta = 0
  gamma = 0
  delta1 = 0.5
  delta2 = 0.5
  if (method == "single"){
    gamma = -0.5
  }
  else if (method == "complete"){
    gamma = 0.5
  }
  else if (method == "median"){
    beta = -0.25
    delta1 = 0.25
    delta2 = 0.25
  }
  #loop
  for (k in 2:n0){
    #case 0: the stack is empty
    if (l.stack==0){
      stack[1] = sample(1:n,1)
      l.stack = 1
    }
    #case 1 : the stack contains 1 element
    if (l.stack==1){
      S.active = S[stack[1],-stack[1]]
      nonzero = S.active != 0
      stack[2] = as.double(names(which.max(S.active[nonzero]-0.5*diag(S)[-stack[1]][nonzero])))
      l.stack = 2
    }
    #case 2 : the stack contains at least 2 elements
    if (l.stack>=2){
      S.active = S[stack[l.stack],-stack[-(l.stack-1)]]
      nonzero = S.active != 0
      NN = as.double(names(which.max(S.active[nonzero]-0.5*diag(S)[-stack[-(l.stack-1)]][nonzero])))
      while (NN != stack[l.stack-1]){
        stack[l.stack+1] = NN
        l.stack = l.stack+1
        S.active = S[stack[l.stack],-stack[-(l.stack-1)]]
        nonzero = S.active != 0
        NN = as.double(names(which.max(S.active[nonzero]-0.5*diag(S)[-stack[-(l.stack-1)]][nonzero])))
      }
      #merge of NN and stack[l.stack-1].
      P[k,] = P[k-1,]
      imin = min(stack[l.stack-1],stack[l.stack])
      imax = max(stack[l.stack-1],stack[l.stack])
      t[k-1] = S[imin,imin]+S[imax,imax]-2*S[imin,imax]
      P[k,P[k,]==imax] = imin
      P[k,P[k,]>imax] = P[k,P[k,]>imax]-1
      #removed from the stack
      stack = stack[-c(l.stack-1,l.stack)]
      l.stack = l.stack-2
      stack[stack>imax] = stack[stack>imax]-1
      #computes parameters of the LW formula
      if (method %in% c("average","centroid")){
        sum = card[imin]+card[imax]
        alpha1 = card[imin]/sum
        alpha2 = card[imax]/sum
        if (method == "centroid"){
          beta = -alpha1*alpha2
          delta1 = alpha1^2
          delta2 = alpha2^2
        }
      }
      else if (method == "Ward"){
        sum = card+card[imin]+card[imax]
        alpha1 = card+card[imin]/sum
        alpha2 = card+card[imax]/sum
        beta = 1-alpha1-alpha2
      }
      card[imin] = card[imin]+card[imax]
      card = card[-imax]
      #compute the new distances
      S.ii = delta1*S[imin,imin]+delta2*S[imax,imax]
      S[imin,] = alpha1*S[imin,]+alpha2*S[imax,]+beta*S[imin,imax]-gamma*abs(S[imin,]-S[imax,])
      S[,imin] = S[imin,]
      S[imin,imin] = S.ii
      S = S[-imax,-imax]
      n = n-1
      if (n != 1){
        colnames(S) = 1:n
      }
    }
  }
  return(list(P=P,coph=t))
}
