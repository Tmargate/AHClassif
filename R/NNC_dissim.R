#' Nearest Neighbors Chain algorithm (dissimilarity)
#'
#' Performs a ascending hierarchical classification using the Nearest-Neighbor Chain algorithm from the dissimilarity matrix.
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
#' P = NNC_dissim(X=X, method = "median")
#'
#'
NNC_dissim = function(D=0,X=0,method="average"){
  if (!is.matrix(D)){
    D = as.matrix(dist(X))
  }
  n0 = nrow(D)
  colnames(D) = 1:n0
  rownames(D) = NULL
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
  if (method == "single"){
    gamma = -0.5
  }
  else if (method == "complete"){
    gamma = 0.5
  }
  else if (method == "median"){
    beta = -0.25
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
      stack[2] = as.double(names(which.min(D[stack[1],-stack[1]])))
      l.stack = 2
    }
    #case 2 : the stack contains at least 2 elements
    if (l.stack>=2){
      NN = as.double(names(which.min(D[stack[l.stack],-stack[-(l.stack-1)]])))
      while (NN != stack[l.stack-1]){
        stack[l.stack+1] = NN
        l.stack = l.stack+1
        NN = as.double(names(which.min(D[stack[l.stack],-stack[-(l.stack-1)]])))
      }
      #merge of NN and stack[l.stack-1].
      P[k,] = P[k-1,]
      imin = min(stack[l.stack-1],stack[l.stack])
      imax = max(stack[l.stack-1],stack[l.stack])
      t[k-1] = D[imin,imax]
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
      D[imin,] = alpha1*D[imin,]+alpha2*D[imax,]+beta*D[imin,imax]+gamma*abs(D[imin,]-D[imax,])
      D[,imin] = D[imin,]
      D[imin,imin] = 0
      D = D[-imax,-imax]
      n = n-1
      if (n != 1){
        colnames(D) = 1:n
      }
    }
  }
  return(list(P=P,coph=t))
}
