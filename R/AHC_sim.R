
#' Naive classification by similarity
#'
#' @param S similarity matrix
#' @param X covariate matrix (if you do not have the similarity matrix, it will be calculated)
#' @param method type of method used for the calculation of similarities (e.g. "centroid", "median" etc...)
#' @param Kernel.FUN kernel function used for similarity calculation (e.g. "polynomial","rbf","sigmoid" etc...)
#' @param ...
#'
#' @return
#' @export
#'
#' @import ggfortify
#' @import ggplot2
#'
#' @examples methods = c("single","complete","average","weighted","centroid","median","Ward")
#' X = as.matrix(iris[-5])
#' iris.pca = prcomp(X,center=TRUE,scale.=TRUE)
#' table.sim = list()
#'
#' for (method in methods){
#'   P = AHC.sim(X=X,method=method)
#'   P3 = P[150-length(levels(iris$Species))+1,]
#'   table.sim[[method]] = table(P3,iris$Species)
#'   print(method)
#'   print(table.sim[[method]])
#'   print(ggplot2::autoplot(iris.pca,data=iris,colour=P3)+ggplot2::ggtitle(method))
#' }
#' print(ggplot2::autoplot(iris.pca,data=iris,colour="Species")+ggplot2::ggtitle("true species"))
#'
AHC.sim = function(S=0,X=0,method="average",Kernel.FUN="linear",...){
  if (!method %in% c("single","complete","average","weighted","centroid","median","Ward")) stop("the method given is unknown")
  if (!is.matrix(X) & !is.matrix(S)) stop("no values have been given to both S and X or the value(s) is (are) not a matrix")
  else if (!is.matrix(S)){
    n0 = dim(X)[1]
    S = diag(1,n0)
    for (i in 1:(n0-1)){
      for (j in (i+1):n0){
        S[i,j] = Kernel(X[i,],X[j,],FUN=Kernel.FUN,...)
        S[j,i] = S[i,j]
      }
    }
  }
  else{
    if (!isSymmetric(S)) stop("S must be square symmetric")
  }
  rownames(S) = 1:n0
  colnames(S) = 1:n0
  n0 = dim(S)[1]
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
  n = n0
  for (k in 2:n){
    # On cherche les 2 classes les plus proches
    arg.max = which.max(S-0.5*(matrix(rep(diag(S),n),n,n)+matrix(rep(diag(S),n),n,n,byrow=TRUE))-diag(Inf,n))
    j = (arg.max-1)%/%n+1
    i = arg.max-(j-1)*n
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
    if (!isSymmetric(S)){
      print(k)
      print(arg.max)
      print(i)
      print(j)
      print(S)
      break
    }
    # On met à jour S
    S.jj = S[j,j]
    S[j,] = alpha.i*S[i,]+alpha.j*S[j,]+beta*S[i,j]-gamma*abs(S[i,]-S[j,])
    S[j,j] = delta.i*S[i,i]+delta.j*S.jj
    S[,j] = S[j,]
    S = S[-i,-i]
    n = n-1
  }
  return(P)
}
