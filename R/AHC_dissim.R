
#' Naive classification by dissimilarity
#'
#' @param D dissimilarity matrix
#' @param X covariate matrix (if you do not have the dissimilarity matrix, it will be calculated)
#' @param method Type of method used for the calculation of dissimilarities (e.g. "centroid", "median" etc...)
#'
#' @return
#' @export
#' @import ggfortify
#' @import ggplot2
#'
#' @examples methods = c("single","complete","average","weighted","centroid","median","Ward")
#' X = as.matrix(iris[-5])
#' iris.pca = prcomp(X,center=TRUE,scale.=TRUE)
#' table.dissim = list()
#'
#' for (method in methods){
#'   P = AHC.dissim(X=X,method=method)
#'   P3 = P[150-length(levels(iris$Species))+1,]
#'   table.dissim[[method]] = table(P3,iris$Species)
#'   print(method)
#'   print(table.dissim[[method]])
#'   print(ggplot2::autoplot(iris.pca,data=iris,colour=P3)+ggplot2::ggtitle(method))
#' }
#' print(ggplot2::autoplot(iris.pca,data=iris,colour="Species")+ggplot2::ggtitle("true species"))
#'
AHC.dissim = function(D=0,X=0,method="average"){
  if (!method %in% c("single","complete","average","weighted","centroid","median","Ward")) stop("the method given is unknown")
  if (!is.matrix(X) & !is.matrix(D)) stop("no values have been given to both D and X or the value(s) is (are) not a matrix")
  else if (!is.matrix(D)){
    n0 = dim(X)[1]
    D = matrix(0,n0,n0)
    for (j in 1:(n0-1)){
      for (i in (j+1):n0){
        D[i,j] = sum((X[i,]-X[j,])**2)
        D[j,i] = D[i,j]
      }
    }
  }
  else{
    if (!isSymmetric(D)) stop("D must be square symmetric")
  }
  diag(D) = Inf
  rownames(D) = 1:n0
  colnames(D) = 1:n0
  n0 = dim(D)[1]
  # On initialise les paramètres du clustering
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
  n = n0
  for (k in 2:n){
    # On cherche les 2 classes les plus proches
    arg.min = which.min(D)
    j = (arg.min-1)%/%n+1
    i = arg.min-(j-1)*n
    # On crée la nouvelle partition
    P[k,] = P[k-1,]
    P[k,P[k,]==i] = j
    P[k,P[k,]>i] = P[k,P[k,]>i]-1
    # On met à jour alpha et beta
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
    # On met à jour D
    D[j,] = alpha.i*D[i,]+alpha.j*D[j,]+beta*D[i,j]+gamma*abs(D[i,]-D[j,])
    D[j,j] = Inf
    D[,j] = D[j,]
    D = D[-i,-i]
    n = n-1
  }
  return(P)
}
