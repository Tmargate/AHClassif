#' Similarities matrix
#'
#' Computes the similarities matrix from a data matrix
#'
#' @param X Data matrix
#' @param Kernel.FUN  kernel function to be used if we want to calculate S from X (default :"linear")
#' @param ... other parameters that will be passed to the Kernel procedure (see Kernel procedure for details)
#'
#' @return similarity matrix calculated from X
#' @export
#'
#' @examples
sim = function(X,Kernel.FUN="linear",...){
  n0 = dim(X)[1]
  S = diag(1,n0)
  for (i in 1:(n0-1)){
    for (j in (i+1):n0){
      S[i,j] = Kernel(X[i,],X[j,],Kernel.FUN,...)/2+0.5
      S[j,i] = S[i,j]
    }
  }
  return(S)
}
