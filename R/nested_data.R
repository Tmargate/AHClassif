#' Nested data simulation
#'
#' Simulates a nested dataset
#'
#' @param n number of points to be simulated
#' @param p size of the points to be simulated
#'
#' @return data matrix
#' @export
#'
#' @examples
nested_data = function(n,p){
  data = matrix(0,1,p)
  vec = diag(1,p)
  log2n = log(n,2)
  for (k in 0:ceiling(log2n)){
    for (i in 1:(2^k)){
      data = rbind(data,data[i,]+vec[k%%p+1,],data[i,]-vec[k%%p+1,])
    }
    vec = vec*0.5
    data = data[-(1:(2^k)),]
  }
  data = data[sample(1:(2^ceiling(log2n)),n),]
  return(data)
}
