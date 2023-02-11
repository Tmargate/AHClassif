#' Kernel functions for similarities
#'
#' Computes the value of a kernel between two points
#'
#' @param x coordinate vector of the first point
#' @param y coordinate vector of the second point
#' @param FUN kernel function (e.g. "chi2", "laplacian" etc..)
#' @param .gamma hyperparameters of the kernel function (default = 1)
#' @param c0 hyperparameters of the kernel function (default = 0)
#' @param d hyperparameters of the kernel function (default = 2)
#' @param normalize Boolean indicating whether the output value should be normalized (default = True)
#'
#' @return kernel value in (x,y)
#' @export
#'
#' @examples
Kernel = function(x,y,FUN,.gamma=1,c0=0,d=2,normalize=TRUE){
  if (FUN == "chi2"){
    return(exp(-.gamma*sum((x-y)**2/(x+y))))
  }
  else if (FUN == "laplacian"){
    return(exp(-.gamma*sum(abs(x-y))))
  }
  else if (FUN == "linear"){
    if (normalize == TRUE){
      return(sum(x*y)/sqrt(sum(x**2)*sum(y**2)))
    }
    else{
      return(sum(x*y))
    }
  }
  else if (FUN == "polynomial"){
    if (normalize == TRUE){
      return((.gamma*sum(x*y)+c0)**d/sqrt(((.gamma*sum(x**2)+c0)*(.gamma*sum(y**2)+c0))**d))
    }
    else{
      return((.gamma*sum(x*y)+c0)**d)
    }
  }
  else if (FUN == "rbf"){
    return(exp(-.gamma*sum((x-y)**2)))
  }
  else if (FUN == "sigmoid"){
    if (normalize == TRUE){
      return(tanh(.gamma*sum(x*y)+c0)/sqrt(tanh(.gamma*sum(x**2)+c0)*tanh(.gamma*sum(y**2)+c0)))
    }
    else{
      return(tanh(.gamma*sum(x*y)+c0))
    }
  }
  else{
    print("Error : the value for FUN is not recognized")
    return()
  }
}
