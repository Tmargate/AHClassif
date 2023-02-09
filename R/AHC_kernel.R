#' Title
#'
#' @param x
#' @param y
#' @param FUN
#' @param .gamma
#' @param c0
#' @param d
#' @param normalize
#'
#' @return
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
