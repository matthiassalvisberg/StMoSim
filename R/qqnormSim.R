myQQNorm <- function(x, nSim){
  n_size <- length(x)
  
  if(n_size < 5){
    stop("too small sample")
  }
  
  n_sd <- mad(x)
  n_mean <- mean(x)
  
  a <- ifelse(n_size <= 10, 3/8, 1/2)
  myseq <- seq(from = (1 - a)/(n_size + (1-a)-a), to = (n_size - a)/(n_size + (1-a)-a), length.out = min(n_size,50))
  
  n_quant <- myQQNormIntern(myseq,n_mean,n_sd,length(x),nSim)
  
  myX <- qnorm((1:n_size - a)/(n_size + (1-a)-a))
  myY <- sort(x)
  
  plot(myX, myY, main = "Normal Q-Q Plot - SIM", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles") 
  matlines(qnorm(myseq), n_quant, lty = 1,pch = 1, lwd = 3, col = "#cdd2d015")
  
  points(myX, myY)
  qqline(x)
  box()
}


#' Plots a QQ-Norm plot with several Gaussian simulations.
#' 
#' Plots a QQ-Norm plot of the variable x with nSim Gaussian simulations.
#'
#' @param x is a lm-object or a numeric vector. If it's a lm-object its residuals are plotted.
#' @param nSim is an optional argument. If you like to have more or less than 500 simulations you can specify this parameter.
#'
#' @return invisible(NULL)
#' 
#' @export
#' @rdname qqnormSim
#'
#' @useDynLib StMoSim
#'
#' @importFrom graphics box matlines plot points
#' @importFrom stats mad qnorm qqline resid
#' @importFrom Rcpp evalCpp
#' @import RcppParallel 
#' @import methods 
#'
#' @examples
#' \dontrun{
#' # The observations should behave like a simulation, 
#' # because the observations are sampled from a Gaussian distribution.
#' qqnormSim(rnorm(100))
#' # On the first glance its obvious that this sample 
#' # doesn't originate from a Gaussian distribution due to the heavy tails.
#' qqnormSim(rt(100,df = 4))
#'
#' Reduce the simulation tracks from 500 to 50. (500 is default).
#' Not recommended unless you have not enough computation power.
#' qqnormSim(rnorm(100), nSim = 50)
#' }
#' @keywords qqnorm
#' 
#' @author Matthias Salvisberg <matthias.salvisberg@@gmail.com>
#' 
setGeneric("qqnormSim", function(x, nSim = 500) 
  standardGeneric("qqnormSim"))


#' @rdname qqnormSim
#' @aliases qqnormSim,lm-method
setMethod("qqnormSim","lm",
          function(x, nSim = 500){
            myQQNorm(resid(x), nSim)
            return(invisible(NULL))
          })

#' @rdname qqnormSim
#' @aliases qqnormSim,numeric-method
setMethod("qqnormSim","numeric",
          function(x, nSim = 500){
            myQQNorm(x, nSim)
            return(invisible(NULL))
          })
