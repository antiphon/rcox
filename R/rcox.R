#' Simulate Cox process in 2- or 3-dimensional box
#' 
#' Simulate the spatial Cox point process in 2- or 3-dimensional box 
#' (rectangular cuboid). Three algorithms are available:
#'
#' @import Rcpp
#' @export
#' @useDynLib rcox

rcox <- function(lambda,
                 n,
                 bbox=cbind(c(0,1), c(0,1), c(0,1)), 
                 iter = 1e4,
                 verb=FALSE) {
  
  if(!"coxintensity"%in%is(lambda)) stop("use lambda(...) to create intensity.")
  #
  d <- ncol(bbox)
  win <- unlist(bbox)
  
  if(!missing(n)) {# use MH
    xyz <- rcox_MH(n, win, lambda, iter, verb)
  }
  # use thinning
  else stop("give n.")
  # 
  xyz <- do.call(cbind, xyz)
  list(x=xyz, bbox=bbox, lambda=lambda)
}
