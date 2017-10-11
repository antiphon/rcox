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
  d <- lambda$d
  win <- unlist(bbox[1:d, 1:d])
  
  if(!missing(n)) {# use MH
    xyz <- rcox_MH(n, win, lambda, iter, verb)
  }
  # use thinning
  else {
    # if sum field, maximum determined in C. 
    if(lambda$type == 0) {
      
    }
    else{ # for product
      #if not set
      # if(lambda$max == 0){
      #   v <- coxintensity2matrix(lambda)
      #   lambda$max <- max(v$Lambda)
      # }
    }
    xyz <- rcox_thin(win, lambda, verb)
  }
  # 
  xyz <- do.call(cbind, xyz)
  list(x=xyz, bbox=bbox, lambda=lambda)
}
