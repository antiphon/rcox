#' Simulate Cox process in 2- or 3-dimensional box
#' 
#' Simulate the spatial Cox point process in 2- or 3-dimensional box 
#' (rectangular cuboid). Three algorithms are available:
#' 
#' @param lambda object from 'lambda'-function
#' @param n points to simulate if fixed count wanted
#' @param bbox bounding box, column matrix giving ranges
#' @param W owin-object rectangular, can be given instead of bbox
#' @param iter iterations of MH algorithm if n given
#' @param verb Print some runtime output
#' 
#' @import Rcpp
#' @export
#' @useDynLib rcox

rcox <- function(lambda,
                 n,
                 bbox,
                 W,
                 iter = 1e4,
                 verb=FALSE) {
  
  if(!"coxintensity"%in%is(lambda)) stop("use lambda(...) to create intensity.")
  #
  d <- lambda$d
  
  if(missing(bbox)){
    if(!missing(W)){if(is(W,"owin")) bbox <- cbind(W$xrange, W$yrange)}
    else stop("give bbox or W")
  }
  
  win <- unlist(bbox)
  
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

#' Helper to Transform to ppp
#' 
#' @export
cox2ppp <- function(x, ...){
  if(ncol(x$x) == 2)
    spatstat::ppp(x$x[,1], x$x[,2], window = spatstat::as.owin(c(x$bbox)))
  else if(ncol(x$x) == 3)
    spatstat::pp3(x$x[,1], x$x[,2], x$x[,3], spatstat::as.box3(c(x$bbox)))
}


