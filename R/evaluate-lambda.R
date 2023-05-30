#' Evaluate lambda
#' 
#' @param x object from 'lambda'-function
#' @param at where to evaluate, matrix of coordinates
#' @param nx if 'at' missing, make a grid with this many x-steps (only in 2d!)
#' @param ny like nx but for y-dimension (can be omitted)
#' @param bbox bounding box for the window in case grid used
#' @param W alternative to bbox, owin-object (2D!)
#' @param ... omitted.
#' 
#' @export
coxintensity2matrix <- function(x, nx=2^8, 
                                ny=NULL, 
                                bbox, at = NULL, W, ...) {
  
  if(missing(bbox)){
    if(!missing(W)) bbox <- cbind(W$xrange, W$yrange)
    else bbox <- apply(x$x, 2, range)
  } 
  gridded <- is.null(at)
  if(gridded) {
    if(x$d == 3) stop("'at' missing, will not evaluate on a 3D grid by default.")
    if(missing(ny)){
      bd <- apply(bbox, 2, diff)
      asp <- bd[2]/bd[1]
      ny <- round(nx * asp)
    }
    dx <- seq(bbox[1,1], bbox[2,1], length=nx)
    dy <- seq(bbox[1,2], bbox[2,2], length=ny)
    at <- as.matrix( expand.grid(dy, dx)[,2:1] )
  }
  # weighting:
  weights <- if(!is.null( x$weights) ) x$weights else 1
  x$x <- cbind(x$x, weights)
  # compute
  out <- evaluate_lambda_c(x, at)
  # 
  if(gridded) {
    V <- t(matrix(out, ncol=nx))
    out <- list(Lambda=V, x=dx, y=dy)
    class(out) <- c("im0", is(out))
  }
  # store parameters
  attr(out, "data") <- x
  attr(out, "pars") <- list(nx = nx, ny = ny, bbox = bbox)
  # done
  out
}

#' print for coxintensitymatrix
#' 
#' @param x coxintensity2matrix result
#' @param ... ignored
#' 
#' @export
print.im0 <- function(x, ...) {
  pars <- attr(x, "pars")
  cat("Evaluated intensity object for simulating Cox-process.\n")
  cat(sprintf("bbox: [%f,%f]x[%f,%f]\n", 
      pars$bbox[1,1], pars$bbox[2,1], pars$bbox[1,2], pars$bbox[2,2]))
  cat(sprintf("res: %i x %i\n", pars$nx, pars$ny))
}

#' Simple Plot Method for Grid Field
#' 
#' @param x result object of coxintensity2matrix()
#' @param ... passed on to image()
#' @param scale scaling function applied to values. Default: identity
#' @param points.cex if >0, add locations
#' @param points.col col of points
#' @param points.pch pch fo points
#' 
#' @export
plot.im0 <- function(x, ..., scale = identity, 
                     points.cex = 0, points.pch = 19, points.col = 1) {
  image(x$x, x$y, scale(x$Lambda), asp = 1, xlab="", ylab="", ...)
  if(points.cex > 0) points(
    attr(x, "data")$x,
    cex = points.cex,
    pch = points.pch,
    col = points.col
  )
}
