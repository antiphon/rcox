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
coxintensity2matrix <- function(x, nx=100, ny=NULL, bbox, at = NULL, W, ...) {
  
  if(missing(bbox)){
    if(!missing(W)) bbox <- cbind(W$xrange, W$yrange)
    else bbox <- apply(x$x, 2, range)
  } 
  gridded <- is.null(at)
  if(gridded) {
    if(x$d == 3) stop("'at' missing, will no evaluate on a 3D grid by default.")
    if(missing(ny)){
      bd <- apply(bbox, 2, diff)
      asp <- bd[2]/bd[1]
      ny <- round(nx * asp)
    }
    dx <- seq(bbox[1,1], bbox[2,1], length=nx)
    dy <- seq(bbox[1,2], bbox[2,2], length=ny)
    at <- as.matrix( expand.grid(dy, dx)[,2:1] )
  }
  out <- evaluate_lambda_c(x, at)
  if(gridded) {
    V <- t(matrix(out, ncol=nx))
    out <- list(Lambda=V, x=dx, y=dy)
    class(out) <- c("im0", is(out))
  }
  out
}



#' Simple Plot Method for Grid Field
#' 
#' @export
plot.im0 <- function(x, ..., scale = identity) {
  image(x$x, x$y, scale(x$Lambda), asp = 1, xlab="", ylab="", ...)
}
