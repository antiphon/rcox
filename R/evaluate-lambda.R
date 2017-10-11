#' Evaluate lambda
#' 
#' @export
coxintensity2matrix <- function(x, nx=100, ny=NULL, bbox, W, ...) {
  if(x$d!=2) stop("Plotting only for 2D.")
  
  if(missing(bbox)){
    if(!missing(W)) bbox <- cbind(W$xrange, W$yrange)
    else bbox <- apply(x$x, 2, range)
  } 
  
  if(missing(ny)){
    bd <- apply(bbox, 2, diff)
    asp <- bd[2]/bd[1]
    ny <- round(nx * asp)
  }
  dx <- seq(bbox[1,1], bbox[2,1], length=nx)
  dy <- seq(bbox[1,2], bbox[2,2], length=ny)
  grid <- as.matrix( expand.grid(dy, dx)[,2:1] )
  v <- evaluate_lambda_c(x, grid)
  V <- t(matrix(v, ncol=nx))
  out <- list(Lambda=V, x=dx, y=dy)
  class(out) <- c("im0", is(out))
  out
}



#' Simple Plot Method 
#' 
#' @export
plot.im0 <- function(x, ..., scale = identity) {
  image(x$x, x$y, scale(x$Lambda), asp = 1, xlab="", ylab="", ...)
}
