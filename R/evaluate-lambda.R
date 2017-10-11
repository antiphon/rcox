#' Evaluate lambda
#' 
#' @export
coxintensity2matrix <- function(x, nx=100, ny=NULL, bbox, ...) {
  if(x$d!=2) stop("Plotting only for 2D.")
  
  if(missing(bbox)) bbox <- apply(x$x, 2, range)
  
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
  list(Lambda=V, x=dx, y=dy)
}

