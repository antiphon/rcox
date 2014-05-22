#' intensity object suitable for use in Cox simulation
#' 
#' @export

lambda <- function(x, max=5, kernel="gauss", sigma=0.5){
  l <- list(type="none")
  # if we have point locations to use for shot noise
  if(is.matrix(x)){
    # ok kernels
    accepted_kernels <- c("gaussian", "step")
    k <- which(match.arg(kernel, accepted_kernels)==accepted_kernels)
    l$kernel <- k-1
    
    d <- ncol(x)
    n <- nrow(x)
    l$max <- max
    l$x <- x
    l$sigma <- sigma
    l$type <- "shotnoise"
  }
  else stop("only shot noise implemented. Give coordinates x")
  class(l) <- "coxintensity"
  l
}

#' print lambda-object
#' 
#' @exportMethod print
print.coxintensity <- function(x, ...){
  cat("Intensity object for simulating Cox-process.\n")
  cat("Type:", x$type)
}