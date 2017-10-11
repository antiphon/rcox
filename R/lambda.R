#' intensity object suitable for use in Cox simulation
#' 
#' @param x matrix of coordinates
#' @param kernel type of kernel to use, "gauss" or "step"
#' @param sigma kernel parameter, sd of gauss or range of step
#' @param type "sum" or "product", see details
#' @param alpha additional scaling parameter
#' @details 
#' For type="sum", the field is 
#' 
#' v(u) = alpha * sum k(u-x)
#' 
#' where k is a kernel, and type="product"
#' 
#' v(u) = exp(alpha[1]) * prod [ 1 + alpha[2]*k(u-x)/k(0) ]
#' 
#' following Jalilian et al. 2015. The kernels are normalised densities in the dimension of x.
#' 
#' This function merely creates an field object. Use "coxintensity2matrix" to evaluate it at arbitrary locations, or "rcox" to simulate uniform/poisson process on it.
#' 
#' 
#' @export

lambda <- function(x, kernel="gauss", sigma=0.5, type = "sum", alpha = 1){
  l <- list(type="none")
  if(is(x, "ppp")) x <- cbind(x$x, x$y)
  # if we have point locations to use for shot noise
  if(is.matrix(x)){
    # ok kernels
    accepted_kernels <- c("gaussian", "step")
    k <- which(match.arg(kernel, accepted_kernels)==accepted_kernels)
    l$kernel <- k-1
    l$d <- ncol(x)
    l$n <- nrow(x)
    l$max <- 0  # we don't know yet
    l$x <- x
    l$type <- pmatch(type, c("sum", "product")) - 1
    if(l$type == 0) if(alpha < 0) stop("alpha>0 required for sum field to ensure positivity.")
    if(l$type == 1){
      if(length(alpha) != 2) stop("product field requires alpha of length 2.")
      if(alpha[2] < -1) stop("alpha[2] >= -1 required for product field to ensure positivity.")
    }
    l$alpha <- alpha
    l$sigma <- sigma
  }
  else stop("only shot noise implemented. Give coordinates x")
  class(l) <- "coxintensity"
  l
}

#' print lambda-object
#' 
#' @export
print.coxintensity <- function(x, ...){
  cat("Intensity object for simulating Cox-process.\n")
  cat("Dimension:", x$d, "\n")
  cat("Type:", c("sum", "product")[1+x$type], "\n")
  cat("Kernel:", c("gaussian", "step")[x$kernel + 1], "\n")
  cat("alpha:", paste(x$alpha, collapse = ", "), "\n")
  cat("sigma:", x$sigma, "\n")
}