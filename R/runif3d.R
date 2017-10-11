#' Uniform points in a bounding box
#' 
#' @export
runif_bbox <- function(n, bbox=cbind(0:1,0:1,0:1)){
  x<-list(x=apply(bbox, 2, function(a)runif(n,a[1],a[2])), bbox=bbox)
  x
}