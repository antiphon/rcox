#' test plotting lambda

devtools::load_all(".")

x0 <- matrix(runif(30),ncol=2)
x <- lambda(x0, sigma=0.05, kernel=k<-"step")


i <- coxintensity2matrix(x, bbox=cbind(0:1,0:1))
image(i$x,i$y,i$Lambda, useRaster = T, asp=1)

range(i$Lambda)

points(x0)