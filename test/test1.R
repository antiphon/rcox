#test parameters
library(rcox)
x <- matrix(runif(30),ncol=2)
l <- lambda(x, sigma=0.01, kernel=k<-"step")

b<-rcox(lambda=l, n=400, bbox=cbind(c(0,1),c(0,1)), iter=it<-1e5)
plot(b$x, asp=1)
points(x, col=2, pch=19)


x <- matrix(runif(30),ncol=3)
l <- lambda(x, sigma=0.01, kernel=k)

b<-rcox(lambda=l, n=400, bbox=cbind(c(0,1),c(0,1), c(0,1)), iter=it)

library(rgl)
plot3d(b$x, asp=1)
points3d(x, col=2, pch=19)

