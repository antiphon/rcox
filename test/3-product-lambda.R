#' test plotting lambda, product kernel

devtools::load_all(".")

x0 <- matrix(runif(30),ncol=2)


x <- lambda(x0, sigma=0.05, kernel=k<-"gauss", type ="product", alpha=1)
y <- lambda(x0, sigma=0.05, kernel=k<-"gauss", type ="sum", alpha=1)

i <- coxintensity2matrix(x, bbox=cbind(0:1,0:1), nx = 2^7)
i2 <- coxintensity2matrix(y, bbox=cbind(0:1,0:1), nx = 2^7)

print(c(range(i$Lambda), range(i2$Lambda)))

par(mfrow=c(2,1))
image(i$x,i$y,i$Lambda, useRaster = T, asp=1, main="product")
points(x0)
image(i2$x,i2$y,i2$Lambda, useRaster = T, asp=1, main="sum")
points(x0)