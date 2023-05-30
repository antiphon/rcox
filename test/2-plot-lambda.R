#' test plotting lambda

devtools::load_all(".")

# draw some random locations
x0 <- matrix(runif(100*2),ncol=2)

#plot(x0)

# create the model
x <- lambda(x0, sigma=0.05, kernel=k<-"gauss")

print(x)

# evaluation:
i <- coxintensity2matrix(x, bbox=cbind(0:1,c(0,1) ))
print(i)
plot(i)

# Should also work directly
plot(x, nx = 2^8, points.cex = .2)

