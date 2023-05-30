#test thinning simulation
devtools::load_all(".")
library(spatstat)
#set.seed(2)
#x <- with(rThomas(10,0.05, 6), cbind(x, y))
#l <- lambda(x, sigma=0.1, kernel=k<-"gauss", alpha = c(log(500), -.4), type = "prod")
l <- lambda(x, sigma=0.1, kernel=k<-"step", alpha = 1, type = "sum")
l2 <- lambda(x, sigma=0.1, kernel=k<-"step", alpha = 1, type = "sum", weights = runif(nrow(x)))

print(l2)


par(mfrow=c(2,1))
plot(l)
plot(l2)

