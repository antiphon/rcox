#test thinning simulation
devtools::load_all(".")
library(spatstat)
#set.seed(2)
x <- with(rThomas(10,0.05, 6), cbind(x, y))
#l <- lambda(x, sigma=0.1, kernel=k<-"gauss", alpha = c(log(500), -.4), type = "prod")
l <- lambda(x, sigma=0.1, kernel=k<-"step", alpha = c(log(50),1), type = "prod")

bb <- cbind(c(0,1),c(0,1)) + c(0,0)

b <- rcox(lambda=l, n=50,bbox=bb, verb = T)

lam <- coxintensity2matrix(l, bbox = bb)

image(lam$x, lam$y, (lam$Lambda), asp=1, col=gray.colors(120))
points(b$x, col=4, pch=20)
points(x, col=1, pch=19)

l$max <- b$lambda$max

# do many, check intensity
xl <- sapply(1:100, function(i) rcox(l, bbox=bb)$x)
print(mean(sapply(xl, nrow))/prod(apply(bb,2,diff)))
