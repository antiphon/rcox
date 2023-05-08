## ---- message=FALSE, warning=FALSE---------------------------------------
library(rcox)
library(spatstat)

## ---- fig.width=6--------------------------------------------------------
set.seed(1) 
x1 <- runifpoint(50)
x2 <- rMaternII(120, 0.08)
generators <- listof(a=x1, b=x2)
plot(generators)

## ------------------------------------------------------------------------
lam1 <- lambda(x1, kernel = "gauss", sigma = 0.05, alpha = 3, type = "sum")
lam2 <- lambda(x2, kernel = "step", sigma = 0.1, alpha = c(log(100), -0.5), type = "prod")

## ------------------------------------------------------------------------
field1 <- coxintensity2matrix(lam1, W = square(1))
field2 <- coxintensity2matrix(lam2, bbox = cbind(0:1, 0:1), nx = 128)

## ---- fig.width=6--------------------------------------------------------
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(field1, axes =F)
points(x1)
plot(field2, axes =F, scale=log)
points(x2)

## ------------------------------------------------------------------------
s1 <- rcox(lambda = lam1, bbox = cbind(0:1, 0:1))
s2 <- rcox(lambda = lam2, n = 100, bbox = cbind(0:1, 0:1), iter = 1e4, verb=T)
# turn to ppp-objects
y1 <- cox2ppp(s1)
y2 <- cox2ppp(s2)

## ---- fig.width=6--------------------------------------------------------
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(y1)
points(x1, pch=3, col = 2)
plot(y2)
points(x2, pch=3, col = 2)

