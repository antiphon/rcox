## ---- message=FALSE, warning=FALSE---------------------------------------
library(rcox)

## ---- fig.width=6--------------------------------------------------------
set.seed(1) 
x1 <- matrix( runif( 10 * 3), ncol = 3 )
x2 <- cbind(0.5, 0.5, 0.5)
bb3 <- cbind(0:1, 0:1, 0:1) # bounding box for later

## ------------------------------------------------------------------------
lam1 <- lambda(x1, kernel = "gauss", sigma = 0.05, alpha = 5, type = "sum")
lam2 <- lambda(x2, kernel = "step", sigma = 0.1, alpha = c(1, 100), type = "prod")

## ------------------------------------------------------------------------
# make a slice z=0.5
ng <- seq(0, 1, l = 50)
grid3 <- as.matrix( expand.grid( ng, ng , 0.5) )
field1 <- coxintensity2matrix(lam1, bbox = bb3, at = grid3)
field2 <- coxintensity2matrix(lam2, bbox = bb3, nx = 128, at = grid3)

## ---- fig.width=6--------------------------------------------------------
par(mfrow=c(1,2))
plot3 <- function(f) image(matrix(f, ncol = length(ng)), x=ng, y=ng, asp=1) 
plot3(field1)
plot3(field2)

## ------------------------------------------------------------------------
s1 <- rcox(lambda = lam1, bbox = bb3)
s2 <- rcox(lambda = lam2, bbox = bb3, n = 100)
# turn to ppp-objects
y1 <- cox2ppp(s1)
y2 <- cox2ppp(s2)

## ---- fig.width=6--------------------------------------------------------
par(mfrow=c(1,2), mar=c(1,1,1,1))
plot(y1)
plot(y2)

