# Test thomas simulation

devtools::load_all()

# Generate mothers
bbox <- cbind(0:1, 0:1) * 100

# Parameters for thomas
omega <- 5
kappa <- 10/V # ~ten per window
alpha <- 10 # ten childern average
#
V <- prod(apply(bbox, 2, diff))
bbox_ex <- bbox + c(-1,1) * (4*omega)
V_ex <- prod(apply(bbox_ex, 2, diff))

nmo <- rpois(1, V_ex * kappa)
mo <- apply(bbox_ex, 2, function(ab) runif(nmo, ab[1], ab[2]))

# Create lambda field object
L <- lambda(x = mo, kernel = "gauss", sigma = omega, alpha = alpha)

# Check it out
e <- coxintensity2matrix(L, bbox = bbox)

# Then simulate daughters
x <- lapply(1:50, function(i) rcox(L, bbox = bbox)$x)
nn<-sapply(x,nrow)
print(mean(nn))