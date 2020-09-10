# Example of using CDFandPDFmvna_v3.R file to compute cumulative probability function of multivariate normal (MVN) distribution

# Load functions and global variables from CDFandPDFmvna_v2.1.R
source("CDFandPDFmvna_v3.R")

# Mean of MVN
mu = c(1, 1.9, 1.1, 2, 0.9)

# Variance matrix of MVN
V = matrix(c( 3,    0.5, -0.4,  0.3, -0.2,
              0.5,  2,    0.2, -0.3,  0.4,
              -0.4,  0.2,  3,    0.5, -0.4,
              0.3, -0.3,  0.5,  2,    0.2,
              -0.2,  0.4, -0.4,  0.2,  1), nrow = 5)

# Random seed to be used in function computation
seed = 1

# Point where distribution function is to be computed
x = c(2, 1, 2, 2, 1)

## Alternative settings
V = diag(4)

mu = c(0,0,0,0)

x = c(0,0,0,0)

## Compute truncation points
pm1 = ncol(V)
#mu = matrix(X %*% beta, nrow = pm1)
#above = rep(0, pm1-1)
#prob = double(pm1)
j = 3
Aj = -diag(pm1-1)
if(j==1){
  Aj = cbind(1,Aj)
}
if(j>1 && j<pm1){
  Aj = cbind(Aj[,1:(j-1)],1,Aj[,j:(pm1-1)])
}
if(j==pm1){
  Aj = cbind(Aj,1)
}

x = as.vector(-Aj %*% mu)
V = -Aj %*% V %*% t(-Aj)
mu = mu[-j]

result = cdfmvnanalytic(mu, V, x, seed)

# result is a list whose first item is the cumulative probability
result[[1]]


.method = "TVBS"  # Choose among "TG", "ME", "OVUS", "OVBS", "TGBME", "BME", "TVBS", "SSJ"

# If only the value of cumulative probability at x is to be calculated, use function cdfmvnanalytic
result = cdfmvnanalytic(mu, V, x, seed)

# result is a list whose first item is the cumulative probability
result[[1]]
# and second item is the random state that may be passed to the next function that uses random functions

# If both the cumulative probability and its gradient with respect to the function arguments is to be computed,
# use function pdfmvnanalytic
result = pdfmvnanalytic(mu, V, x, seed)

# result is a list whose first item is the cumulative probability
result[[1]]

# Second item is gradient vector with respect to mean
result[[2]]

# Third item is gradient vector with respect to variance
# The vector contains gradients with respect to elements of the upper triangular matrix of variance in row major order
result[[3]]

# Fourth item is gradient vector with respect to the point where cumulative probability is computed
result[[4]]

# and fifth item is the random state that may be passed to the next function that uses random functions
