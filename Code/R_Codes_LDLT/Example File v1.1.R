# Example of using CDFandPDFmvna_v3.R file to compute cumulative probability function of multivariate normal (MVN) distribution

# Load functions and global variables from CDFandPDFmvna_v2.1.R
source("C:/Users/gs27556/Box Sync/Conversion/Code/CDFandPDFmvna_v2.1.R")

# Mean of MVN
mu = c(1, 1.9, 1.1, 2, 0.9)

# Point where distribution function is to be computed
x = c(2, 1, 2, 2, 1)

# Variance matrix of MVN
V = matrix(c( 3,    0.5, -0.4,  0.3, -0.2,
              0.5,  2,    0.2, -0.3,  0.4,
              -0.4,  0.2,  3,    0.5, -0.4,
              0.3, -0.3,  0.5,  2,    0.2,
              -0.2,  0.4, -0.4,  0.2,  1), nrow = 5)

# Random seed to be used in function computation
seed = 1

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
