# The following code was originally written for GAUSS by Bhat (2018); see http://www.caee.utexas.edu/prof/bhat/CodeRepository/CODES/LDLT/GAUSS_Codes_LDLT.zip.
# The GAUSS code has been translated below to R by Bhat's PhD student, Gopindra Nair.
  
#*****************************************************************************************************************************************************
#                   Please reference Bhat, 2018 and provide due credit if you use this source code, in part or in the whole; THANKS!
#            Bhat, C.R. (2018), "New Matrix-Based Methods for the Analytic Evaluation of the Multivariate Cumulative Normal Distribution Function," 
#                                                  Transportation Research Part B, 109, 238-256.
#****************************************************************************************************************************************************/

# Defining new functions and operators that function similar to the functions and operators in GAUSS

library(mvtnorm)
library(Matrix)

`%|%` = function(x, y) {
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  if (length(x) == 1) x = matrix(x)
  if (length(y) == 1) y = matrix(y)
  return(rbind(x, y))
}

`%~%` = function(x, y) {
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  if (length(x) == 1) x = matrix(x)
  if (length(y) == 1) y = matrix(y)
  return(cbind(x, y))
}

`%*~%` = function(x, y) { 
  if (length(x) == 1) x = matrix(x)
  if (length(y) == 1) y = matrix(y)
  n = nrow(x) 
  out <- matrix(0, n, ncol(x) * ncol(y)) 
  for(i in 1:n) out[i, ] <- c(y[i, ] %o% x[i, ]) 
  return(out) 
} 

`%**%` = function(x, y) {
  if (length(x) == 1 & length(y) != 1) {
    x = matrix(x)
    x = x[1, 1, drop=T]
    return(x*y)
  } else if (length(y) == 1 & length(x) != 1) {
    y = matrix(y)
    y = y[1, 1, drop=T]
    return(x*y)
  } else {
    return(x%*%y)
  }
}

`%//%` = function(x, y) {
  if (length(x) == 1) {
    x = matrix(x)
    x = x[1, 1, drop=T]
  }
  if (length(y) == 1) {
    y = matrix(y)
    y = y[1, 1, drop=T]
  }
  return(x / y)
}

`%.*%` = function(x, y) {
  
  if (length(x) == 1 & length(y) != 1) {
    x = matrix(x)
    x = x[1, 1, drop=T]
    return(x*y)
  }
  
  if (length(y) == 1 & length(x) != 1) {
    y = matrix(y)
    y = y[1, 1, drop=T]
    return(x*y)
  }
  
  if (all(dim(x) == dim(y)))
    return(x*y)
  
  if (dim(x)[1] == 1 & dim(y)[1] > 1) {
    return(y %*% diag(x[,]))
  }
  
  if (dim(y)[1] == 1 & dim(x)[1] > 1) {
    return(x %*% diag(y[,]))
  }
  
  if (dim(x)[2] == 1 & dim(y)[2] > 1) {
    return(diag(x[,]) %*% y)
  }
  
  if (dim(y)[2] == 1 & dim(x)[2] > 1) {
    return(diag(y[,]) %*% x)
  }
  
}

gss.reshape = function(x, row, col) {
  return(matrix(t(x), row, col, byrow=T))
}

gss.eye = function(x) {
  return(diag(x[1]))
}

gss.submat = function(x, rows, cols) {
  if (length(rows) == 1) if (rows == 0) return(x[, t(cols), drop=F])
  if (length(cols) == 1) if (cols == 0) return(x[t(rows),, drop=F])
  return(x[t(rows), t(cols), drop = F])
}

gss.seqa = function(start, inc, n) {
  return(matrix(seq.int(from=start, by=inc, length.out = n), ncol = 1))
}

gss.zeros = function(row, col) {
  return(matrix(0, row, col))
}

gss.ones = function(row, col) {
  return(matrix(1, row, col))
}

gss.vec = function(x) {
  return(matrix(x, ncol = 1))
}

gss.vecr = function(x) {
  return(matrix(t(x), ncol = 1))
}

gss.diag = function(x) {
  # This may work with only square matrices unlike in GAUSS
  return(matrix(x[as.logical(diag(nrow(x)))], ncol = 1))
}

gss.diagrv = function(x, v) {
  # This may work with only square matrices unlike in GAUSS
  x[as.logical(diag(nrow(x)))] = v
  return(x)
}

gss.minc = function(x) {
  return(matrix(apply(x, 2, min), ncol = 1))
}

gss.maxc = function(x) {
  return(matrix(apply(x, 2, max), ncol = 1))
}

gss.prodc = function(x) {
  return(matrix(apply(x, 2, prod), ncol = 1))
}

gss.sumc = function(x) {
  return(matrix(colSums(x)))
}

gss.rev = function(x) {
  return(x[nrow(x):1,, drop=F])
}

gss.sortc = function(x, key) {
  ordering = order(x[,key])
  return(x[ordering,,drop=F])
}

gss.minindc = function(x) {
  return(matrix(apply(x, 2, which.min), ncol = 1))
}

gss.maxindc = function(x) {
  return(matrix(apply(x, 2, which.max), ncol = 1))
}

gss.cdfbvn = function(limit1, limi2, rho) {
  # This function is not vectorized like in GAUSS
  cur.seed = .Random.seed
  corr = diag(2)
  corr[1, 2] = rho
  corr[2, 1] = rho
  p = pmvnorm(upper = c(limit1, limi2), corr = corr)[[1]]
  .Random.seed <<- cur.seed
  return(p)
}

gss.cdftvn = function(limit1, limi2, limit3, rho12, rho23, rho13) {
  # This function is not vectorized like in GAUSS
  cur.seed = .Random.seed
  corr = diag(3)
  corr[1, 2] = rho12
  corr[2, 1] = rho12
  corr[1, 3] = rho13
  corr[3, 1] = rho13
  corr[2, 3] = rho23
  corr[3, 2] = rho23
  p = pmvnorm(upper = c(limit1, limi2, limit3), corr = corr, abseps=10^-8)[[1]]
  .Random.seed <<- cur.seed
  return(p)
}

gss.upmat = function(x) {
  return(x * as.integer(upper.tri(x, diag = T)))
}

gss.selif = function(x, e) {
  return(x[as.logical(e),,drop=F])
}

gss.delif = function(x, e) {
  return(x[!as.logical(e),,drop=F])
}

gss.delrows = function(x, r) {
  return(x[-r,,drop=F])
}

gss.ldl = function(x) {
  U = chol(x)   ## start from factor given by chol 
  D = diag(U)   ## extract sqrt(D) 
  L = t(U/D)    ## get unit lower triangular factor 
  D = diag(D^2) ## and diagonal 
  return(list(L, D))
}

gss.rndu = function(r, c, state=NULL) {
  if (is.null(state)) {
    return(matrix(runif(r*c), r, c, byrow = T))
  } else {
    if (length(state) == 1) set.seed(state)
    else .Random.seed <<- state
    rands = matrix(runif(r*c), r, c, byrow = T)
    return(list(rands, .Random.seed))
  }
}

gss.rndn = function(r, c, state=NULL) {
  if (is.null(state)) {
    return(matrix(rnorm(r*c), r, c, byrow = T))
  } else {
    if (length(state) == 1) set.seed(state)
    else .Random.seed <<- state
    rands = matrix(rnorm(r*c), r, c, byrow = T)
    return(list(rands, .Random.seed))
  }
}

#/* setting global variables */
.x1symmetric = 0
.x2symmetric = 0
.x1diagonal = 0
.x2diagonal = 0
.x2correlation = 0
.xinvsymmetric = 0
.xinvcorrelation = 0
.xinvdiagonal = 0
.omsymmetric = 0
.omdiagonal = 0
.omegacorr = 0
.condcov = 0
.condspecialcov = 0
.cholesky = 0
.cholcov = 0
.condcovspecial = 0
.condcovmeantrunc = 0
.condcovsigtrunc = 0
.optimal = 3
.covarr = 1
.perms = 20
.method = "OVUS"


#/* **
#    Format                              Purpose       
#** ===================================================================
#**  cdfmvnanalytic(mu,cov,x,s)          Analytically approximated CDF of multivariate normal            
#**  pdfmvnanalytic(mu,cov,x,s)          Gradient of approximated CDF w.r.t parameters
#
#
#Global:  .covarr = 1     means cov is a covariance matrix; .covar=0 means cov is actually a correlation matrix
#         .perms = n      means n permutations of abscissae will be used in the Switzer, Solow, Joe analytic approach, n=1 means only one permutation will be used.
#                         For all non-SSJ methods, only one sequence (based on global .optimal) will be used, and .perms is irrelevant
#         .optimal = 0    means non-SSJ methods will be based on simple ascending order of abscissae in evaluation; .optimal=1 means non-SJJ methods will be based
#                         on an ordering following the GGE approach (see Gibson et al., 1994 and Bhat, 2018 paper); .optimal=2 means that the abscissae are randomly ordered 
#                         before implementing the non-SSJ method; .optimal=3 means the abscissae are used in the same order as given. The value of .optimal does not matter for the SSJ method;
#                         IMPORTANT NOTE: For the TVBS method, only the .optimal=0,2,or 3 values have to be used to get correct analytic gradients, because the TVBS method already is 
#                         based on an approximate quadrivariate CDF function evaluation; this requirement is already taken care of in the code,
#                         by putting .optimal=0 (if .optimal is provided as 1) if the TVBS method is used.  
#
#         .method =       "SSJ"                     - Switzer, Solow, and Joe Method
#                         "TG"                      - Trinh and Genz's univariate conditioning approximation procedure
#                         "ME"                      - The traditional ME approach, implemented in a new matrix-based and LDLT-based manner
#                         "OVUS"                    - One-variate univariate screening approach
#                         "OVBS"                    - One-variate bivariate screening approach
#                         "TGBME"                   - Trinh and Genz's bivariate conditioning approximation procedure
#                         "BME"                     - Bivariate ME approach
#                         "TVBS"                    - Two-variate bivariate screening approach
#
#**  mu = c( mu1,mu2,mu3,mu4 )                        Kx1 (K>=2) vector of means
#
#**  cov = { cov11  cov12  cov13   cov14,             cov is KxK covariance matrix (if _covarr=1)
#            cov12  cov22  cov23   cov24,
#            cov13  cov23  cov33   cov34,
#            cov14  cov24  cov34   cov44    };
#
#**  cov = {   1    rho12  rho13   rho14,             cov is KxK correlation matrix (if _covarr=0)
#            rho12    1    rho23   rho24,
#            rho13  rho23    1     rho34,
#            rho14  rho24  rho34     1    };
#            
#**  x = c( x1,x2,x3,x4 )                             Kx1 vector of abscissae
#             
#**  s is a seed value that is relevant only for the SSJ method; This argument is optional
#
#OUTPUT
#
#**  pdfmvnanalytic(x,mu,cov)[[1]] = Distribution function                   
#                                                     1x1 scalar
#
#**  pdfmvnanalytic(x,mu,cov)[[2]] = Gradient of distribution function with respect to mu
#                         c( dP/dmu1,                 Kx1 vector
#                           dP/dmu2,
#                           dP/dmu3,
#                           dP/dmu4 )                           
#
#**  pdfmvnanalytic(x,mu,cov)[[3]] = Gradient of distribution function with respect to elements of covariance/correlation matrix
#                        c( dP/dcov11                 [K*(K+1)/2 x 1] vector if .covar=1
#                           dP/dcov12
#                           dP/dcov13
#                           dP/dcov14
#                           dP/dcov22
#                           dp/dcov23
#                           dP/dcov24
#                           dP/dcov33
#                           dP/dcov34
#                           dP/dcov44 )   
#
#                        c( dP/drho12                 [K*(K-1)/2 x 1] vector if .covar=0
#                           dP/drho13
#                           dP/drho14
#                           dP/drho23
#                           dP/drho24
#                           dP/drho34  )
#
##**  pdfmvnanalytic(x,mu,cov)[[4]] = Gradient of distribution function with respect to x
#                        c( dP/dx1,                   Kx1 vector
#                           dP/dx2,
#                           dP/dx3,
#                           dP/dx4 )   
#
#**  pdfmvnanalytic(x,mu,cov)[[5]] = 
#           Vector seed coming out of SSJ method or from random ordering of other methods (having this, and using this as seed for 
#           next call of an MVNCD evaluation is helpful for model estimation where multiple MVNCD evaluations have to be undertaken); 
#           if not SSJ method and not random ordering for other methods, s1 is the same as the input scalar seed s   */


cdfmvnanalytic = function(mu,cov,x,s=.Random.seed) {
  
  mu = matrix(mu)
  x = matrix(x)
  .method = tolower(.method)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  
  kk1 = (x-mu)%//%sqrtom
  
  #/* this next line is needed, because truncated values go bizarre if(truncation happens on abscissa less than -5.8 */) {
  kk1 = -5.8%**%(kk1 < -5.8)+kk1%.*%(kk1>=-5.8)      
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #       correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(nrow(kk1),1))   
  
  if(.method == "ssj") {
    output = cdfmvnassj(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "tg") {
    output = cdfmvnatg(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "me") {
    output = cdfmvname(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "ovus") {
    output = cdfmvnaovus(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "ovbs") {
    output = cdfmvnaovbs(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "tgbme") {
    output = cdfmvnatgbme(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "bme") {
    output = cdfmvnabme(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)
  } else if(.method == "tvbs") {
    if(.optimal==1) {
      .optimal <<-0
    }
    output = cdfmvnatvbs(t(kk1),kk2,s)
    p = output[[1]]
    s1 = output[[2]]
    output = list(p,(s1))
    return(output)        
  }
}   


pdfmvnanalytic = function(mu,cov,x,s) {
  
  mu = matrix(mu)
  x = matrix(x)
  .method = tolower(.method)
  
  s1=s
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  #/* this next line is needed, because truncated values go bizarre if(truncation happens on abscissa less than -6 */) {
  kk1 = -5.8%**%(kk1 < -5.8)+kk1%.*%(kk1>=-5.8)    
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #       correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(nrow(kk1),1))      
  if(.method == "ssj") {
    output = pdfmvnassj(t(kk1),kk2,s)
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "tg") {
    output = pdfmvnatg(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "me") {
    output = pdfmvname(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "ovus") {
    output = pdfmvnaovus(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "ovbs") {
    output = pdfmvnaovbs(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "tgbme") {
    output = pdfmvnatgbme(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "bme") {
    output = pdfmvnabme(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  } else if(.method == "tvbs") {
    if(.optimal==1) {
      .optimal <<-0
    }
    output = pdfmvnatvbs(t(kk1),kk2,s) 
    p = output[[1]]
    gw = output[[2]]
    grho = output[[3]]
    s1 = output[[4]]
  }
  if(.covarr) {
    output = gradcorcov(kk1,sqrtom,kk2)
    gbcorcov = output[[1]]
    gomegacorcov = output[[2]]
    gcov = gbcorcov%**%(gw)+gomegacorcov%**%grho
  } else {
    gcov = grho
  }
  gmu = -(gw)%//%sqrtom
  gx = -gmu    
  output = list(p,gmu,gcov,gx,(s1))
  return(output)    
}   

#/*****************************************************************************************************************************************************
#                             SET 1: THE CODES BELOW CONSTITUTE THE MAIN ROUTINES FOR THE DIFFERENT METHODS
#               cdfmvnaXXXXXX - approximating multivariate CDF (and corresponding gradient functions) using different approaches 
#*****************************************************************************************************************************************************/

#/*****************************************************************************************************************************************************
#                                                         The TG Method
#*****************************************************************************************************************************************************/

cdfmvnatg = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  p = gss.zeros(m,1) 
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =pnorm(atemp[1, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  for(h in seq.int(1,m-1,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    l=l[2:(m-h+1),2:(m-h+1), drop = F]
    d=d[2:(m-h+1),2:(m-h+1), drop = F]
    p[h+1] =pnorm((atemp[h+1, drop = F]-mutemp[1, drop = F])%//%sqrt(d[1,1, drop = F]))
  }
  output = list(gss.prodc(p),s1)
  return(output)
}


pdfmvnatg = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  kkk = 1 #/* dimension of truncation = 1 for ME */
  p = gss.zeros(m,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal == 1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  
  p[1] =pnorm(atemp[1, drop = F])
  gp1 = dnorm(atemp[1, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  gc = gss.zeros(m,1)
  gcumulfrompc={}
  
  for(h in seq.int(1,m-1,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])        
    mutilde = output[[1]]
    omega = output[[2]]
    if(h==1) {
      .condcov<<-0
      .condcovmeantrunc<<-0
    } else {
      .condcov<<-1
      .condcovmeantrunc<<-1
    }
    
    output = gcondmeantrunc(gss.eye(m-h),mutemp,sigtemp,atemp[h, drop = F])
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcov(gss.eye(m-h),sigtemp)
    gy1 = output[[1]]
    gxcov = output[[2]]
    gmucov=gss.zeros(m-h+1,((m-h)%**%(m-h+1)%//%2))
    gccov = gss.zeros(1,((m-h)%**%(m-h+1)%//%2))
    
    if(h==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-1,(m-1+((m-1)%**%(m)%//%2)))
    } else if(h > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[h,] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    l=l[2:(m-h+1),2:(m-h+1), drop = F]
    d=d[2:(m-h+1),2:(m-h+1), drop = F]
    sigtemp=l%**%d%**%t(l)
    p[h+1] =noncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[h+1, drop = F])
    output = gradnoncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[h+1, drop = F])
    gfrompmu = output[[1]]
    gfrompcov = output[[2]]
    gfrompc = output[[3]]
    grhorho = grhorho+(1%//%p[h+1, drop = F])%**%(gcumulmusig[,1, drop = F]%**%gfrompmu+gcumulmusig[,m-h+1, drop = F]%**%gfrompcov)
    
    gcumulfrompc = gcumulfrompc%|%gfrompc
    gc = gc+(1%//%p[h+1, drop = F])%**%(gcumulc[,1:kkk, drop = F]%**%gfrompmu+gcumulc[, (m-h+1):((m-h+kkk%**%(kkk+1)%//%2)), drop = F]%**%gfrompcov)
    
  }
  gcumulfrompc = gp1%|%gcumulfrompc
  gc = gc+(gcumulfrompc%//%p)
  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}


#/*****************************************************************************************************************************************************
#                                                         The ME Method
#*****************************************************************************************************************************************************/

cdfmvname = function(a,rr,s) {
  
  a = matrix(a, 1)
  
  m=ncol(a)
  p = gss.zeros(m,1) 
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =pnorm(atemp[1, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  for(h in seq.int(1,m-1,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)
    l = output[[1]]
    d = output[[2]]
    p[h+1] =pnorm((atemp[h+1, drop = F]-mutemp[1, drop = F])%//%sqrt(d[1,1, drop = F]))
  }
  output = list(gss.prodc(p),s1)
  return(output)
}

pdfmvname = function(a,rr,s) {
  
  a = matrix(a, 1)
  
  m=ncol(a)
  kkk = 1 #/* dimension of truncation = 1 for me */
  p = gss.zeros(m,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)        
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  
  p[1] =pnorm(atemp[1, drop = F])
  gp1 = dnorm(atemp[1, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  gc = gss.zeros(m,1)
  gcumulfrompc={}
  for(h in seq.int(1,m-1,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])        
    mutilde = output[[1]]
    omega = output[[2]]
    if(h==1) {
      .condcovsigtrunc<<-0
      .condcovmeantrunc<<-0
    } else {
      .condcovsigtrunc<<-1
      .condcovmeantrunc<<-1
    }
    
    output = gcondmeantrunc(gss.eye(m-h),mutemp,sigtemp,atemp[h, drop = F])
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcovtrunc(mutemp,sigtemp,atemp[h, drop = F])
    gmucov = output[[1]]
    gxcov = output[[2]]
    gccov = output[[3]]
    
    if(h==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-1,(m-1+((m-1)%**%(m)%//%2)))
    } else if(h > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[h,] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)
    l = output[[1]]
    d = output[[2]]
    sigtemp = l%**%d%**%t(l)
    p[h+1] =noncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[h+1, drop = F])
    output = gradnoncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[h+1, drop = F])
    gfrompmu = output[[1]]
    gfrompcov = output[[2]]
    gfrompc = output[[3]]
    grhorho = grhorho+(1%//%p[h+1, drop = F])%**%(gcumulmusig[,1, drop = F]%**%gfrompmu+gcumulmusig[,m-h+1, drop = F]%**%gfrompcov)
    
    gcumulfrompc = gcumulfrompc%|%gfrompc
    gc = gc+(1%//%p[h+1, drop = F])%**%(gcumulc[,1:kkk, drop = F]%**%gfrompmu+gcumulc[, (m-h+1):((m-h+kkk%**%(kkk+1)%//%2)), drop = F]%**%gfrompcov)
    
  }
  gcumulfrompc = gp1%|%gcumulfrompc
  gc = gc+(gcumulfrompc%//%p)
  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}

#/*****************************************************************************************************************************************************
#                                                         The OVUS Method
#*****************************************************************************************************************************************************/

cdfmvnaovus = function(a,rr,s) { 
  
  a = matrix(a, nrow=1)
  
  m=ncol(a)
  p = gss.zeros(m-1,1) 
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output  = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  
  for(h in seq.int(1,m-2,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)        
    l = output[[1]]
    d = output[[2]]
    
    sigtemp1 = l[1:2,1:2, drop = F]%**%d[1:2,1:2, drop = F]%**%(t((l[1:2,1:2, drop = F])))
    p[h+1] = noncdfbvn(mutemp[1:2, drop = F],sigtemp1,(t((atemp[(h+1):(h+2), drop = F]))))%//%(noncdfn(mutemp[1, drop = F],sigtemp1[1,1, drop = F],atemp[h+1, drop = F]))   
    
  }
  
  output = list(gss.prodc(p),s1)
  return(output)
}


pdfmvnaovus = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  kkk = 2 #/* denotes the dimensionality of the numerator in p[h+1] computation */
  p = gss.zeros(m-1,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  rr10=sigtemp
  
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]

  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  gc = gss.zeros(m,1)
  gcumulfrompc=gss.zeros(m,1)
  
  for(h in seq.int(1,m-2,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])        
    mutilde = output[[1]]
    omega = output[[2]]
    if(h==1) {
      .condcovsigtrunc <<-0
      .condcovmeantrunc <<-0
    } else {
      .condcovsigtrunc <<-1
      .condcovmeantrunc <<-1
    }
    
    output = gcondmeantrunc(diag(m-h),mutemp,sigtemp,atemp[h, drop = F])
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcovtrunc(mutemp,sigtemp,atemp[h, drop = F])
    gmucov = output[[1]]
    gxcov = output[[2]]
    gccov = output[[3]]
    
    if(h==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-1,(m-1+((m-1)%**%(m)%//%2)))
    } else if(h > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[h,] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)
    l = output[[1]]
    d = output[[2]]
    sigtemp = l%**%d%**%t(l) 
    
    p[h+1] = noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],(t((atemp[(h+1):(h+2), drop = F]))))%//%(noncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[h+1, drop = F]))
    
    output = gradnoncdfbvnbycdfn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],(t((atemp[(h+1):(h+2), drop = F]))))
    gfrompmu = output[[1]]
    gfrompcov = output[[2]]
    gfrompc = output[[3]]
    grhorho = grhorho+(1%//%p[h+1, drop = F])%**%(gcumulmusig[,1:2, drop = F]%**%gfrompmu+(gcumulmusig[, (m-h+1):(m-h+2), drop = F]%~%gcumulmusig[,2%**%m-2%**%h+1, drop = F])%**%gfrompcov)
    
    gcumulfrompc[(h+1):(h+2)] = gcumulfrompc[(h+1):(h+2), drop = F]+(gfrompc%//%p[h+1, drop = F])  #/* getting gradient contribution directly from noncdfbvn%//%noncdfn for absiccae except the first two */ 
    gc = gc+(1%//%p[h+1, drop = F])%**%(gcumulc[,1:2, drop = F]%**%gfrompmu+(gcumulc[, (m-h+1):(m-h+2), drop = F]%~%gcumulc[,2%**%m-2%**%h+1, drop = F])%**%gfrompcov)
    
  }
  output = gradcdfbvn(atemp[1, drop = F],atemp[2, drop = F],rr10[1,2, drop = F])
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  
  grhorho[1] =grhorho[1, drop = F]+(grho%//%p[1, drop = F])  #/* adding contribution of rho12 from initial cdfbvn function */   
  gc = gc+(gcumulfrompc)  #/* adding contribution of all abscissa originating from the probability function */
  gc[1:2] = gc[1:2, drop = F]+((gw1%|%gw2)%//%p[1, drop = F])   #/* inserting gradient contribution of first two abscissae directly from the cdfbvn function */
  
  
  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}   


#/*****************************************************************************************************************************************************
#                                                         The OVBS Method
#*****************************************************************************************************************************************************/


cdfmvnaovbs = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  p = gss.zeros(m-2,1) 
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s        
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  for(h in seq.int(1,m-3,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)
    l = output[[1]]
    d = output[[2]]
    sigtemp1 = l[1:3,1:3, drop = F]%**%d[1:3,1:3, drop = F]%**%(t((l[1:3,1:3, drop = F])))
    p[h+1] = noncdftvn(mutemp[1:3, drop = F],sigtemp1,(t((atemp[(h+1):(h+3), drop = F]))))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp1[1:2,1:2, drop = F],(t((atemp[(h+1):(h+2), drop = F])))))
  }
  output = list(gss.prodc(p),s1)
  return(output)
}


pdfmvnaovbs = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  p = gss.zeros(m-2,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m) 
    s1 = s
  }
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  rr10=sigtemp
  p[1] =gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F])
  output = ldltblock(sigtemp,1)
  l = output[[1]]
  d = output[[2]]
  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  gc = gss.zeros(m,1)
  gcumulfrompc=gss.zeros(m,1)
  
  for(h in seq.int(1,m-3,1)) {
    output = univariatenormaltrunc(mutemp[1, drop = F],d[1,1, drop = F],atemp[h, drop = F])        
    mutilde = output[[1]]
    omega = output[[2]]
    if(h==1) {
      .condcovsigtrunc <<- 0
      .condcovmeantrunc <<- 0
    } else {
      .condcovsigtrunc <<- 1
      .condcovmeantrunc <<- 1
    }
    
    output = gcondmeantrunc(diag(m-h),mutemp,sigtemp,atemp[h, drop = F])
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcovtrunc(mutemp,sigtemp,atemp[h, drop = F])
    gmucov = output[[1]]
    gxcov = output[[2]]
    gccov = output[[3]]
    
    if(h==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-1,(m-1+((m-1)%**%(m)%//%2)))
    } else if(h > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[h,] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[2:(m-h+1), drop = F]+l[2:(m-h+1),1, drop = F]%**%(mutilde-mutemp[1, drop = F])
    output = ldltupspecial(l,d,omega,1)
    l = output[[1]]
    d = output[[2]]
    sigtemp = l%**%d%**%t(l) 
    
    p[h+1] = noncdftvn(mutemp[1:3, drop = F],sigtemp[1:3,1:3, drop = F],(t((atemp[(h+1):(h+3), drop = F]))))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],(t((atemp[(h+1):(h+2), drop = F])))))
    
    output = gradnoncdftvnbycdfbvn(mutemp[1:3, drop = F],sigtemp[1:3,1:3, drop = F],(t((atemp[(h+1):(h+3), drop = F]))))
    gfrompmu = output[[1]]
    gfrompcov = output[[2]]
    gfrompc = output[[3]]
    grhorho = grhorho+(1%//%p[h+1, drop = F])%**%((gcumulmusig[,1:3, drop = F]%**%gfrompmu)+((gcumulmusig[, (m-h+1):(m-h+3), drop = F]%~%gcumulmusig[, (2%**%m-2%**%h+1):(2%**%m-2%**%h+2), drop = F]%~%gcumulmusig[,3%**%m-3%**%h, drop = F])%**%gfrompcov))
    
    gcumulfrompc[(h+1):(h+3)] = gcumulfrompc[(h+1):(h+3), drop = F]+(gfrompc%//%p[h+1, drop = F])  #/* getting gradient contribution directly from noncdfbvn%//%noncdfn for absiccae except the first two */ 
    gc = gc+(1%//%p[h+1, drop = F])%**%(gcumulc[,1:3, drop = F]%**%gfrompmu+(gcumulc[, (m-h+1):(m-h+3), drop = F]%~%gcumulc[, (2%**%m-2%**%h+1):(2%**%m-2%**%h+2), drop = F]%~%gcumulc[,3%**%m-3%**%h, drop = F])%**%gfrompcov)
    
  }
  
  output = gradcdftvn(t((atemp[1:3, drop = F])),(rr10[1,2, drop = F]%|%rr10[1,3, drop = F]%|%rr10[2,3, drop = F]))
  gw = output[[1]]
  grho = output[[2]]
  
  grhorho[1:2] =grhorho[1:2, drop = F]+((grho[1:2, drop = F])%//%p[1, drop = F])  #/* adding contribution of rho12, rho13 from initial cdftvn function */  
  grhorho[m] =grhorho[m, drop = F]+(grho[3, drop = F]%//%p[1, drop = F]) #/* adding contribution of rho 23 from initial cdftvn function */
  gc = gc+(gcumulfrompc)  #/* adding contribution of all abscissae originating from the probability function */
  gc[1:3] = gc[1:3, drop = F]+((gw)%//%p[1, drop = F])   #/* inserting gradient contribution of first two absiccae directly from the cdfbvn function */

  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}   


#/*****************************************************************************************************************************************************
#                                                         The TGBME Method
#*****************************************************************************************************************************************************/

cdfmvnatgbme = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)   
    s1 = s
  }
  atemp = a[,temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  output = ldltblock(sigtemp,2)
  
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  p = gss.zeros(ktilde+1,1) 
  
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  
  for(k1 in seq.int(1,ktilde,1)) {
    output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],(t((atemp[(2%**%k1-1):(2%**%k1), drop = F]))))
    mutilde = output[[1]]
    omega = output[[2]]
    mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
    l=l[3:(m-2%**%k1+2),3:(m-2%**%k1+2), drop = F]
    d=d[3:(m-2%**%k1+2),3:(m-2%**%k1+2), drop = F]
    
    if(m>=2%**%k1+2) {
      sigtemp1 = d[1:2,1:2, drop = F]
      newtemp = 1%//%(sqrt(gss.diag(sigtemp1)))
      mutemp1 = newtemp%.*%((t((atemp[,(2%**%k1+1):(2%**%k1+2), drop = F])))-mutemp[1:2, drop = F])
      cortemp = cov2cor(sigtemp1) 
      p[k1+1] =gss.cdfbvn(mutemp1[1, drop = F],mutemp1[2, drop = F],cortemp[1,2, drop = F])  
    } else {
      p[k1+1] = pnorm((atemp[2%**%k+1, drop = F]-mutemp[1, drop = F])%//%sqrt(d[1,1, drop = F]))
    }
  }
  output = list(gss.prodc(p),s1)
  return(output)
}


pdfmvnatgbme = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  
  p = gss.zeros(m-1,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)   
    s1 = s
  }
  
  atemp = a[,temp1, drop = F]
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  rr10=sigtemp
  output = ldltblock(sigtemp,2)
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  p = gss.zeros(ktilde+1,1) 
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  
  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  
  gc = gss.zeros(m,1)
  gcumulfrompc=gss.zeros(m,1)
  
  for(k1 in seq.int(1,ktilde,1)) {
    output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],(t((atemp[(2%**%k1-1):(2%**%k1), drop = F]))))
    mutilde = output[[1]]
    omega = output[[2]]
    
    if(k1==1) {
      .condcov <<- 0
      .condcovmeantrunc <<- 0
    } else {
      .condcov <<- 1
      .condcovmeantrunc <<- 1
    }

    output = gcondmeantrunc(diag((m-2%**%k1)[,]),mutemp,sigtemp,(t((atemp[(2%**%k1-1):(2%**%k1), drop = F]))))
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcov(diag((m-2%**%k1)[,]),sigtemp)
    gy1 = output[[1]]
    gxcov = output[[2]]
    gmucov=gss.zeros(nrow(gmumean),((m-2%**%k1)%**%(m-2%**%k1+1)%//%2))

    gccov = gss.zeros(2,((m-2%**%k1)%**%(m-2%**%k1+1)%//%2))
    
    if(k1==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-2,(m-2+((m-2)%**%(m-1)%//%2)))
    } else if(k1 > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[(2%**%(k1-1)+1):(2%**%k1),] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
    l=l[3:(m-2%**%k1+2),3:(m-2%**%k1+2), drop = F]
    d=d[3:(m-2%**%k1+2),3:(m-2%**%k1+2), drop = F]
    sigtemp = l%**%d%**%t(l)
    
    if(m>=2%**%k1+2) {
      p[k1+1] = noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],(t((atemp[(2%**%k1+1):(2%**%k1+2), drop = F]))))
      output = gradnoncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],(t((atemp[(2%**%k1+1):(2%**%k1+2), drop = F]))))
      gfrompmu = output[[1]]
      gfrompcov = output[[2]]
      gfrompc = output[[3]]
      grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1:2, drop = F]%**%gfrompmu+(gcumulmusig[, (m-2%**%k1+1):(m-2%**%k1+2), drop = F]%~%gcumulmusig[,2%**%m-4%**%k1+1, drop = F])%**%gfrompcov)
      gcumulfrompc[(2%**%k1+1):(2%**%k1+2)] = gcumulfrompc[(2%**%k1+1):(2%**%k1+2), drop = F]+(gfrompc%//%p[k1+1, drop = F])  #/* getting gradient contribution directly from noncdfbvn for absiccae except the first two */       
      gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1:2, drop = F]%**%gfrompmu+(gcumulc[, (m-2%**%k1+1):(m-2%**%k1+2), drop = F]%~%gcumulc[,2%**%m-4%**%k1+1, drop = F])%**%gfrompcov)
      
    } else {
      p[k1+1] = noncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[2%**%k+1, drop = F])
      output = gradnoncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[2%**%k1+1, drop = F])
      gfrompmu = output[[1]]
      gfrompcov = output[[2]]
      gfrompc = output[[3]]
      grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1, drop = F]%**%gfrompmu+gcumulmusig[,m-2%**%k1+1, drop = F]%**%gfrompcov)
      gcumulfrompc[2%**%k1+1] = gcumulfrompc[2%**%k1+1, drop = F]+(gfrompc%//%p[k1+1, drop = F])
      gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1, drop = F]%**%gfrompmu+gcumulc[,m-2%**%k1+1, drop = F]%**%gfrompcov)           
      
    }
    
  }
  
  output = gradcdfbvn(atemp[1, drop = F],atemp[2, drop = F],rr10[1,2, drop = F])
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  
  grhorho[1] =grhorho[1, drop = F]+(grho%//%p[1, drop = F])  #/* adding contribution of rho12 from initial cdfbvn function */   
  gc = gc+(gcumulfrompc)  #/* adding contribution of all abscissa originating from the probability function */
  gc[1:2] = gc[1:2, drop = F]+((gw1%|%gw2)%//%p[1, drop = F])   #/* inserting gradient contribution of first two abscissae directly from the cdfbvn function */
  
  
  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}     

#/*****************************************************************************************************************************************************
#                                                         The BME Method
#*****************************************************************************************************************************************************/

cdfmvnabme = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s        
  }
  atemp = t((a[temp1, drop = F]))
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  output = ldltblock(sigtemp,2)
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  p = gss.zeros(ktilde+1,1) 
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  for(k1 in seq.int(1,ktilde,1)) {
    output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],atemp[(2%**%k1-1):(2%**%k1), drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
    output = ldltupspecial(l,d,omega,2)
    l = output[[1]]
    d = output[[2]]
    if(m>=2%**%k1+2) {
      p[k1+1] = noncdfbvn(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])
    } else {
      p[k1+1] = noncdfn(mutemp[1, drop = F],d[1,1, drop = F],atemp[2%**%k+1, drop = F])
    }
  }
  output = list(gss.prodc(p),s1)
  return(output)
}


pdfmvnabme = function(a,rr,s) { 
  
  a = matrix(a, 1)
  m=ncol(a)
  
  p = gss.zeros(m-1,1)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal ==1) {
    temp1 = ordering(a,rr)
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s        
  }
  
  atemp = t((a[temp1, drop = F]))
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  rr10=sigtemp
  output = ldltblock(sigtemp,2)
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  p = gss.zeros(ktilde+1,1) 
  p[1] =gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  
  grhorho = gss.zeros(m%**%(m-1)%//%2,1)
  gc = gss.zeros(m,1)
  gcumulfrompc=gss.zeros(m,1)
  
  for(k1 in seq.int(1,ktilde,1)) {
    output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],atemp[(2%**%k1-1):(2%**%k1), drop = F])
    mutilde = output[[1]]
    omega = output[[2]]
    
    if(k1==1) {
      .condcovsigtrunc<<-0
      .condcovmeantrunc<<-0
    } else {
      .condcovsigtrunc<<-1
      .condcovmeantrunc<<-1
    }
    
    output = gcondmeantrunc(gss.eye(m-2%**%k1),mutemp,sigtemp,atemp[(2%**%k1-1):(2%**%k1), drop = F])
    gy = output[[1]]
    gmumean = output[[2]]
    gxmean = output[[3]]
    gcmean = output[[4]]
    output = gcondcovtrunc(mutemp,sigtemp,atemp[(2%**%k1-1):(2%**%k1), drop = F])
    gmucov = output[[1]]
    gxcov = output[[2]]
    gccov = output[[3]]
    
    if(k1==1) {
      gcumulmusig = gxmean%~%gxcov
      gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-2,(m-2+((m-2)%**%(m-1)%//%2)))
    } else if(k1 > 1) {
      gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
      gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
      gcumulc[(2%**%(k1-1)+1):(2%**%k1),] = (gcmean%~%gccov)
    }
    
    mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
    output = ldltupspecial(l,d,omega,2)
    l = output[[1]]
    d = output[[2]]
    sigtemp = l%**%d%**%t(l)
    
    if(m>=2%**%k1+2) {
      p[k1+1] = noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])
      output = gradnoncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])
      gfrompmu = output[[1]]
      gfrompcov = output[[2]]
      gfrompc = output[[3]]
      grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1:2, drop = F]%**%gfrompmu+(gcumulmusig[, (m-2%**%k1+1):(m-2%**%k1+2), drop = F]%~%gcumulmusig[,2%**%m-4%**%k1+1, drop = F])%**%gfrompcov)
      gcumulfrompc[(2%**%k1+1):(2%**%k1+2)] = gcumulfrompc[(2%**%k1+1):(2%**%k1+2), drop = F]+(gfrompc%//%p[k1+1, drop = F])  #/* getting gradient contribution directly from noncdfbvn for absiccae except the first two */       
      gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1:2, drop = F]%**%gfrompmu+(gcumulc[, (m-2%**%k1+1):(m-2%**%k1+2), drop = F]%~%gcumulc[,2%**%m-4%**%k1+1, drop = F])%**%gfrompcov)
    } else {
      p[k1+1] = noncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[2%**%k+1, drop = F])
      output = gradnoncdfn(mutemp[1, drop = F],sigtemp[1,1, drop = F],atemp[2%**%k1+1, drop = F])
      gfrompmu = output[[1]]
      gfrompcov = output[[2]]
      gfrompc = output[[3]]
      grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1, drop = F]%**%gfrompmu+gcumulmusig[,m-2%**%k1+1, drop = F]%**%gfrompcov)
      gcumulfrompc[2%**%k1+1] = gcumulfrompc[2%**%k1+1, drop = F]+(gfrompc%//%p[k1+1, drop = F])
      gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1, drop = F]%**%gfrompmu+gcumulc[,m-2%**%k1+1, drop = F]%**%gfrompcov)           
    }
  }
  output = gradcdfbvn(atemp[1, drop = F],atemp[2, drop = F],rr10[1,2, drop = F])
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  
  grhorho[1] =grhorho[1, drop = F]+(grho%//%p[1, drop = F])  #/* adding contribution of rho12 from initial cdfbvn function */   
  gc = gc+(gcumulfrompc)  #/* adding contribution of all abscissa originating from the probability function */
  gc[1:2] = gc[1:2, drop = F]+((gw1%|%gw2)%//%p[1, drop = F])   #/* inserting gradient contribution of first two abscissae directly from the cdfbvn function */
  
  
  p = gss.prodc(p)
  gc=p%**%gc
  grhorho = p%**%grhorho
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  
  output = list(p,gc,grhorho,s1)
  return(output)  
}     

#/*****************************************************************************************************************************************************
#                                                         The TVBS Method
#*****************************************************************************************************************************************************/

cdfmvnatvbs = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)  
    s1 = s
  }
  atemp = t((a[temp1, drop = F]))
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  output = ldltblock(sigtemp,2)
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  if(m==2) {
    output = list(gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F]),s1)
    return(output) 
  } else if(m==3) {
    output = list(gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F]),s1)
    return(output) 
  } else if(m==4) {
    output = list(cdfqvn(atemp[1:4, drop = F],sigtemp[1:4,1:4, drop = F]),s1)
    return(output)
  } else {
    p = gss.zeros(ktilde,1) 
    p[1] =cdfqvn(atemp,sigtemp)
    for(k1 in seq.int(1,ktilde-1,1)) {
      output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],atemp[(2%**%k1-1):(2%**%k1), drop = F])
      mutilde = output[[1]]
      omega = output[[2]]
      mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
      output = ldltupspecial(l,d,omega,2)
      l = output[[1]]
      d = output[[2]]
      if(m>=2%**%k1+4) {
        sigtemp1 = l[1:4,1:4, drop = F]%**%d[1:4,1:4, drop = F]%**%(t((l[1:4,1:4, drop = F])))
        p[k1+1] = (noncdfqvn(mutemp[1:4, drop = F],sigtemp1,atemp[(2%**%k1+1):(2%**%k1+4), drop = F]))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp1[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])) 
      } else {   
        sigtemp1 = l%**%d%**%t(l)  
        p[k1+1] =(noncdftvn(mutemp[1:3, drop = F],sigtemp1[1:3,1:3, drop = F],atemp[(2%**%k1+1):(2%**%k1+3), drop = F]))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp1[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])) 
      }
    }
    output = list(gss.prodc(p),s1)
    return(output)
  }
}


pdfmvnatvbs = function(a,rr,s) { 
  
  a = matrix(a, 1)
  
  m=ncol(a)
  
  if(.optimal == 0) {
    temp1= vecindascending(t(a))
    s1 = s
  } else if(.optimal == 2) {
    output = randomordering(a,rr,s)
    temp1 = output[[1]]
    s1 = output[[2]]
  } else if(.optimal == 3) {
    temp1 = gss.seqa(1,1,m)
    s1 = s        
  }
  
  atemp = t((a[temp1, drop = F]))
  mutemp=gss.zeros(m,1) 
  sigtemp = rr[temp1,temp1, drop = F]
  rr10=sigtemp
  output = ldltblock(sigtemp,2)
  l = output[[1]]
  d = output[[2]]
  k = floor(m%//%2)
  if(m-2%**%k==0) {
    ktilde = k-1
  } else {
    ktilde = k
  }
  if(m==2) {
    output = (gradcdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])) 
    gwfin1 = output[[1]]
    gwfin2 = output[[2]]
    grhorho = output[[3]]
    gc=gwfin1%|%gwfin2
    p = gss.cdfbvn(atemp[1, drop = F],atemp[2, drop = F],sigtemp[1,2, drop = F])
  } else if(m==3) {
    output = (gradcdftvn(atemp[1:3, drop = F],(sigtemp[1,2, drop = F]%|%sigtemp[1,3, drop = F]%|%sigtemp[2,3, drop = F]))) 
    gc = output[[1]]
    grhorho = output[[2]]
    p=gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F])
  } else if(m==4) {
    output = gradcdfqvn(atemp[1:4, drop = F],sigtemp[1:4,1:4, drop = F])
    gc = output[[1]]
    grhorho = output[[2]]
    p=cdfqvn(atemp[1:4, drop = F],sigtemp[1:4,1:4, drop = F])        
  } else {    
    p = gss.zeros(ktilde,1) 
    p[1] =cdfqvn(atemp,sigtemp)
    
    grhorho = gss.zeros(m%**%(m-1)%//%2,1)
    gc = gss.zeros(m,1)
    gcumulfrompc=gss.zeros(m,1)
    
    for(k1 in seq.int(1,ktilde-1,1)) {
      output = bivariatenormaltrunc(mutemp[1:2, drop = F],d[1:2,1:2, drop = F],atemp[(2%**%k1-1):(2%**%k1), drop = F])
      mutilde = output[[1]]
      omega = output[[2]]
      
      if(k1==1) {
        .condcovsigtrunc<<-0
        .condcovmeantrunc<<-0
      } else {
        .condcovsigtrunc<<-1
        .condcovmeantrunc<<-1
      }
      
      output = gcondmeantrunc(gss.eye(m-2%**%k1),mutemp,sigtemp,atemp[(2%**%k1-1):(2%**%k1), drop = F])
      gy = output[[1]]
      gmumean = output[[2]]
      gxmean = output[[3]]
      gcmean = output[[4]]
      output = gcondcovtrunc(mutemp,sigtemp,atemp[(2%**%k1-1):(2%**%k1), drop = F])
      gmucov = output[[1]]
      gxcov = output[[2]]
      gccov = output[[3]]
      
      if(k1==1) {
        gcumulmusig = gxmean%~%gxcov
        gcumulc = (gcmean%~%gccov)%|%gss.zeros(m-2,(m-2+((m-2)%**%(m-1)%//%2)))
      } else if(k1 > 1) {
        gcumulmusig = gcumulmusig%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov)) 
        gcumulc = gcumulc%**%((gmumean%~%gmucov)%|%(gxmean%~%gxcov))     
        gcumulc[(2%**%(k1-1)+1):(2%**%k1),] = (gcmean%~%gccov)
      }
      
      mutemp = mutemp[3:(m-2%**%k1+2), drop = F]+l[3:(m-2%**%k1+2),1:2, drop = F]%**%(mutilde-mutemp[1:2, drop = F])
      output = ldltupspecial(l,d,omega,2)
      l = output[[1]]
      d = output[[2]]
      sigtemp = l%**%d%**%t(l)
      if(m>=2%**%k1+4           ) {
        p[k1+1] = (noncdfqvn(mutemp[1:4, drop = F],sigtemp[1:4,1:4, drop = F],atemp[(2%**%k1+1):(2%**%k1+4), drop = F]))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])) 
        output = gradnoncdfqvnbycdfbvn(mutemp[1:4, drop = F],sigtemp[1:4,1:4, drop = F],atemp[(2%**%k1+1):(2%**%k1+4), drop = F])
        gfrompmu = output[[1]]
        gfrompcov = output[[2]]
        gfrompc = output[[3]]
        grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1:4, drop = F]%**%gfrompmu+(gcumulmusig[, (m-2%**%k1+1):(m-2%**%k1+4), drop = F]%~%gcumulmusig[, (2%**%m-4%**%k1+1):(2%**%m-4%**%k1+3), drop = F]%~%gcumulmusig[, (3%**%m-6%**%k1):(3%**%m-6%**%k1+1), drop = F]%~%gcumulmusig[,4%**%m-8%**%k1-2, drop = F])%**%gfrompcov)          
        gcumulfrompc[(2%**%k1+1):(2%**%k1+4)] = gcumulfrompc[(2%**%k1+1):(2%**%k1+4), drop = F]+(gfrompc%//%p[k1+1, drop = F])  #/* getting gradient contribution directly from noncdfbvn for absiccae except the first two */       
        gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1:4, drop = F]%**%gfrompmu+(gcumulc[, (m-2%**%k1+1):(m-2%**%k1+4), drop = F]%~%gcumulc[, (2%**%m-4%**%k1+1):(2%**%m-4%**%k1+3), drop = F]%~%gcumulc[, (3%**%m-6%**%k1):(3%**%m-6%**%k1+1), drop = F]%~%gcumulc[,4%**%m-8%**%k1-2, drop = F])%**%gfrompcov)
      } else {
        p[k1+1] =(noncdftvn(mutemp[1:3, drop = F],sigtemp[1:3,1:3, drop = F],atemp[(2%**%k1+1):(2%**%k1+3), drop = F]))%//%(noncdfbvn(mutemp[1:2, drop = F],sigtemp[1:2,1:2, drop = F],atemp[(2%**%k1+1):(2%**%k1+2), drop = F])) 
        output = gradnoncdftvnbycdfbvn(mutemp[1:3, drop = F],sigtemp[1:3,1:3, drop = F],atemp[(2%**%k1+1):(2%**%k1+3), drop = F])
        gfrompmu = output[[1]]
        gfrompcov = output[[2]]
        gfrompc = output[[3]]
        grhorho = grhorho+(1%//%p[k1+1, drop = F])%**%(gcumulmusig[,1:3, drop = F]%**%gfrompmu+(gcumulmusig[, (m-2%**%k1+1):(m-2%**%k1+3), drop = F]%~%gcumulmusig[, (2%**%m-4%**%k1+1):(2%**%m-4%**%k1+2), drop = F]%~%gcumulmusig[,3%**%m-6%**%k1, drop = F])%**%gfrompcov)
        gcumulfrompc[(2%**%k1+1):(2%**%k1+3)] = gcumulfrompc[(2%**%k1+1):(2%**%k1+3), drop = F]+(gfrompc%//%p[k1+1, drop = F])
        gc = gc+(1%//%p[k1+1, drop = F])%**%(gcumulc[,1:3, drop = F]%**%gfrompmu+(gcumulc[, (m-2%**%k1+1):(m-2%**%k1+3), drop = F]%~%gcumulc[, (2%**%m-4%**%k1+1):(2%**%m-4%**%k1+2), drop = F]%~%gcumulc[,3%**%m-6%**%k1, drop = F])%**%gfrompcov)           
      }
    }
    
    output = gradcdfqvn(atemp[1:4, drop = F],rr10[1:4,1:4, drop = F])
    gw = output[[1]]
    grho = output[[2]]
    
    grhorho[1:3] =grhorho[1:3, drop = F]+(grho[1:3, drop = F]%//%p[1, drop = F])  #/* adding contribution of rho12,rho13, and rho14 from initial cdfqvn function */  
    grhorho[m:(m+1)] =grhorho[m:(m+1), drop = F]+(grho[4:5, drop = F]%//%p[1, drop = F])  #/* adding contribution of rho23 & rho24 from initial cdfqvn function */
    grhorho[2%**%m-2] = grhorho[2%**%m-2, drop = F]+(grho[6, drop = F]%//%p[1, drop = F])   #/* adding contribution of rho34 from initial cdfqvn function */
    gc = gc+(gcumulfrompc)  #/* adding contribution of all abscissa originating from the probability function */
    gc[1:4] = gc[1:4, drop = F]+((gw)%//%p[1, drop = F])   #/* inserting gradient contribution of first four abscissae directly from the cdfqvn function */
    
    p = gss.prodc(p)
    gc=p%**%gc
    grhorho = p%**%grhorho
  }
  
  tempnew = gss.seqa(1,1,m)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  output = list(p,gc,grhorho,s1)
  return(output)  
}       

#/*****************************************************************************************************************************************************
#                                                         The SSJ Method
#*****************************************************************************************************************************************************/

cdfmvnassj = function(a,r,s) {
  
  a = matrix(a, 1)
  
  m = ncol(a)
  w = {}
  for(z in seq.int( 1, .perms, 1)) {
    output = gss.rndu(1,m,s)  
    d = output[[1]]
    s = output[[2]]
    aa = gss.seqa(1,1,m)%~%t(d)
    aa = gss.sortc(aa,2)
    w = w%|%(t((aa[,1, drop = F])))
  }
  s1=s
  w=t(w)
  n1 = ncol(w)
  p = 0
  j = 1
  
  a = a%.*%(a <5.7)+5.7%**%(a>=5.7)       
  ab = t(a)
  while(!(j > n1)) {
    x = ab[w[,j, drop = F],, drop = F]                  # 4x1 Matrix
    rho = gss.submat(r,w[,j, drop = F],w[,j, drop = F]) # 4x4 matrix
    y = gss.reshape((t((gss.ones(m,1)%x%x))),m,m)       # 4x4 Matrix
    z = gss.zeros(m, m)
    for (i in 1:m) {
      for (k in 1:m) {
        z[i, k] = gss.cdfbvn(x[i], y[i, k], rho[i, k])
      }
    }
    # z = gss.cdfbvn(x,y,rho)
    z = gss.diagrv(z,pnorm(x))
    z1 = pnorm(x)%.*%pnorm(y)
    z3 = pnorm(x)
    z2 = 1-z3
    cond = 1
    k = 3
    while(!(k > m)) {
      omega21 = z[k,1:(k-1), drop = F]-z1[k,1:(k-1), drop = F]
      omega11 = z[1:(k-1),1:(k-1), drop = F]-z1[1:(k-1),1:(k-1), drop = F]
      cm=solve(omega11)
      condk = z3[k, drop = F]+omega21%**%(cm)%**%z2[1:(k-1), drop = F]   
      cond=cond%**%condk
      k = k+1
    }
    pcomb = z[1,2, drop = F]%**%cond
    p = p+pcomb  
    j=j+1
  }
  
  if((p%//%n1) >0 & (p%//%n1) < 1) {
    output = list(p%//%n1,(s1))
    return(output)
  } else {
    print("yessss"    )
    output = list(cdfmvname(a,r),s1)
    return(output)
  }    
}    


pdfmvnassj = function(a,r,s) {
  
  a = matrix(a, 1)
  
  m = ncol(a)
  w = {}
  for(z in seq.int( 1, .perms, 1)) {
    output = gss.rndu(1,m,s)              
    d = output[[1]]
    s = output[[2]]
    aa = gss.seqa(1,1,m)%~%t(d)
    aa = gss.sortc(aa,2)
    w = w%|%(t((aa[,1, drop = F])))
  }
  s1=s
  w=t(w)
  
  n1 = ncol(w)
  p = 0
  j = 1
  a = a%.*%(a <5.7)+5.7%**%(a>=5.7)      
  ab = t(a)
  condpass={}
  while(!(j > n1)) {
    x = ab[w[,j, drop = F],, drop = F]                  # 4x1 Matrix
    rho = gss.submat(r,w[,j, drop = F],w[,j, drop = F]) # 4x4 matrix
    y = gss.reshape((t((gss.ones(m,1)%x%x))),m,m)       # 4x4 Matrix
    z = gss.zeros(m, m)
    for (i in 1:m) {
      for (k in 1:m) {
        z[i, k] = gss.cdfbvn(x[i], y[i, k], rho[i, k])
      }
    }
    # z = gss.cdfbvn(x,y,rho)
    z = gss.diagrv(z,pnorm(x))
    z1 = pnorm(x)%.*%pnorm(y)
    z3 = pnorm(x)
    z2 = 1-z3
    cond = 1
    k = 3
    while(!(k > m)) {
      omega21 = z[k,1:(k-1), drop = F]-z1[k,1:(k-1), drop = F]
      omega11 = z[1:(k-1),1:(k-1), drop = F]-z1[1:(k-1),1:(k-1), drop = F]
      cm=solve(omega11)
      condk = z3[k, drop = F]+omega21%**%(cm)%**%z2[1:(k-1), drop = F]  
      cond=cond%**%condk
      k = k+1
    }
    pcomb = z[1,2, drop = F]%**%cond
    condpass = condpass%~%cond
    p = p+pcomb    
    j=j+1
  }
  
  gwfinal = gss.zeros(1,m)
  grfinal = gss.zeros(1,((m-1)%**%(m)%//%2))
  j=1
  while(!(j > n1)) {
    x = ab[w[,j, drop = F], drop = F]
    rho = gss.submat(r,w[,j, drop = F],w[,j, drop = F])
    rhovec = {}
    gss.c=1
    while(!(gss.c==ncol(rho))) {
      rhovec = rhovec%~%rho[gss.c, (gss.c+1):ncol(rho), drop = F]
      gss.c=gss.c+1
    }
    y = gss.reshape((t((gss.ones(m,1)%x%x))),m,m)
    z = gss.zeros(m, m)
    for (i in 1:m) {
      for (k in 1:m) {
        z[i, k] = gss.cdfbvn(x[i], y[i, k], rho[i, k])
      }
    }
    z = gss.diagrv(z,pnorm(x))
    z1 = pnorm(x)%.*%pnorm(y)
    z3 = pnorm(x)
    z2 = 1-z3
    rho1 = gss.diagrv(rho,gss.zeros(m,1))
    rho2 = sqrt(1-rho1^2)
    g3 = dnorm(x)
    g5 = (y-rho1%.*%x)%//%rho2
    g10 = g3%.*%pnorm(g5)
    g11 = g3%.*%pnorm(y)
    g12 = g10-g11
    g13 = g12
    g14 = g3%.*%(1-2%**%pnorm(x))
    g15 = gss.diagrv(g13,g14)
    g20 = -g3
    
    g25 = (1%//%rho2)%.*%g3%.*%dnorm(g5)
    g25 = gss.diagrv(g25,gss.zeros(m,1))
    
    
    g30 = g10[1,2, drop = F]%~%g10[2,1, drop = F]%~%gss.zeros(1,m-2)
    g53 = g25[1,2, drop = F]%~%gss.zeros(1,ncol(rhovec)-1)
    
    
    k = 3
    gw1 = gss.zeros(1,m)
    gr1 = gss.zeros(1,ncol(rhovec))
    while(!(k > m)) {
      omega21 = z[k,1:(k-1), drop = F]-z1[k,1:(k-1), drop = F]
      omega11 = z[1:(k-1),1:(k-1), drop = F]-z1[1:(k-1),1:(k-1), drop = F]
      invomg11=solve(omega11)
      condk = z3[k, drop = F]+omega21%**%invomg11%**%z2[1:(k-1), drop = F]
      g31 = gss.zeros(1,m)
      g31[k] = g3[k, drop = F]
      l = 1
      g40={}
      g46={}
      g51={}
      g81 = z2[1:(k-1), drop = F]    
      while(!(l == k)) {
        g35 = gss.zeros(k-1,k-1)
        g35[l,1:(k-1)] = g15[l,1:(k-1), drop = F]
        g35=g35+t(g35)
        g35 = gss.diagrv(g35,(gss.diag(g35)%//%2))
        g36 = -invomg11%**%g35%**%invomg11
        g36 = omega21%**%g36%**%g81
        g40=g40%~%g36
        
        g45 = gss.zeros(1,k-1)
        g45[l] = g15[l,k, drop = F]
        g46 = g46%~%(g45%**%invomg11%**%g81)
        
        g50 = gss.zeros(k-1,1)
        g50[l] =g20[l, drop = F]
        g51 = g51%~%(omega21%**%invomg11%**%g50)
        l=l+1
      }
      
      g40 = g40%~%gss.zeros(1,m-(k-1))
      
      g49 = g15[k,1:(k-1), drop = F]
      g46 = g46%~%(g49%**%invomg11%**%g81)
      
      g47 = gss.zeros(1,m)
      g47[1:ncol(g46)] = g46
      
      g51 =g51%~%gss.zeros(1,m-(k-1))
      
      gw = (g31+g40+g47+g51)%**%((condpass[j, drop = F])%//%condk)%.*%z[1,2, drop = F]
      
      #/* start here for gradients with respect to rho parameters */
      
      l=1
      kk = ncol(rhovec)
      g60={}
      g65={}
      while(!(l>kk)) {
        g55 = gss.zeros(1,kk)
        g55[l] =1 
        mm = 1
        sk=0
        g56 = gss.zeros(m,m)
        while(!(mm>m-1)) {
          g56[mm, (mm+1):m] = (g55[(sk+1):(sk+m-mm), drop = F])
          sk=sk+m-mm
          mm=mm+1
        }
        g57=g56+t(g56)
        g59 = g57%.*%g25
        g59 = g59[1:(k-1),1:(k-1), drop = F]
        g58 = -invomg11%**%g59%**%invomg11
        g58 = omega21%**%g58%**%z2[1:(k-1), drop = F]
        g60 = g60%~%g58
        
        g61 = g57[k,1:(k-1), drop = F]%.*%g25[k,1:(k-1), drop = F]
        g62 = g61%**%invomg11%**%z2[1:(k-1), drop = F]
        g65 = g65%~%g62      
        l=l+1
      }
      
      gr = (g60+g65)%**%((condpass[j, drop = F])%//%condk)%.*%z[1,2, drop = F]
      gw1 = gw+gw1
      gr1 = gr+gr1
      k=k+1
    }
    gw2 = (g30%.*%condpass[j, drop = F])+gw1
    gr2 = (g53%.*%condpass[j, drop = F])+gr1
    
    mm=1
    sk=0
    gr3 = gss.zeros(m,m)
    while(!(mm>m-1)) {
      gr3[mm, (mm+1):m] = (gr2[(sk+1):(sk+m-mm), drop = F])
      sk=sk+m-mm
      mm=mm+1
    }
    
    #/* commands below to resequence gradients based on permutation */
    
    aaa = match(gss.seqa(1,1,m),w[,j, drop = F])
    gwf = gw2[aaa, drop = F]
    grf = gss.submat(gr3,aaa,aaa) 
    
    res = {}
    ir = 1
    while(!(ir == m)) {
      jr = ir+1
      res = res%|%((gss.ones(m-ir,1)%x%w[ir,j, drop = F])%~%w[jr:m,j, drop = F])
      ir=ir+1
    }
    
    res=t(res)
    res2={}
    gss.t=1
    while(!(gss.t>ncol(res))) {
      res2 = res2%~%gss.sortc(res[,gss.t, drop = F],1)
      gss.t=gss.t+1
    }
    
    res2=t(res2)
    res1 = t(combn(m,2))
    res1 = res1[,1, drop = F] %.*%(10^(trunc (log(res1[,2, drop = F])%//%log(10))+1)) + res1[,2, drop = F]
    res2 = res2[,1, drop = F] %.*%(10^(trunc (log(res2[,2, drop = F])%//%log(10))+1)) + res2[,2, drop = F]
    
    aaa1 = match(res1,res2)
    grf = gr2[aaa1, drop = F]
    gwfinal = gwfinal+gwf
    grfinal = grfinal+grf
    j=j+1
  }
  
  
  if((p%//%n1) >0 & p%//%n1 < 1) {
    output = list((p%//%n1),(t((gwfinal%//%n1))),(t((grfinal%//%n1))),s1)
    return(output)
  } else {
    print("yessssp"    )
    output = pdfmvname(a,r)
    gwfinal = output[[1]]
    grfinal = output[[2]]
    output = list(cdfmvname(a,r),gwfinal,grfinal,s1)
    return(output)
  }    
}


#/*****************************************************************************************************************************************************
#                            Re-arranging abscissae based on ascending order of (atemp-mutemp)/sqrt(D[1,1])
#*****************************************************************************************************************************************************/ 

ordering = function(a,rr) { 
  
  a = matrix(a, nrow = 1)
  
  m=ncol(a)
  for(h in seq.int(1,m-2,1)) {
    if(h==1) {
      temp4 = vecindascending(t(a))
      tempper = temp4
      mutemp=gss.zeros(m,1) 
      sigtemp = rr[temp4,temp4, drop = F]
      atemp = t((a[temp4, drop = F]))
    }
    output = multrunc(mutemp,sigtemp,atemp[1, drop = F])
    mutemp = output[[1]]
    sigtemp = output[[2]]
    mutemp = mutemp[2:(m-h+1), drop = F]
    sigtemp = sigtemp[2:(m-h+1),2:(m-h+1), drop = F]
    atemp=atemp[2:(m-h+1), drop = F]
    temp2 = tempper[(h+1):m, drop = F]
    temp3 = vecindascending((atemp-mutemp)%//%(sqrt(gss.diag(sigtemp))))
    tempper[(h+1):m] =temp2[temp3, drop = F]
  }
  return(tempper)
}

#/*****************************************************************************************************************************************************
#                                          Re-arranging abscissae based on random ordering
#*****************************************************************************************************************************************************/ 

randomordering = function(a,rr,s) { 
  
  a = matrix(a, nrow = 1)
  
  m=ncol(a)
  output = gss.rndu(1,m,s)              
  d = output[[1]]
  s1 = output[[2]]
  aa = gss.seqa(1,1,m)%~%t(d)
  aa = gss.sortc(aa,2)
  output = list(aa[,1, drop = F],s1)
  return(output)
}

#/*********************************************************************************************************************************************/
#/* All the procedures below have been written by Bhat to facilitate matrix manipulations and develop gradients for a variety of functions */
#/* These are all auxiliary procedures called by the main procedures above */
#/*********************************************************************************************************************************************/
#
#/************************************************************************************************************************************************
#Next few procedures deal with mean and variance of truncated normal, and mean and variance of untruncated elements of a multivariate normal given
#some other elements are truncated; this set of procedures also is useful to compute the gradients of the mean and variance of the truncated elements
#with respect to the unconditional mean and variance, and the gradients of the mean and variance of the untruncated elements (given other elements are truncated)
#with respect to the unconditional mean and variance; other auxiliary gradient procedures needed for the gradient computations are also included in 
#this set of procedures.*/
#
#/*****************************************************************************************************************************************************
#                SET 2: THE CODES BELOW CONSTITUTE THE ROUTINES FOR TRUNCATION-BASED MEANS/VARIANCE AND THEIR GRADIENTS
#*****************************************************************************************************************************************************/


univariatenormaltrunc = function(muuntrunc,siguntrunc,trpoint) { 
  sig = sqrt(siguntrunc)
  w = (trpoint-muuntrunc)%//%sig
  lam = (-dnorm(w))%//%pnorm(w)
  mutrunc = muuntrunc+(sig%.*%lam)
  sigtrunc = siguntrunc%.*%(1+(lam%.*%(w-lam)))  
  output = list(mutrunc,sigtrunc)
  return(output)  
}


bivariatenormaltrunc = function(muuntrunc,cov,trpoint) { 
  
  muuntrunc = matrix(muuntrunc)
  trpoint = matrix(trpoint)
  
  m = nrow(trpoint) 
  newtrpoint = (trpoint-muuntrunc)%//%(sqrt(gss.diag(cov)))   
  cor = cov2cor(cov)
  rho = cor[1,2, drop = F]
  p = gss.cdfbvn(newtrpoint[1, drop = F],newtrpoint[2, drop = F],rho)
  newtrpointrev = gss.rev(newtrpoint)
  rhotilde = sqrt((1-rho^2))
  tr1 = (newtrpoint-rho%**%newtrpointrev)%//%rhotilde
  tr2 = gss.rev(tr1)
  pd1 = dnorm(newtrpoint)
  pd2 = gss.rev(pd1)
  cd1 = pnorm(tr1)
  cd2 = gss.rev(cd1)
  mu = (-(rho%**%pd2)%.*%cd1-pd1%.*%cd2)%//%p 
  sig = ((p-(newtrpoint%.*%pd1%.*%cd2)-((rho^2)%**%newtrpointrev%.*%pd2%.*%cd1)+(rhotilde%**%rho%**%pd2%.*%dnorm(tr1)))%//%p)-mu^2
  sig12 = ((rho%**%p-rho%**%newtrpoint[1, drop = F]%**%pd1[1, drop = F]%**%cd2[1, drop = F]+rhotilde%**%pd1[1, drop = F]%**%dnorm(tr2[1, drop = F])-rho%**%newtrpointrev[1, drop = F]%**%pd2[1, drop = F]%**%cd1[1, drop = F])%//%p)-gss.prodc(mu)
  diagcovroot = sqrt(gss.diag(cov))
  diag1 = gss.diagrv(gss.zeros(m,m),diagcovroot)
  mu = muuntrunc+(diagcovroot%.*%mu)  
  omg = gss.zeros(2,2)
  omg[1,2] = sig12
  omg[2,1] =sig12
  omg = gss.diagrv(omg,sig)  
  omg = diag1%**%omg%**%diag1  
  output = list(mu,omg)
  return(output)
}  


gradunivariatenormaltrunc = function(muuntrunc,siguntrunc,trpoint) { 
  sig = sqrt(siguntrunc)
  w = (trpoint-muuntrunc)%//%sig
  lam = (-dnorm(w))%//%pnorm(w)
  mutrunc = muuntrunc+(sig%.*%lam)
  sigtrunc = siguntrunc%.*%(1+(lam%.*%(w-lam)))
  dlamdw = lam%.*%(lam-w)
  dmutrunc = (t((1-dlamdw)))%|%((t(((lam-dlamdw%.*%w)%//%(2%**%sig)))))%|%(t((dlamdw)))
  dsigtemp = dlamdw%.*%(w-2%**%lam)
  dsigtruncdw = (sig^2)%.*%(lam+dsigtemp)
  dsigtemp1 = (1%//%sig)%.*%dsigtruncdw
  dsigtemp2 = 1+(lam%.*%(w-lam))-((w%//%2)%.*%(lam+dsigtemp))
  dsigtrunc = (t((-dsigtemp1)))%|%(t(dsigtemp2))%|%(t(dsigtemp1))     
  output = list(dmutrunc,dsigtrunc)
  return(output)  
}    


gradbivariatenormaltrunc = function(muuntrunc,cov,trpoint) { 
  
  muuntrunc = matrix(muuntrunc)
  trpoint = matrix(trpoint)
  
  m = nrow(trpoint) 
  sigorig = sqrt(gss.diag(cov))
  sigorigsq = sigorig^2
  sigoriginv = 1%//%sigorig
  sigoriginvrev = gss.rev(sigoriginv)
  newtrpoint = (trpoint-muuntrunc)%//%sigorig   
  cor = cov2cor(cov)
  rho = cor[1,2, drop = F]
  p = gss.cdfbvn(newtrpoint[1, drop = F],newtrpoint[2, drop = F],rho)
  newtrpointrev = gss.rev(newtrpoint)
  rhotilde = sqrt((1-rho^2))
  rhotilde2 = rhotilde^2
  tr1 = (newtrpoint-rho%**%newtrpointrev)%//%rhotilde
  tr2 = gss.rev(tr1)
  trcomp = rho%**%newtrpoint-newtrpointrev
  trcomprev = gss.rev(trcomp)
  pd1 = dnorm(newtrpoint)
  pd2 = gss.rev(pd1)
  cd1 = pnorm(tr1)
  cd2 = gss.rev(cd1)
  del1 = pd1%.*%cd2
  del2 = pd2%.*%cd1
  pdf2 = (1%//%rhotilde)%**%(pd1[1, drop = F]%**%dnorm(tr1[2, drop = F]))
  mu = (-rho%**%del2-del1)%//%p
  
  dmutilde1dw1=(1%//%p)%**%sigorig%.*%((newtrpoint%.*%del1)+(del1%.*%(-mu)))
  dmutilde1dw2 = (1%//%p)%**%(-sigorig)%.*%((rhotilde2%**%pdf2)[,]+(mu-rho%**%(newtrpointrev))%.*%(del2))
  dmutildedrho = (1%//%p)%**%(-sigorig)%.*%(del2+pdf2%.*%(mu-newtrpoint))
  dmutilde1da1 = 1-(sigoriginv)%.*%dmutilde1dw1
  dmutilde1da2 = -(sigoriginvrev)%.*%dmutilde1dw2
  dmutilde1dsig1 = -(sigoriginv%.*%(newtrpoint)%.*%dmutilde1dw1-mu+rho%**%sigoriginv%.*%dmutildedrho)
  dmutilde1dsig2 = -(sigoriginvrev%.*%(newtrpointrev)%.*%dmutilde1dw2+rho%**%sigoriginvrev%.*%dmutildedrho)
  dmutilde1dsig1sq=0.5%**%(sigoriginv)%.*%dmutilde1dsig1
  dmutilde1dsig2sq=0.5%**%(sigoriginvrev)%.*%dmutilde1dsig2  
  dmutildedsig12 = (gss.prodc(sigoriginv))%.*%dmutildedrho 
  dmutilde1dtr1 = (sigoriginv)%.*%dmutilde1dw1
  dmutilde1dtr2 = (sigoriginvrev)%.*%dmutilde1dw2
  
  dmuderiv1 = dmutilde1da1[1, drop = F]%|%dmutilde1da2[1, drop = F]%|%dmutilde1dsig1sq[1, drop = F]%|%dmutildedsig12[1, drop = F]%|%dmutilde1dsig2sq[1, drop = F]%|%dmutilde1dtr1[1, drop = F]%|%dmutilde1dtr2[1, drop = F]
  dmuderiv2 = dmutilde1da2[2, drop = F]%|%dmutilde1da1[2, drop = F]%|%dmutilde1dsig2sq[2, drop = F]%|%dmutildedsig12[2, drop = F]%|%dmutilde1dsig1sq[2, drop = F]%|%dmutilde1dtr2[2, drop = F]%|%dmutilde1dtr1[2, drop = F]  
  dmuderiv = dmuderiv1%~%dmuderiv2
  
  
  sig = ((p-(newtrpoint%.*%pd1%.*%cd2)-((rho^2)%**%newtrpointrev%.*%pd2%.*%cd1)+(rhotilde%**%rho%**%pd2%.*%dnorm(tr1)))%//%p)-mu^2
  sig12 = ((rho%**%p-rho%**%newtrpoint[1, drop = F]%**%pd1[1, drop = F]%**%cd2[1, drop = F]+rhotilde%**%pd1[1, drop = F]%**%dnorm(tr2[1, drop = F])-rho%**%newtrpointrev[1, drop = F]%**%pd2[1, drop = F]%**%cd1[1, drop = F])%//%p)-gss.prodc(mu)
  diagcovroot = sqrt(gss.diag(cov))
  diag1 = gss.diagrv(gss.zeros(m,m),diagcovroot)
  
  
  dsigtilde1dw1 = -(del1%//%p)%.*%((-newtrpoint^2)+sig+mu^2)-2%**%mu%.*%(dmutilde1dw1%//%sigorig)
  dsigtilde1dw2 = -(1%//%p)%**%((rhotilde2%**%(rho%**%pdf2%**%newtrpointrev+newtrpoint%**%pdf2-del2))-(rho^2)%**%(newtrpointrev^2)%.*%del2+(sig+mu^2)%.*%del2)-2%**%mu%.*%(dmutilde1dw2%//%sigorig)
  dsigtildedrho = -(1%//%p)%**%(pdf2%**%(-newtrpoint^2)+2%**%rho%**%newtrpointrev%.*%del2-(2%**%pdf2%**%rhotilde2)[,]+((sig+mu^2)%**%pdf2))-2%**%mu%.*%(dmutildedrho%//%sigorig)
  dsigtilde = (1%//%p)%**%((rho%**%newtrpoint)%.*%((rho%**%pdf2)[,]+newtrpoint%.*%del1)-pdf2%**%newtrpoint-del1%**%(sig12+gss.prodc(mu)))
  dsigtilde12dw1 = dsigtilde[1, drop = F]-(gss.sumc((mu%.*%((dmutilde1dw2[2, drop = F]%|%dmutilde1dw1[1, drop = F])%//%(gss.rev(sigorig)))))) 
  dsigtilde12dw2 = dsigtilde[2, drop = F]-(gss.sumc((mu%.*%((dmutilde1dw1[2, drop = F]%|%dmutilde1dw2[1, drop = F])%//%(gss.rev(sigorig))))))    
  dsigtilde12drho = 1-(1%//%p)%**%((-gss.prodc(newtrpoint))%**%pdf2+(gss.sumc((newtrpoint%.*%del1)))+(gss.prodc(mu)+sig12)%**%pdf2)  
  dsigtilde12drho = dsigtilde12drho-gss.sumc((mu%.*%(gss.rev(dmutildedrho))%.*%(sigoriginvrev))) 
  
  domg1dtr1 = (sigorig)%.*%dsigtilde1dw1
  domg12dtr = (gss.rev(sigorig))%.*%(dsigtilde12dw1%|%dsigtilde12dw2)
  domg1dtr2 = ((sigorig)^2)%.*%sigoriginvrev%.*%dsigtilde1dw2
  domg1da1 = - domg1dtr1
  domg12da=-domg12dtr
  domg1da2 = - domg1dtr2
  domg1s1 = sig+dsigtilde1dw1%.*%(newtrpoint)%.*%(-0.5)+dsigtildedrho%.*%rho%**%(-0.5)
  domg1s12 = (sigorig%//%(gss.rev(sigorig)))%.*%dsigtildedrho
  domg1s2 = (-0.5)%**%(sigorigsq%.*%(1%//%gss.rev(sigorigsq)))%.*%(dsigtilde1dw2%.*%newtrpointrev+dsigtildedrho%**%rho)
  domg12s1 = (0.5)%**%(gss.rev(sigorig)%//%sigorig)%.*%(sig12[,]-newtrpoint%.*%(dsigtilde12dw1%|%dsigtilde12dw2)-(rho%**%dsigtilde12drho)[,])  
  domg12s12 = dsigtilde12drho
  
  domgderiv11 = domg1da1[1, drop = F]%|%domg1da2[1, drop = F]%|%domg1s1[1, drop = F]%|%domg1s12[1, drop = F]%|%domg1s2[1, drop = F]%|%domg1dtr1[1, drop = F]%|%domg1dtr2[1, drop = F]
  domgderiv12 = (domg12da)%|%domg12s1[1, drop = F]%|%domg12s12%|%domg12s1[2, drop = F]%|%(domg12dtr)
  domgderiv22 = domg1da2[2, drop = F]%|%domg1da1[2, drop = F]%|%domg1s2[2, drop = F]%|%domg1s12[2, drop = F]%|%domg1s1[2, drop = F]%|%domg1dtr2[2, drop = F]%|%domg1dtr1[2, drop = F]
  domgderiv = domgderiv11%~%domgderiv12%~%domgderiv22
  output = list(dmuderiv,domgderiv)
  return(output)
}  


gcondmeantrunc = function(y,mu,x,gss.c) { 
  
  if (length(y) == 1) y = matrix(y)
  mu = matrix(mu)
  gss.c = matrix(gss.c)
  
  dim1 = nrow(y)
  dim2 = nrow(x)
  dimdiff = dim2-dim1
  
  ddcov = (dimdiff%**%(dimdiff+1))%//%2   #/* these next four lines are to get a way of sequencing the gradients later to put in the right order of elements in matrix x */
  ddcor = (dimdiff%**%(dimdiff-1))%//%2
  dd1cov = (dim1%**%(dim1+1))%//%2
  dd1cor = (dim1%**%(dim1-1))%//%2
  x11 = x[1:dimdiff,1:dimdiff, drop = F]
  x12 = x[(dimdiff+1):nrow(x),1:dimdiff, drop = F]
  x22 = x[(dimdiff+1):nrow(x), (dimdiff+1):nrow(x), drop = F]
  invx11 = solve(x11)
  
  gmutilde = (t((y%**%x12%**%invx11)))
  gmu = -gmutilde%|%diag(dim1)
  
  if(dimdiff==1) {
    output = univariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    mutrunc = output[[1]]
    sigtrunc = output[[2]]
    output = gradunivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    gmutrunc = output[[1]]
    gsigtrunc = output[[2]]
  } else if(dimdiff==2) {
    output = bivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)        
    mutrunc = output[[1]]
    sigtrunc = output[[2]]
    output = gradbivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    gmutrunc = output[[1]]
    gsigtrunc = output[[2]]
  }
  
  gy1 = x12%**%invx11%**%(mutrunc-mu[1:dimdiff, drop = F])
  ggy = gss.diagrv(gss.zeros(dim1,dim1),gy1)
  
  gmutildenew = gmutrunc%**%gmutilde
  gmu[1:dimdiff,] = (gmu[1:dimdiff,, drop = F])+gmutildenew[1:dimdiff,, drop = F]
  
  gx12 = y%x%(invx11%**%(mutrunc-mu[1:dimdiff, drop = F]))
  indic1 = vectranspose(dim1,dimdiff)
  gx12 = gss.submat(gx12,indic1,0)
  .omsymmetric <<- 1
  .omdiagonal<<- 0
  if(.condcovmeantrunc) {
    .xinvsymmetric <<- 1
    .xinvdiagonal <<- 0
    .xinvcorrelation <<- 0
    indic2 = matdup(gss.seqa(1,1,ddcov))%~%(gss.reshape((gss.seqa(ddcov+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    indic3 = gss.zeros(dim1,dimdiff)%~%matdup(gss.seqa(ddcov+dimdiff%**%dim1+1,1,dd1cov))
    indic = vecdup(indic2%|%indic3)
    gx11 = ginverse(x11)%**%gaomegab((y%**%x12),(mutrunc-mu[1:dimdiff, drop = F]))
    gx11 = gx11+gmutildenew[(dimdiff+1):(dimdiff+dimdiff%**%((dimdiff+1)%//%2)),, drop = F]
    ggx = gss.submat((gx11%|%gx12%|%gss.zeros(dd1cov,dim1)),indic,0)
    
  } else if(.condcovmeantrunc==0) {
    .xinvsymmetric <<- 1
    .xinvdiagonal <<- 0
    .xinvcorrelation <<- 1
    gx11 = ginverse(x11)%**%gaomegab((y%**%x12),(mutrunc-mu[1:dimdiff, drop = F]))

    if(dimdiff==2) {
      gx11 = gx11 + gmutildenew[dimdiff+2,, drop = F]
    }
    if(ddcor==0  ) {
      indic2 = 0%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    } else {
      indic2 = matndupdiagzero(gss.seqa(1,1,ddcor))%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    }
    if(dd1cor==0) {
      indic3 = gss.zeros(dim1,dimdiff)%~%0
    } else {
      indic3 = gss.zeros(dim1,dimdiff)%~%matndupdiagzero(gss.seqa(ddcor+dimdiff%**%dim1+1,1,dd1cor))
    }
    indic = vecndup(indic2%|%indic3)       
    if(ddcor==0 & dd1cor==0) {
      ggx=gx12
    } else if(ddcor==0 & dd1cor!=0) {
      ggx = gx12%|%(gss.zeros(dd1cor,dim1))
    } else if(ddcor!=0 & dd1cor==0) {
      ggx = gss.submat((gx11%|%gx12),indic,0)
    } else {
      ggx = gss.submat((gx11%|%gx12%|%(gss.zeros(dd1cor,dim1))),indic,0)
    } 
    
  }
  if(.cholesky & .condcovmeantrunc) {
    ggx = gcholeskycov(x)%**%ggx
  } else if(.cholesky & .condcovmeantrunc==0) {
    ggx = gcholeskycor(x)%**%ggx 
  }
  
  gc = gmutildenew[(dimdiff+dimdiff%**%((dimdiff+1)%//%2)+1):nrow(gmutildenew),, drop = F]
  output = list(ggy,gmu,ggx,gc)
  return(output)
}


multrunc = function(mu,cov,trpoint) { 
  
  mu = matrix(mu)
  
  m = nrow(mu)
  rr = chol(cov)
  lam = (trpoint-mu[1, drop = F])%//%rr[1,1, drop = F]
  mutilde = (-dnorm(lam))%//%pnorm(lam)
  munew = mutilde%|%gss.zeros(m-1,1)
  sigtilde2 = 1+(mutilde%**%(lam-mutilde))
  signew2 = (sigtilde2%|%gss.zeros(m-1,1))%~%(gss.zeros(1,m-1)%|%diag(m-1))  
  output = list((mu+t(rr)%**%munew),(t(rr)%**%signew2%**%rr))
  return(output)  
}  

gcondcov = function(y,x) { 
  
  if (length(y) == 1) y = matrix(y)
  
  dim1 = nrow(y)
  dim2 = nrow(x)
  dimdiff = dim2-dim1
  ddcov = (dimdiff%**%(dimdiff+1))%//%2   #/* these next four lines are to get a way of sequencing the gradients later to put in the right order of elements in matrix x */
  ddcor = (dimdiff%**%(dimdiff-1))%//%2
  dd1cov = (dim1%**%(dim1+1))%//%2
  dd1cor = (dim1%**%(dim1-1))%//%2
  x11 = x[1:dimdiff,1:dimdiff, drop = F]
  x12 = x[(dimdiff+1):nrow(x),1:dimdiff, drop = F]
  x22 = x[(dimdiff+1):nrow(x), (dimdiff+1):nrow(x), drop = F]
  invx11 = solve(x11)
  x22condcov = x22-x12%**%invx11%**%t(x12)
  
  .x2symmetric <<- 1
  .x2diagonal <<- 0
  .x1symmetric <<- 1
  .x1diagonal <<- 1 
  .x2correlation <<- 0
  
  output = gbothxomegax(y,x22condcov)
  ggy = output[[1]]
  gg2 = output[[2]]
  
  .x2symmetric <<- 1
  .x2diagonal <<- 0
  .x1symmetric <<- 0
  .x1diagonal <<- 0
  .x2correlation <<- 0
  
  output = gbothxomegax(x12,invx11)
  gg22 = output[[1]]
  gg23 = output[[2]]
  gg22 =-gg22%**%gg2
  indic1 = vectranspose(dim1,dimdiff)
  gg22new = gss.submat(gg22,indic1,0)
  
  if(.condcov) {
    
    .xinvsymmetric <<- 1
    .xinvdiagonal <<- 0
    .xinvcorrelation <<- 0
    
    indic2 = matdup(gss.seqa(1,1,ddcov))%~%(gss.reshape((gss.seqa(ddcov+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    indic3 = gss.zeros(dim1,dimdiff)%~%matdup(gss.seqa(ddcov+dimdiff%**%dim1+1,1,dd1cov))
    indic = vecdup(indic2%|%indic3)
    gg23 = -ginverse(x11)%**%gg23%**%gg2
    ggx = gss.submat((gg23%|%gg22new%|%gg2),indic,0)
  } else if(.condcov==0) {
    
    .xinvsymmetric <<- 1
    .xinvdiagonal <<- 0
    .xinvcorrelation <<- 1
    
    if(ddcor == 0) {
      indic2 = 0%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    } else {
      indic2 = matndupdiagzero(gss.seqa(1,1,ddcor))%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    }
    if(dd1cor==0) {
      indic3 = gss.zeros(dim1,dimdiff)%~%0
    } else {
      indic3 = gss.zeros(dim1,dimdiff)%~%matndupdiagzero(gss.seqa(ddcor+dimdiff%**%dim1+1,1,dd1cor))
    }
    
    indic = vecndup(indic2%|%indic3)       
    gg2new = gss.delif(gg2,vecdup(diag(dim1)))
    gg23 = -ginverse(x11)%**%gg23%**%gg2
    
    if(ddcor==0 & dd1cor==0) {
      ggx=gg22new
    } else if(ddcor==0 & dd1cor!=0) {
      ggx = gg22new%|%gg2new
    } else if(ddcor!=0 & dd1cor==0) {
      ggx = gss.submat((gg23%|%gg2new),indic,0)
    } else {
      ggx = gss.submat((gg23%|%gg22new%|%gg2new),indic,0)
    }    
  }
  
  if(.cholesky & .condcov) {
    ggx = gcholeskycov(x)%**%ggx
  } else if(.cholesky & .condcov==0) {
    ggx = gcholeskycor(x)%**%ggx 
  }
  output = list(ggy,ggx)
  return(output)
}

gcondcovtrunc = function(mu,x,gss.c) { 
  
  mu = matrix(mu)
  gss.c = matrix(gss.c)
  
  dim1 = nrow(x)-nrow(gss.c) 
  dim2 = nrow(x)
  dimdiff = dim2-dim1   
  ddcov = (dimdiff%**%(dimdiff+1))%//%2 
  ddcor = (dimdiff%**%(dimdiff-1))%//%2
  dd1cov = (dim1%**%(dim1+1))%//%2
  dd1cor = (dim1%**%(dim1-1))%//%2
  
  x11 = x[1:dimdiff,1:dimdiff, drop = F]
  x12 = x[(dimdiff+1):nrow(x),1:dimdiff, drop = F]
  x22 = x[(dimdiff+1):nrow(x), (dimdiff+1):nrow(x), drop = F]    
  
  if(dimdiff==1) {
    output = univariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    mutrunc = output[[1]]
    sigtrunc = output[[2]]
    output = gradunivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    gmutrunc = output[[1]]
    gsigtrunc = output[[2]]
  } else if(dimdiff==2) {
    output = bivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)        
    mutrunc = output[[1]]
    sigtrunc = output[[2]]
    output = gradbivariatenormaltrunc(mu[1:dimdiff, drop = F],x[1:dimdiff,1:dimdiff, drop = F],gss.c)
    gmutrunc = output[[1]]
    gsigtrunc = output[[2]]
  }
  invx11 = solve(x11)
  b = invx11%**%sigtrunc%**%invx11
  .cholesky=0
  if(.condcovsigtrunc==0) {
    .condcov<<-0
    .x2correlation<<-1
    .xinvsymmetric <<- 1
    .xinvdiagonal<<-0
    .xinvcorrelation<<-1
    
  } else {
    .condcov<<-1
    .x2correlation<<-0
    .xinvsymmetric <<- 1
    .xinvdiagonal<<-0
    .xinvcorrelation<<-0
    
  }
  output = gcondcov(diag(dim1),x)   
  gy = output[[1]]
  gd1 = output[[2]]
  
  
  .x2symmetric<<-1
  .x2diagonal<<-0
  .x1symmetric<<-1
  .x1diagonal<<-0
  .x2correlation <<- 0
  output = gbothxomegax(invx11,sigtrunc)
  ginvx11 = output[[1]]
  gbomega = output[[2]]
  
  gbpsi11= ginverse(x11)%**%ginvx11 
  .x2symmetric<<-1
  .x2diagonal<<-0
  .x1symmetric<<-0
  .x1diagonal<<-0
  .x2correlation <<- 0
  
  output = gbothxomegax(x12,b)
  gpsi12 = output[[1]]
  gb = output[[2]]
  gpsi11 = gbpsi11%**%gb
  gomega = gbomega%**%gb
  gsigtildenew = gsigtrunc%**%gomega
  
  gmu = gsigtildenew[1:dimdiff,, drop = F]%|%gss.zeros(dim1,ncol(gsigtildenew))
  gc = gsigtildenew[(nrow(gsigtildenew)-dimdiff+1):nrow(gsigtildenew),, drop = F]

  indic1 = vectranspose(dim1,dimdiff)
  gpsi12 = gss.submat(gpsi12,indic1,0)
  gpsi22 = gss.zeros(nrow(gd1),ncol(gd1))
  if(.condcov) {
    gpsi11 = gpsi11+gsigtildenew[(dimdiff+1):(nrow(gsigtildenew)-dimdiff),, drop = F]        
    indic2 = matdup(gss.seqa(1,1,ddcov))%~%(gss.reshape((gss.seqa(ddcov+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    indic3 = gss.zeros(dim1,dimdiff)%~%matdup(gss.seqa(ddcov+dimdiff%**%dim1+1,1,dd1cov))
    indic = vecdup(indic2%|%indic3)
    ggx = gd1+gss.submat((gpsi11%|%gpsi12%|%gpsi22),indic,0)
  } else if(.condcov==0) {
    if(dimdiff==2) {
      gpsi11 = gpsi11+gsigtildenew[dimdiff+2,, drop = F]
    }
    if(ddcor == 0) {
      indic2 = 0%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    } else {
      indic2 = matndupdiagzero(gss.seqa(1,1,ddcor))%~%(gss.reshape((gss.seqa(ddcor+1,1,dimdiff%**%dim1)),dimdiff,dim1))
    }
    if(dd1cor==0) {
      indic3 = gss.zeros(dim1,dimdiff)%~%0
    } else {
      indic3 = gss.zeros(dim1,dimdiff)%~%matndupdiagzero(gss.seqa(ddcor+dimdiff%**%dim1+1,1,dd1cor))
    }
    indic = vecndup(indic2%|%indic3) 
    if(ddcor==0 & dd1cor==0) {
      ggx=gd1+gpsi12
    } else if(ddcor==0 & dd1cor!=0) {
      ggx = gd1+(gpsi12%|%gss.zeros(dd1cor,dd1cov))
    } else if(ddcor!=0 & dd1cor==0) {
      ggx = gd1+gss.submat((gpsi11%|%gpsi12),indic,0)
    } else {
      ggx = gd1+gss.submat((gpsi11%|%gpsi12%|%gpsi22),indic,0)
    }           
  }
  
  
  if(.cholesky & .condcov) {
    ggx = gcholeskycov(x)%**%ggx
  } else if(.cholesky & .condcov==0) {
    ggx = gcholeskycor(x)%**%ggx 
  }
  output = list(gmu,ggx,gc)
  return(output)
}  

#/*****************************************************************************************************************************************************
#                 SET 3: THE CODES BELOW CONSTITUTE PROCEDURES FOR VARIOUS MULTIVARIATE CUMULATIVE NORMAL DISTRIBUTIONS AND GRADIENTS
#*****************************************************************************************************************************************************/

noncdfn = function(mu,sig,x) { 
  w = (x-mu)%//%(sqrt(sig))
  return(pnorm(w))
}   


gradnoncdfn = function(mu,sig2,x) { 
  invsig2 = 1%//%sig2
  invsig = sqrt(invsig2)
  w = (x-mu)%.*%invsig
  output = list(-invsig%.*%dnorm(w),dnorm(w)%.*%(-0.5%**%((invsig2)%.*%w)),dnorm(w)%.*%invsig)
  return(output)
}


noncdfbvn = function(mu,cov,x) { 
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(2,1))     
  rho = kk2[1,2, drop = F]
  return(gss.cdfbvn(kk1[1, drop = F],kk1[2, drop = F],rho))
}


gradcdfbvn = function(w1,w2,rho) { 
  rhotilde = sqrt((1-rho^2))
  tr1 = (w2-rho%.*%w1)%//%rhotilde 
  tr2 = (w1-rho%.*%w2)%//%rhotilde
  pdf2 = (1%//%rhotilde)%.*%(dnorm(w1)%.*%dnorm(tr1))
  gw1 = dnorm(w1)%.*%pnorm(tr1)
  gw2 = dnorm(w2)%.*%pnorm(tr2)
  grho = pdf2
  output = list(gw1,gw2,grho)
  return(output)
}    

gradnoncdfbvn = function(mu,cov,x) { 
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(2,1))     
  rho = kk2[1,2, drop = F]
  output = gradcdfbvn(kk1[1, drop = F],kk1[2, drop = F],rho) 
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  output = gradcorcov(kk1,sqrtom,kk2)
  gbcorcov = output[[1]]
  gomegacorcov = output[[2]]
  
  gmu = -(gw1%|%gw2)%//%sqrtom
  gcov = gbcorcov%**%(gw1%|%gw2)+gomegacorcov%**%grho
  gx = -gmu
  if(.cholesky==1) {
    gcov = gcholeskycov(cov)%**%gcov
  }  
  output = list(gmu,gcov,gx)
  return(output)
}    


gradcdfbvnbycdfn = function(w1,w2,rho) { 
  bivarcdf = gss.cdfbvn(w1,w2,rho)
  univarcdf = pnorm(w1)
  output  = gradcdfbvn(w1,w2,rho)
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  gw1new = (univarcdf%.*%gw1-bivarcdf%.*%dnorm(w1))%//%(univarcdf^2)
  output = list(gw1new,gw2%//%univarcdf,grho%//%univarcdf)
  return(output)
}    

gradnoncdfbvnbycdfn = function(mu,cov,x) { 
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(2,1))     
  rho = kk2[1,2, drop = F]
  output = gradcdfbvnbycdfn(kk1[1, drop = F],kk1[2, drop = F],rho) 
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  output = gradcorcov(kk1,sqrtom,kk2)
  gbcorcov = output[[1]]
  gomegacorcov = output[[2]]
  gmu = -(gw1%|%gw2)%//%sqrtom
  gcov = gbcorcov%**%(gw1%|%gw2)+gomegacorcov%**%grho
  gx = -gmu
  if(.cholesky==1) {
    gcov = gcholeskycov(cov)%**%gcov
  }
  output = list(gmu,gcov,gx)
  
  return(output)
}    


#/* cumulative distribution function of non-standard trivariate cumulative distribution function */
#/* mu - 3x1 vector, cov - 3x3 matrix, x - 3x1 vector of abscissa */

noncdftvn = function(mu,cov,x) { 
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(3,1))     
  return(gss.cdftvn(kk1[1, drop = F],kk1[2, drop = F],kk1[3, drop = F],kk2[1,2, drop = F],kk2[2,3, drop = F],kk2[1,3, drop = F]))
}  

gradcdftvnbycdfbvn = function(w,cor) {
  w1 = w[1, drop = F]
  w2 = w[2, drop = F]
  w3 = w[3, drop = F]
  rho12 = cor[1, drop = F]
  rho13 = cor[2, drop = F]
  rho23 = cor[3, drop = F]
  trivarcdf = gss.cdftvn(w1,w2,w3,rho12,rho23,rho13)
  bivarcdf = gss.cdfbvn(w1,w2,rho12)
  output  = gradcdfbvn(w1,w2,rho12)
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  output = gradcdftvn(w,cor)
  gwtri = output[[1]]
  grhotri = output[[2]]
  gw1w2new = (bivarcdf%**%gwtri[1:2, drop = F]-trivarcdf%**%(gw1%|%gw2))%//%(bivarcdf^2)
  grho12new = (bivarcdf%**%grhotri[1, drop = F]-trivarcdf%**%(grho))%//%(bivarcdf^2)
  
  output = list((gw1w2new%|%(gwtri[3, drop = F]%//%bivarcdf)),(grho12new%|%((grhotri[2:3,, drop = F])%//%bivarcdf)))
  return(output)
}

gradnoncdftvnbycdfbvn = function(mu,cov,x) {
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(3,1))     
  output = gradcdftvnbycdfbvn(kk1,kk2[1,2, drop = F]%|%kk2[1,3, drop = F]%|%kk2[2,3, drop = F]) 
  gw = output[[1]]
  grho = output[[2]]
  output = gradcorcov(kk1,sqrtom,kk2)
  gbcorcov = output[[1]]
  gomegacorcov = output[[2]]
  gmu = -(gw)%//%sqrtom
  
  gcov = gbcorcov%**%(gw)+gomegacorcov%**%grho
  gx = -gmu
  if(.cholesky==1) {
    gcov = gcholeskycov(cov)%**%gcov
  }
  output = list(gmu,gcov,gx)
  return(output)
} 

gradcdftvn = function(h,r) { 
  
  h = matrix(h)
  r = matrix(r)
  
  #//h = parm[1:3, drop = F]
  
  epst = 1e-7    
  pt = pi%//%2    
  h1 = h[1, drop = F]      
  h2 = h[2, drop = F]      
  h3 = h[3, drop = F]
  r12 = r[1, drop = F]     
  r13 = r[2, drop = F]     
  r23 = r[3, drop = F]
  
  gg = gss.zeros(1,6)
  enteredloop1 = 0 
  enteredloop2 = 0
  if(( abs(r12) - abs(r13) > 0.0000001 ) ) {
    h2 = h3    
    h3 = h[2, drop = F]
    r12 = r13  
    r13 = r[1, drop = F]
    enteredloop1=1
  }
  if(( abs(r13) - abs(r23) > 0.0000001 ) ) {
    h1 = h2    
    h2 = h[1, drop = F]
    r23 = r13  
    r13 = r[3, drop = F]
    enteredloop2=1
  }
  
  
  if((abs(h1) + abs(h2) + abs(h3) < epst) ) {
    gg[4] = ( 2%**% 1%//%sqrt(1-r12^2) %//%pi )%//%8        
    gg[5] = ( 2%**% 1%//%sqrt(1-r13^2) %//%pi )%//%8        
    gg[6] = ( 2%**% 1%//%sqrt(1-r23^2) %//%pi )%//%8
    d.tvn = gg
    
  } else if((abs(r12) + abs(r13) < epst)) {
    tempval1 = pnorm(h1)
    tempval2 =gss.cdfbvn( h2, h3, r23 ) 
    gg[1] = dnorm(h1)%**%tempval2
    temp.bvn.grad = gradcdfbvn( h2, h3, r23 )
    gg[2] = tempval1%**%temp.bvn.grad[1, drop = F]
    gg[3] = tempval1%**%temp.bvn.grad[2, drop = F]
    gg[6] = tempval1%**%temp.bvn.grad[3, drop = F]
    d.tvn = gg
    
  } else if((abs(r13) + abs(r23) < epst)) {
    tempval1 = pnorm(h3)
    tempval2 =gss.cdfbvn( h1, h2, r12 )
    gg[3] = dnorm(h3)%**%tempval2
    temp.bvn.grad = gradcdfbvn( h1, h2, r12 ) 
    gg[1] = tempval1%**%temp.bvn.grad[1, drop = F]
    gg[2] = tempval1%**%temp.bvn.grad[2, drop = F]
    gg[4] = tempval1%**%temp.bvn.grad[3, drop = F]
    d.tvn = gg
    
  } else if((abs(r12) + abs(r23) < epst)) {
    tvn = pnorm(h2)%**%gss.cdfbvn( h1, h3, r13 )
    tempval1 = pnorm(h2)
    tempval2 =gss.cdfbvn( h1, h3, r13 )
    gg[2] = dnorm(h2)%**%tempval2
    temp.bvn.grad = gradcdfbvn( h1, h3, r13 )
    gg[1] = tempval1%**%temp.bvn.grad[1, drop = F]
    gg[3] = tempval1%**%temp.bvn.grad[2, drop = F]
    gg[5] = tempval1%**%temp.bvn.grad[3, drop = F]
    d.tvn = gg
    
  } else if((1 - r23 < epst)) {
    tvn = gss.cdfbvn( h1, gss.minc( h2%|%h3 ), r12 )
    temp.bvn.grad = gradcdfbvn( h1, gss.minc( h2%|%h3 ), r12 )
    gg[1] = temp.bvn.grad[1, drop = F]
    if (h2 == h3) {
      gg[2] = temp.bvn.grad[2, drop = F]
      gg[3] = temp.bvn.grad[2, drop = F]
    } else if(h2 < h3) {
      gg[2] = temp.bvn.grad[2, drop = F]
    } else {
      gg[3] = temp.bvn.grad[2, drop = F]
    }
    gg[4] = temp.bvn.grad[3, drop = F]
    d.tvn = gg
    
  } else if((r23 + 1 < epst & h2 > -h3 )) {
    temp.bvn.grad  = gradcdfbvn( h1, h2, r12 ) 
    temp.bvn.grad2 = gradcdfbvn( h1, -h3, r12 )
    gg[1] = temp.bvn.grad[1, drop = F] - temp.bvn.grad2[1, drop = F]
    gg[2] = temp.bvn.grad[2, drop = F] 
    gg[3] = temp.bvn.grad2[2, drop = F]
    gg[4] = temp.bvn.grad[3, drop = F] - temp.bvn.grad2[3, drop = F]        
    d.tvn = gg
    
  } else {
    
    #//compute singular tvn value
    cdfbvn.h2.h3 = gss.cdfbvn( h2, h3, r23 )
    cdfn.h1      = pnorm(h1)
    tvn = cdfn.h1%**%cdfbvn.h2.h3
    
    d.a     = dnorm(h2) %.*% pnorm((h3- r23 %.*% h2) %//% sqrt(1-r23 %.*% r23) )
    d.b     = dnorm(h3) %.*% pnorm((h2- r23 %.*% h3) %//% sqrt(1-r23 %.*% r23) )
    d.corr  = (exp(-0.5%**%(h2^2 + h3^2 - 2%**%r23 %.*% h2 %.*% h3   ) %//% (1-r23 %.*% r23) )) %//% sqrt(1-r23 %.*% r23)%//%(2%**%pi)
    
    d.tvn = (dnorm(h1)%**%cdfbvn.h2.h3)%~%(cdfn.h1%**%d.a)%~%(cdfn.h1%**%d.b)%~%0%~%0%~%(cdfn.h1%**%d.corr)
    
    wg = matrix(c(0.1279381953467518,0.1258374563468280,0.1216704729278031,0.1155056680537265,0.1074442701159659,.09761865210411358,.08619016153195296,.07334648141108081,.05929858491543594,.04427743881742087,.02853138862893389,.01234122979998693))
    xg = matrix(c(.06405689286260559,0.1911188674736164,0.3150426796961635,0.4337935076260450,0.5454214713888396,0.6480936519369754,0.7401241915785546,0.8200019859739028,0.8864155270044012,0.9382745520027329,0.9747285559713096,0.9951872199970214))
    rua = asin( r12 )
    rub = asin( r13 )
    res = 0
    d.rua.r12 = 1%//%sqrt(1-r12^2)
    d.rub.r13 = 1%//%sqrt(1-r13^2)
    d.res.h1 = 0
    d.res.h2 = 0
    d.res.h3 = 0
    d.res.r12=0
    d.res.r13=0
    d.res.r23=0
    for(j in seq.int(1,12,1)) {
      fc = 0
      d.fc.h1 = 0
      d.fc.h2 = 0
      d.fc.h3 = 0
      d.fc.r12=0
      d.fc.r13=0
      d.fc.r23=0
      temp.sincs.grad = sincs.grad( rua%**%( 1 - xg[j, drop = F] )%//%2 ) 
      r12t = temp.sincs.grad[3, drop = F]
      rr2  = temp.sincs.grad[4, drop = F]
      d.r12t.r12 = temp.sincs.grad[1, drop = F]%**%( 1 - xg[j, drop = F] )%//%2%**%d.rua.r12
      d.rr2.r12  = temp.sincs.grad[2, drop = F]%**%( 1 - xg[j, drop = F] )%//%2%**%d.rua.r12
      
      temp.sincs.grad = sincs.grad( rub%**%( 1 - xg[j, drop = F] )%//%2 ) 
      r13t = temp.sincs.grad[3, drop = F]
      rr3  = temp.sincs.grad[4, drop = F]
      d.r13t.r13 = temp.sincs.grad[1, drop = F]%**%( 1 - xg[j, drop = F] )%//%2%**%d.rub.r13
      d.rr3.r13  = temp.sincs.grad[2, drop = F]%**%( 1 - xg[j, drop = F] )%//%2%**%d.rub.r13
      
      if(( abs(rua) > 0 )) {
        temp.pntgnd.grad = pntgnd.grad( h1, h2, h3, r13t, r23, r12t, rr2 )
        fc = fc + rua%**%temp.pntgnd.grad[8, drop = F]
        d.fc.h1 = d.fc.h1 + rua%**%temp.pntgnd.grad[1, drop = F]
        d.fc.h2 = d.fc.h2 + rua%**%temp.pntgnd.grad[2, drop = F]
        d.fc.h3 = d.fc.h3 + rua%**%temp.pntgnd.grad[3, drop = F]
        d.fc.r12 = d.fc.r12 + d.rua.r12 %**% temp.pntgnd.grad[8, drop = F] + rua %**%(temp.pntgnd.grad[6, drop = F]%**%d.r12t.r12+temp.pntgnd.grad[7, drop = F]%**%d.rr2.r12)
        d.fc.r13 = d.fc.r13 + rua %**%(temp.pntgnd.grad[4, drop = F]%**%d.r13t.r13)
        d.fc.r23 = d.fc.r23 + rua%**%temp.pntgnd.grad[5, drop = F]
      }
      if(( abs(rub) > 0 )) {
        temp.pntgnd.grad = pntgnd.grad( h1, h3, h2, r12t, r23, r13t, rr3 )
        fc = fc + rub%**%temp.pntgnd.grad[8, drop = F]
        d.fc.h1 = d.fc.h1 + rub%**%temp.pntgnd.grad[1, drop = F]
        d.fc.h2 = d.fc.h2 + rub%**%temp.pntgnd.grad[3, drop = F]
        d.fc.h3 = d.fc.h3 + rub%**%temp.pntgnd.grad[2, drop = F]
        d.fc.r12 = d.fc.r12 + rub %**%(temp.pntgnd.grad[4, drop = F]%**%d.r12t.r12)
        d.fc.r13 = d.fc.r13 + d.rub.r13 %**% temp.pntgnd.grad[8, drop = F] + rub %**%(temp.pntgnd.grad[6, drop = F]%**%d.r13t.r13+temp.pntgnd.grad[7, drop = F]%**%d.rr3.r13)
        d.fc.r23 = d.fc.r23 + rub%**%temp.pntgnd.grad[5, drop = F]
      }
      
      temp.sincs.grad = sincs.grad( rua%**%( 1 + xg[j, drop = F] )%//%2 ) 
      r12t = temp.sincs.grad[3, drop = F]
      rr2  = temp.sincs.grad[4, drop = F]
      
      d.r12t.r12 = temp.sincs.grad[1, drop = F]%**%( 1 + xg[j, drop = F] )%//%2%**%d.rua.r12
      d.rr2.r12  = temp.sincs.grad[2, drop = F]%**%( 1 + xg[j, drop = F] )%//%2%**%d.rua.r12
      
      temp.sincs.grad = sincs.grad( rub%**%( 1 + xg[j, drop = F] )%//%2 )
      r13t = temp.sincs.grad[3, drop = F]
      rr3  = temp.sincs.grad[4, drop = F]
      d.r13t.r13 = temp.sincs.grad[1, drop = F]%**%( 1 + xg[j, drop = F] )%//%2%**%d.rub.r13
      d.rr3.r13  = temp.sincs.grad[2, drop = F]%**%( 1 + xg[j, drop = F] )%//%2%**%d.rub.r13
      
      if(( abs(rua) > 0 ) ) {
        temp.pntgnd.grad = pntgnd.grad(h1, h2, h3, r13t, r23, r12t, rr2 )
        fc = fc + rua%**%temp.pntgnd.grad[8, drop = F]
        d.fc.h1 = d.fc.h1 + rua%**%temp.pntgnd.grad[1, drop = F]
        d.fc.h2 = d.fc.h2 + rua%**%temp.pntgnd.grad[2, drop = F]
        d.fc.h3 = d.fc.h3 + rua%**%temp.pntgnd.grad[3, drop = F]
        d.fc.r12 = d.fc.r12 + d.rua.r12 %**% temp.pntgnd.grad[8, drop = F] + rua %**%(temp.pntgnd.grad[6, drop = F]%**%d.r12t.r12+temp.pntgnd.grad[7, drop = F]%**%d.rr2.r12)
        d.fc.r13 = d.fc.r13 + rua %**%(temp.pntgnd.grad[4, drop = F]%**%d.r13t.r13)
        d.fc.r23 = d.fc.r23 + rua%**%temp.pntgnd.grad[5, drop = F]
      }
      if(( abs(rub) > 0 )) {
        temp.pntgnd.grad = pntgnd.grad( h1, h3, h2, r12t, r23, r13t, rr3 )
        fc = fc + rub%**%temp.pntgnd.grad[8, drop = F]
        d.fc.h1 = d.fc.h1 + rub%**%temp.pntgnd.grad[1, drop = F]
        d.fc.h2 = d.fc.h2 + rub%**%temp.pntgnd.grad[3, drop = F]
        d.fc.h3 = d.fc.h3 + rub%**%temp.pntgnd.grad[2, drop = F]
        d.fc.r12 = d.fc.r12 + rub %**%(temp.pntgnd.grad[4, drop = F]%**%d.r12t.r12)
        d.fc.r13 = d.fc.r13 + d.rub.r13 %**% temp.pntgnd.grad[8, drop = F] + rub %**%(temp.pntgnd.grad[6, drop = F]%**%d.r13t.r13+temp.pntgnd.grad[7, drop = F]%**%d.rr3.r13)
        d.fc.r23 = d.fc.r23 + rub%**%temp.pntgnd.grad[5, drop = F]
      }
      res = res + wg[j, drop = F]%**%fc 
      d.res.h1  = d.res.h1  + wg[j, drop = F]%**%d.fc.h1
      d.res.h2  = d.res.h2  + wg[j, drop = F]%**%d.fc.h2
      d.res.h3  = d.res.h3  + wg[j, drop = F]%**%d.fc.h3
      d.res.r12 = d.res.r12 + wg[j, drop = F]%**%d.fc.r12
      d.res.r13 = d.res.r13 + wg[j, drop = F]%**%d.fc.r13
      d.res.r23 = d.res.r23 + wg[j, drop = F]%**%d.fc.r23
    }
    tvn = tvn + res%//%( 4%**%pi )   
    d.tvn = d.tvn + (d.res.h1%~%d.res.h2%~%d.res.h3%~%d.res.r12%~%d.res.r13%~%d.res.r23)%//%( 4%**%pi )  
  }
  
  if(( enteredloop2 == 1 ) ) {
    d.tvn.unsorted = d.tvn
    d.tvn.unsorted[1] = d.tvn[2, drop = F]
    d.tvn.unsorted[2] = d.tvn[1, drop = F]
    d.tvn.unsorted[5] = d.tvn[6, drop = F]
    d.tvn.unsorted[6] = d.tvn[5, drop = F]
    d.tvn = d.tvn.unsorted
  }   
  if(( enteredloop1 == 1 ) ) {
    d.tvn.unsorted = d.tvn
    d.tvn.unsorted[2] = d.tvn[3, drop = F]
    d.tvn.unsorted[3] = d.tvn[2, drop = F]
    d.tvn.unsorted[4] = d.tvn[5, drop = F]
    d.tvn.unsorted[5] = d.tvn[4, drop = F]
    d.tvn = d.tvn.unsorted
  }
  
  output = list((t((d.tvn[,1:3, drop = F]))),(t((d.tvn[,4:6, drop = F]))))
  return(output)
}

#/****************************************************************************************************************************************************
#                            proc sinc_grad(x) and pntgnd_grad are auxiliary procedures for the gradient of CDFTVN
#****************************************************************************************************************************************************/

sincs.grad = function(x) {
  ee = ( pi%//%2 - abs(x) )^2
  
  if(x<0) {
    d.ee = 2%**%( pi%//%2 + x )
  } else {
    d.ee = 2%**%( x - pi%//%2  )
  }
  
  if((ee < 5e-5)) {
    cs = ee%**%( 1 - ee%**%( 1 - 2%**%ee%//%15 )%//%3 )
    d.cs = ( 1 - ee%**%( 1 - 2%**%ee%//%15 )%//%3 )     +       ee%**%(-  (    ( 1 - 2%**%ee%//%15 )%//%3   -   ee %**% ( 2%//%45)    )         )
    d.cs = d.ee%**%d.cs
    
    signx = x%//%abs(x)
    sx = ( 1 - ee%**%( 1 - ee%//%12 )%//%2 )%**%signx
    if(x > 0) {
      d.sx = - (  ( 1 - ee%//%12 )%//%2 +ee%**%( - 1%//%12 )%//%2 )
    } else {
      d.sx =  (  ( 1 - ee%//%12 )%//%2 +ee%**%( - 1%//%12 )%//%2 )
    }
    d.sx = d.sx%**%d.ee        
  } else {
    sx = sin(x)
    d.sx = cos(x)
    cs = 1 - sx%**%sx
    d.cs = -2%**%sx%**%d.sx
  }
  
  return(d.sx%~%d.cs%~%sx%~%cs)
}

pntgnd.grad = function(ba, bb, bc, ra, rb, r, rr) {
  
  f = 0
  d.ba = 0
  d.bb=0
  d.bc=0
  d.ra=0
  d.rb=0
  d.r = 0
  d.rr=0
  dt = rr%**%( rr - ( ra - rb )^2 - 2%**%ra%**%rb%**%( 1 - r ) )
  
  d.dt.ra = rr%**%(  - 2%**%( ra - rb ) - 2%**%rb%**%( 1 - r ) )
  d.dt.rb = rr%**%(   2%**%( ra - rb ) - 2%**%ra%**%( 1 - r ) )
  d.dt.r  = rr%**%(  2%**%ra%**%rb )
  d.dt.rr = ( rr - ( ra - rb )^2 - 2%**%ra%**%rb%**%( 1 - r ) ) + rr%**%( 1  )
  
  if((dt > 0  )) {
    bt = ( bc%**%rr + ba%**%( r%**%rb - ra ) + bb%**%( r%**%ra -rb ) )%//%sqrt(dt) 
    ft = ( ba - r%**%bb )^2%//%rr + bb%**%bb
    
    d.bt.dt = -0.5%**%bt%//%dt
    
    d.bt.ba = (  ( r%**%rb - ra ) )%//%sqrt(dt) 
    d.ft.ba = 2%**%( ba - r%**%bb )%//%rr
    d.bt.bb = (  ( r%**%ra - rb ) )%//%sqrt(dt) 
    d.ft.bb = 2%**%( ba - r%**%bb )%**%(-r)%//%rr + 2%**%bb
    d.bt.bc = ( rr  )%//%sqrt(dt)       
    
    d.bt.ra = ( ba%**%( -1 ) + bb%**%( r ) )%//%sqrt(dt) +  d.bt.dt%**%d.dt.ra 
    d.bt.rb = ( ba%**%( r ) + bb%**%( -1 ) )%//%sqrt(dt) +  d.bt.dt%**%d.dt.rb   
    d.bt.r  = ( ba%**%( rb) + bb%**%( ra) )%//%sqrt(dt)  +  d.bt.dt%**%d.dt.r 
    d.ft.r  = 2%**%( ba - r%**%bb )%**%(-bb)%//%rr 
    d.bt.rr = ( bc  )%//%sqrt(dt)                  +  d.bt.dt%**%d.dt.rr 
    d.ft.rr = -( ba - r%**%bb )^2%//%(rr^2) 
    
    if((bt > -10 & ft < 100)) {
      f = exp( -ft%//%2 )
      d.ba =  f%**%(-1%//%2)%**%d.ft.ba
      d.bb =  f%**%(-1%//%2)%**%d.ft.bb
      d.r  =  f%**%(-1%//%2)%**%d.ft.r
      d.rr =  f%**%(-1%//%2)%**%d.ft.rr
      if(bt < 10) { 
        d.ba =  d.ba%**%pnorm(bt) + f%**%dnorm(bt)%**%d.bt.ba
        d.bb =  d.bb%**%pnorm(bt) + f%**%dnorm(bt)%**%d.bt.bb
        d.bc =                  f%**%dnorm(bt)%**%d.bt.bc
        d.ra =                  f%**%dnorm(bt)%**%d.bt.ra
        d.rb =                  f%**%dnorm(bt)%**%d.bt.rb
        d.r  =   d.r%**%pnorm(bt) + f%**%dnorm(bt)%**%d.bt.r
        d.rr =  d.rr%**%pnorm(bt) + f%**%dnorm(bt)%**%d.bt.rr
        f = f%**%pnorm(bt) 
        
      }
    }
  }
  
  gg = d.ba%~%d.bb%~%d.bc%~%d.ra%~%d.rb%~%d.r%~%d.rr%~%f
  return(gg)
}

cdfqvn = function(a,rr) { 
  
  a = matrix(a)
  
  p = gss.zeros(2,1)
  #/* this next line is needed, because truncated values go bizzarre if(truncation happens on abscissa less than -6 */) {
  a = -6.0%**%(a < -6.0)+a%.*%(a>=-6.0)     
  temp1 = gss.seqa(1,1,nrow(a))
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(nrow(a),1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F])
  output = multruncbivariate(mutemp,sigtemp,atemp[1:2, drop = F]) 
  mutemp = output[[1]]
  sigtemp = output[[2]]
  p[2] = (noncdfbvn(mutemp[3:4, drop = F],sigtemp[3:4,3:4, drop = F],atemp[3:4, drop = F]))%//%(noncdfn(mutemp[3, drop = F],sigtemp[3,3, drop = F],atemp[3, drop = F]))    
  return(p[1, drop = F]%**%p[2, drop = F])
}

#/* cumulative distribution function of non-standard approximated quadrivariate cumulative distribution function */
#/* mu - 4x1 vector, cov - 4x4 matrix, x - 4x1 vector of abscissa */

noncdfqvn = function(mu,cov,x) { 
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(4,1))     
  output = list(cdfqvn(kk1,kk2))
  return(output)
}   

gradcdfqvn = function(a,rr) { 
  
  a = matrix(a)
  
  p = gss.zeros(2,1)
  #/* temp1 = vecindascending(a)  */
  temp1 = gss.seqa(1,1,nrow(a))
  atemp = a[temp1, drop = F]
  mutemp=gss.zeros(nrow(a),1) 
  sigtemp = rr[temp1,temp1, drop = F]
  p[1] =gss.cdftvn(atemp[1, drop = F],atemp[2, drop = F],atemp[3, drop = F],sigtemp[1,2, drop = F],sigtemp[2,3, drop = F],sigtemp[1,3, drop = F])
  output = multruncbivariate(mutemp,sigtemp,atemp[1:2, drop = F])
  mutemp1 = output[[1]]
  sigtemp1 = output[[2]]
  p[2] = (noncdfbvn(mutemp1[3:4, drop = F],sigtemp1[3:4,3:4, drop = F],atemp[3:4, drop = F]))%//%(noncdfn(mutemp1[3, drop = F],sigtemp1[3,3, drop = F],atemp[3, drop = F]))
  .condcovsigtrunc<<-0
  .condcovmeantrunc<<-0
  output = gcondmeantrunc(gss.eye(2),mutemp,sigtemp,atemp[1:2, drop = F])
  gy = output[[1]]
  gmumean = output[[2]]
  gxmean = output[[3]]
  gcmean = output[[4]]
  output = gcondcovtrunc(mutemp,sigtemp,atemp[1:2, drop = F])    
  gmucov = output[[1]]
  gxcov = output[[2]]
  gccov = output[[3]]
  
  gcumulmusig = gxmean%~%gxcov
  gcumulc = (gcmean%~%gccov)%|%gss.zeros(2,5)
  gc = gss.zeros(4,1)
  output = gradnoncdfbvnbycdfn(mutemp1[3:4, drop = F],sigtemp1[3:4,3:4, drop = F],atemp[3:4, drop = F])
  gfrompmu = output[[1]]
  gfrompcov = output[[2]]
  gfrompc = output[[3]]
  grhorho = (p[1, drop = F])%**%(gcumulmusig[,1:2, drop = F]%**%gfrompmu+(gcumulmusig[,3:5, drop = F])%**%gfrompcov)
  gc = (p[1, drop = F])%**%(gcumulc[,1:2, drop = F]%**%gfrompmu+gcumulc[,3:5, drop = F]%**%gfrompcov)
  
  gc[3:4] = (p[1, drop = F])%**%(gfrompc)  #/* getting gradient contribution directly from noncdfbvn for absiccae except the first two */       
  
  output = gradcdftvn(atemp[1:3, drop = F],(sigtemp[1,2, drop = F]%|%sigtemp[1,3, drop = F]%|%sigtemp[2,3, drop = F]))
  gc1 = output[[1]]
  grho1 = output[[2]]
  
  grhorho = grhorho+p[2, drop = F]%**%(grho1[1:2,, drop = F]%|%0%|%grho1[3,, drop = F]%|%0%|%0)
  gc = gc+p[2, drop = F]%**%(gc1%|%0)
  
  tempnew = gss.seqa(1,1,4)%~%temp1
  tempnew = gss.sortc(tempnew,2)
  tempnew = tempnew[,1, drop = F]
  gc = (gc[tempnew,, drop = F])
  grhorho = matndupdiagzerofull((grhorho)) 
  grhorho = gss.submat(grhorho,tempnew,tempnew)
  grhorho = (vecndup(grhorho))
  output = list(gc,grhorho)
  return(output)
  
}


gradcdfqvnbycdfbvn = function(w,cor) { 
  
  w = matrix(w)
  
  w1 = w[1,, drop = F]
  w2 = w[2,, drop = F]
  rho12 = cor[1,2, drop = F]
  quadvarcdf = cdfqvn(w,cor)
  bivarcdf = gss.cdfbvn(w1,w2,rho12)
  output  = gradcdfbvn(w1,w2,rho12)
  gw1 = output[[1]]
  gw2 = output[[2]]
  grho = output[[3]]
  output = gradcdfqvn(w,cor)
  gwquad = output[[1]]
  grhoquad = output[[2]]
  
  gw1w2new = (bivarcdf%**%gwquad[1:2, drop = F]-quadvarcdf%**%(gw1%|%gw2))%//%(bivarcdf^2)
  grho12new = (bivarcdf%**%grhoquad[1, drop = F]-quadvarcdf%**%(grho))%//%(bivarcdf^2)
  
  output = list((gw1w2new%|%(gwquad[3:4,, drop = F]%//%bivarcdf)),(grho12new%|%(grhoquad[2:6,, drop = F]%//%bivarcdf)))
  return(output)
}

gradnoncdfqvnbycdfbvn = function(mu,cov,x) { 
  
  
  mu = matrix(mu)
  x = matrix(x)
  
  om = gss.diag(cov)
  sqrtom = sqrt(om)
  kk1 = (x-mu)%//%sqrtom
  kk2 = cov2cor(cov)  
  #/* this next line is needed because some diagonal elements of kk2 computed above are slightly greater than 1 because of the algebra in computing
  #correlations from corrvc, so ensuring that diagonal elements are no more than 1 */
  kk2 = gss.diagrv(kk2,gss.ones(4,1))     
  output = gradcdfqvnbycdfbvn(kk1,kk2) 
  gw = output[[1]]
  grho = output[[2]]
  output = gradcorcov(kk1,sqrtom,kk2)
  gbcorcov = output[[1]]
  gomegacorcov = output[[2]]
  gmu = -(gw)%//%sqrtom
  gcov = gbcorcov%**%(gw)+gomegacorcov%**%grho
  gx = -gmu
  if(.cholesky==1) {
    gcov = gcholeskycov(cov)%**%gcov
  }
  output = list(gmu,gcov,gx)
  return(output)
} 


#/*****************************************************************************************************************************************************
#                SET 4: THE CODES BELOW ARE FOR LDLT BLOCK DECOMPOSITIONS AND LDLT RANK-1/RANK-2 UPDATING
#*****************************************************************************************************************************************************/

ldltblock = function(a,m) { 
  n=nrow(a)
  l = diag(n) 
  d = gss.zeros(n,n)
  for(i in seq.int(1,n-m,m)) {
    j=i+m-1
    h=j+1
    
    d[i:j,i:j] =a[i:j,i:j, drop = F]
    l[h:n,i:j] =a[h:n,i:j, drop = F]%**%solve(d[i:j,i:j, drop = F])
    a[h:n,h:n] =a[h:n,h:n, drop = F]-l[h:n,i:j, drop = F]%**%d[i:j,i:j, drop = F]%**%(t((l[h:n,i:j, drop = F])))
  }
  d[h:n,h:n] =a[h:n,h:n, drop = F]   
  output = list(l,d)
  return(output) 
}

ldltupspecial = function(l,d,omega,m) { 
  n=nrow(l)
  n1 = n-m
  m1 = omega
  d2 = d[(m+1):n, (m+1):n, drop = F]
  b = gss.zeros(n1,m)
  l22 = l[(m+1):n, (m+1):n, drop = F]
  z = (solve(l22))%**%l[(m+1):n,1:m, drop = F]
  newl = diag(n1)
  newd = gss.zeros(n1,n1)
  if(n1<=m) {
    newl = diag(n1)
    newd = z[1:n1,, drop = F]%**%omega%**%(t((z[1:n1,, drop = F])))+d2
  } else {
    for(i in seq.int(1,n1,m)) {
      
      if(i+m-1<=n1) {
        
        gss.t = m1%**%(t((z[i:(i+m-1),, drop = F])))
        newd[i:(i+m-1),i:(i+m-1)] = d2[i:(i+m-1),i:(i+m-1), drop = F] +z[i:(i+m-1),, drop = F]%**%gss.t
        dinv = solve(newd[i:(i+m-1),i:(i+m-1), drop = F])
        b[i:(i+m-1),] =(t((gss.t%**%dinv)))
        m1 = m1-gss.t%**%dinv%**%t(gss.t)
        
        if(i<n1-m+1) {
          if(i!=n1-m) {
            newl[(i+m):(i+2%**%m-1),1:(i+m-1)] =z[(i+m):(i+2%**%m-1),, drop = F]%**%(t((b[1:(i+m-1),, drop = F])))
          } else if(i==n1-m) {
            newl[i+m,1:(i+m-1)] =z[i+m,, drop = F]%**%(t((b[1:(i+m-1),, drop = F]))) 
          }
        }  
      } else if(i==n1) {
        gss.t = m1%**%(t((z[i,, drop = F])))
        newd[i,i] = d2[i,i, drop = F] +z[i,, drop = F]%**%gss.t
      }
    }
  }
  output = list(l22%**%newl,newd)
  return(output) 
}  


#/*****************************************************************************************************************************************************
#                         SET 5: THE CODES BELOW CONSTITUTE PROCEDURES FOR MISCELLANEOUS MATRIX MANIPULATIONS
#*****************************************************************************************************************************************************/

vecindascending = function(r) { 
  
  r = matrix(r)
  
  c2=r
  c3={}
  for(i in seq.int(1,nrow(r),1)) {
    c1 = gss.minindc(c2)
    c2[c1]=gss.maxc(c2)+1
    c3=c3%|%c1
  }
  return(c3)
}


vecndup = function(r) { 
  if(ncol(r) != nrow(r)) {
    print("error: rows & columns not equal")
  }
  h = {}
  gss.c=1
  while(!(gss.c >ncol(r)-1)) {
    h = h%~%r[gss.c, (gss.c+1):ncol(r), drop = F]
    gss.c=gss.c+1
  }
  return(t(h))
}


vectranspose = function(r,gss.c) { 
  c1 = gss.seqa(1,1,r%**%gss.c)
  c2 = gss.reshape(c1,r,gss.c)
  c3 = gss.vecr((t(c2)))
  return(c3)
}

vecsymmetry = function(r) { 
  temp1={}
  for(gss.c in seq.int(1,r,1)) {
    for(j in seq.int(gss.c,r,1)) {
      temp = gss.zeros(r,r)
      temp[gss.c,j]=1
      temp[j,gss.c]=1
      temp1 = temp1%|%gss.reshape(temp,1,r%**%r)     #/* temp1 is mechanism to add gradients of symmetric elements & get derivatives with respect to only upper triangular elements of capomega */         
    }
  }
  return(temp1)
}


gasymtosym = function(x) { 
  k = sqrt(nrow(x))
  n = sqrt(ncol(x))
  temp1 = vecsymmetry(k)    
  tempselmatrix = gss.upmat(gss.ones(n,n))
  tempselmatrix = t((gss.reshape(tempselmatrix,1,n^2)))  #/* tempselmatrix is because t(x)%x%t(x) gives derivatives including a12 & a21, a13 
  #and a31 etc., but should focus on only one since a12=a21 & a13=a31 tempsel deletes lower triangular elements of a */
  newx = t((gss.selif(t(x),tempselmatrix)))
  gsym = temp1%**%newx
  return(gsym)
}


vecdup = function(r) { 
  
  if(ncol(r) != nrow(r)) {
    print("error: rows & columns not equal")
  }
  h = {}
  gss.c=1
  while(!(gss.c >ncol(r))) {
    h = h%~%r[gss.c,gss.c:ncol(r), drop = F]
    gss.c=gss.c+1
  }
  return(t(h))
}


matdup = function(r) { 
  
  r = matrix(r, ncol = 1)
  
  p = (-1+(1+8%**%nrow(r))^0.5)%//%2
  if((p%**%(p+1)%//%2) != nrow(r)) {
    print("error in input vector")
    stop()
  } else { 
    l = gss.zeros(p, p)
    gss.c = 1
    sk=1
    while(!(gss.c>p)) {
      l[gss.c,gss.c:p] = t((r[sk:(sk+p-gss.c), drop = F]))
      sk=sk+p-gss.c+1
      gss.c=gss.c+1
    }
    return(l)
  }
}


matdupfull = function(r) { 
  
  r = matrix(r, ncol = 1)
  
  p = (-1+(1+8%**%nrow(r))^0.5)%//%2
  if((p%**%(p+1)%//%2) != nrow(r)) {
    print("error in input vector")
    stop()
  } else { 
    l = gss.zeros(p, p)
    gss.c = 1
    sk=1
    while(!(gss.c>p)) {
      l[gss.c,gss.c:p] = t((r[sk:(sk+p-gss.c), drop = F]))
      sk=sk+p-gss.c+1
      gss.c=gss.c+1
    }
    ldiag=gss.diag(l)
    l=l+t(l)
    l=gss.diagrv(l,ldiag)
    return(l)
  }
}


matndupdiagzero = function(r) { 
  
  r = matrix(r, ncol = 1)
  
  p = (1+(1+8%**%nrow(r))^0.5)%//%2
  if((p%**%(p-1)%//%2) != nrow(r)) {
    print("error in input vector")
    stop()
  } else {
    l = gss.zeros(p, p)
    gss.c=1
    sk=1
    while(!(gss.c==p)) {
      l[gss.c, (gss.c+1):p] = t((r[sk:(sk+p-gss.c-1), drop = F]))
      sk=sk+p-gss.c
      gss.c=gss.c+1
    }
    
  }  
  return(l)
}


matndupdiagzerofull = function(r) { 
  
  r = matrix(r, ncol = 1)
  
  p = (1+(1+8%**%nrow(r))^0.5)%//%2
  if((p%**%(p-1)%//%2) != nrow(r)) {
    print("error in input vector")
    stop()
  } else {
    l = gss.zeros(p, p)
    gss.c=1
    sk=1
    while(!(gss.c==p)) {
      l[gss.c, (gss.c+1):p] = t((r[sk:(sk+p-gss.c-1), drop = F]))
      sk=sk+p-gss.c
      gss.c=gss.c+1
    }
  }  
  return(l+t(l))
}


gcholeskycov = function(capomega) { 
  litomega = sqrt(gss.diag(capomega))
  parmdim = nrow(litomega)
  parmcov = parmdim%**%(parmdim+1)%//%2
  omegastar.chol = chol(capomega)
  omegastar.chol.diag = gss.diag(omegastar.chol)   
  dg = gss.zeros(parmcov,parmcov)
  xsig = vecdup(omegastar.chol)
  j=1
  mm = 1
  if(parmdim>1) {
    while(!(j > parmcov  )) {
      gss.c=j
      l=1
      while(!(gss.c> parmcov)) {
        dg[(j+l-1):(j+parmdim-mm),gss.c:(gss.c+parmdim-l+1-mm)] = diag(parmdim-l+2-mm)%**%xsig[j+l-1, drop = F]
        dg[j+l-1,gss.c:(gss.c+parmdim-l+1-mm)] = dg[j+l-1,gss.c:(gss.c+parmdim-l+1-mm), drop = F]+t((xsig[(j+l-1):(j+parmdim-mm), drop = F]))
        gss.c=gss.c+parmdim-l+1-mm+1
        l=l+1
      }
      j = j+parmdim-mm+1
      mm=mm+1
    }
  } else if(parmdim==1) {
    dg=1 
  } 
  return(dg)
}

.test.gcholeskycov = function() {
  V = matrix(c( 1, -0.2,  0.1,  0, 
                -0.2,  2,  0.2, -0.1, 
                0.1,  0.2,  2, -0.2, 
                0,   -0.1, -0.2,  1), nrow = 4)
  gcholeskycov(V)
}
 
gcholeskycor = function(capomega) { 
  litomega = sqrt(gss.diag(capomega))
  parmdim = nrow(litomega)
  parmcor = parmdim%**%(parmdim-1)%//%2
  omegastar.chol = chol(capomega)
  omegastar.chol.diag = gss.diag(omegastar.chol)   
  dg = gss.zeros(parmcor,parmcor)
  j=1
  mm = 1
  if(parmdim >2) {
    for(mm in seq.int(1,parmdim,1)) {
      gss.c=j
      l=1        
      if(j<parmcor) {
        ss = (gss.upmat(gss.ones(parmdim-mm-1,parmdim-mm-1)))%.*%omegastar.chol[mm, (mm+2):parmdim, drop = F]
        ss1=ss-(omegastar.chol[(mm+1):(parmdim-1), (mm+2):parmdim, drop = F]%.*%(t((omegastar.chol[mm, (mm+1):(parmdim-1), drop = F])))%.*%(1%//%(omegastar.chol.diag[(mm+1):(parmdim-1),1, drop = F]))%.*%gss.upmat(gss.ones(parmdim-mm-1,parmdim-mm-1)))
      }
      while(!(gss.c > parmcor)) {
        dg[(j+l-1):(j+parmdim-mm-1),gss.c:(gss.c+parmdim-l-mm)] = diag(parmdim-l+1-mm)%**%omegastar.chol[mm,mm+l-1, drop = F]
        if(gss.c<parmcor) {
          dg[j+l-1, (gss.c+parmdim-l-mm+1):(gss.c+2%**%parmdim-2%**%l-2%**%mm)] = ss1[l,l:ncol(ss), drop = F]
        }
        gss.c=gss.c+parmdim-l+1-mm
        l=l+1
      }
      j = j+parmdim-mm
    }
  } else if(parmdim==2) {
    dg = 1   
  }
  return(dg)
}

.test.gcholeskycor = function() {
    V = matrix(c( 1, -0.2,  0.1,  0, 
                  -0.2,  1,  0.2, -0.1, 
                  0.1,  0.2,  1, -0.2, 
                  0,   -0.1, -0.2,  1), nrow = 4)
    gcholeskycor(V)
}


gbothxomegax = function(x1,x2) { 
  n = nrow(x1)
  k = ncol(x1)
  temp = t(x2)%**%t(x1)    
  t1={}
  for(i in seq.int(1,n,1)) {
    t1 = t1%~%(diag(n)%x%temp[,i, drop = F])
  }
  gx1cov = (diag(n)%x%(x2%**%t(x1)))+t1
  gx2cov = t(x1)%x%t(x1)      
  if(.x2symmetric) {
    tempselmatrix = gss.upmat(gss.ones(n,n))
    tempselmatrix = t((gss.reshape(tempselmatrix,1,n^2)))
    #/* tempselmatrix is because t(x)%x%t(x) gives derivatives including a12 & a21, a13                                                           & a31 etc., but should focus on only one since a12=a21 & a13=a31
    #tempsel deletes lower triangular elements of a */
    gx1cov = t((gss.selif(t(gx1cov),tempselmatrix)))
    gx2cov = gasymtosym(gx2cov)
    tempselmatrix = vecdup(diag(k))
    if(.x1symmetric) {
      temp1 = vecsymmetry(n) 
      gx1cov=temp1%**%gx1cov
    }
    if(.x2diagonal==0) {
      if(.x1diagonal) {
        gx1cov = gss.selif(gx1cov,tempselmatrix)
      }    
    } else if(.x2diagonal) {
      gx2cov = gss.selif(gx2cov,tempselmatrix)
      if(.x1diagonal) {
        gx1cov = gss.selif(gx1cov,tempselmatrix)
        gx1cov = t((gss.selif(t(gx1cov),tempselmatrix)))               
        gx2cov = t((gss.selif(t(gx2cov),tempselmatrix))) 
      }
    }
    if(.x2correlation & .x2diagonal) {
      gx2cov=0
    } else if(.x2correlation & .x2diagonal==0) {
      gx2cov=gss.delif(gx2cov,tempselmatrix)
    }
  }
  output = list(gx1cov,gx2cov)
  return(output)
}

.test.gbothxomegax = function() {
  V = matrix(c( 0.5, -0.2,  0.1,  0, 
                -0.2,  0.3,  0.2, -0.1, 
                0.1,  0.2,  0.5, -0.2, 
                0,   -0.1, -0.2,  0.3), nrow = 4)
  x = matrix(1:8, nrow = 2)
  gbothxomegax(x, V)
}

ginverse = function(x) { 
  k = nrow(x)
  if(.xinvsymmetric==0) {
    g1 = -solve((t(x)))%x%solve(x)  
  } else if(.xinvsymmetric==1) {
    tempselmatrix = vecdup(diag(k))
    g1 = -solve((t(x)))%x%solve(x)
    g1 = gasymtosym(g1)
    if(.xinvdiagonal==1) {
      g1 = gss.selif(g1,tempselmatrix)
      g1 = (t((gss.selif(t(g1),tempselmatrix))))
    }       
    if(.xinvcorrelation & .xinvdiagonal==0) {
      g1=gss.delif(g1,tempselmatrix)
    }
    if(.xinvcorrelation & .xinvdiagonal) {
      g1=0           
    }
  }
  return(g1)
}   

.test.ginverse = function(x) {
  V = matrix(c( 0.5, -0.2,  0.1,  0, 
                -0.2,  0.3,  0.2, -0.1, 
                0.1,  0.2,  0.5, -0.2, 
                0,   -0.1, -0.2,  0.3), nrow = 4)
  ginverse(V)
}


gradcorcov = function(gss.c,litomega,omegastar) { 
  
  litomega = matrix(litomega, ncol = 1)
  gss.c = matrix(gss.c, ncol=1)
  
  dim.lcl = nrow(litomega)
  row.lcl = dim.lcl%**%(dim.lcl+1)/2
  col.lcl = dim.lcl%**%(dim.lcl-1)/2
  
  nu1.lcl=litomega^2
  nu2.lcl=litomega
  temp = gss.diagrv(gss.zeros(dim.lcl,dim.lcl),litomega)
  nuomega.lcl = temp%**%omegastar%**%temp
  
  grm.lcl = gss.zeros(row.lcl,dim.lcl)
  cc1.lcl=0.5%**%((gss.c)/(nu1.lcl))
  l.lcl = 0
  for(idimvn in seq.int(1,dim.lcl,1)) {
    grm.lcl[l.lcl+1,idimvn]=cc1.lcl[idimvn, drop = F]
    l.lcl = l.lcl + (dim.lcl+1-idimvn)
  }
  
  grc.lcl = gss.zeros(row.lcl,col.lcl)
  j=1
  mm = 1
  c1=1
  while(!(j== row.lcl  )) {
    ndiag1 = 1/(nu2.lcl[mm, drop = F]%**%nu2.lcl[(mm+1):dim.lcl, drop = F])
    ndiag1 = gss.diagrv(diag(dim.lcl-mm),ndiag1)
    nhor1 = (-0.5)%**%(omegastar[mm, (mm+1):dim.lcl, drop = F]/nuomega.lcl[mm,mm, drop = F])
    nhor2 = (-0.5)%**%(omegastar[mm, (mm+1):dim.lcl, drop = F]/(t((nu1.lcl[(mm+1):dim.lcl, drop = F]))))
    grc.lcl[(j+1):(j+dim.lcl-mm),c1:(c1+dim.lcl-mm-1)]=ndiag1
    grc.lcl[j,c1:(c1+dim.lcl-1-mm)]=nhor1
    k=j
    l=1
    while(!(l==dim.lcl+1-mm)) {
      grc.lcl[k+dim.lcl+2-l-mm,l+c1-1]=nhor2[l, drop = F]
      k=k+dim.lcl+2-l-mm
      l=l+1
    }
    c1=c1+dim.lcl-mm
    j = j+dim.lcl+1-mm
    mm=mm+1
  }
  output = list(-grm.lcl,grc.lcl)
  return(output)
  #/* note that you need a negative in the first returned term
  #be careful with this  */
}

gradcorcov = function(gss.c,litomega,omegastar) { 
  
  litomega = matrix(litomega, ncol = 1)
  gss.c = matrix(gss.c, ncol=1)
  
  dim.lcl = nrow(litomega)
  row.lcl = dim.lcl%**%(dim.lcl+1)/2
  col.lcl = dim.lcl%**%(dim.lcl-1)/2
  
  nu1.lcl=litomega^2
  nu2.lcl=litomega
  temp = gss.diagrv(gss.zeros(dim.lcl,dim.lcl),litomega)
  nuomega.lcl = temp%**%omegastar%**%temp
  
  grm.lcl = gss.zeros(row.lcl,dim.lcl)
  cc1.lcl=0.5%**%((gss.c)%//%(nu1.lcl))
  l.lcl = 0
  for(idimvn in seq.int(1,dim.lcl,1)) {
    grm.lcl[l.lcl+1,idimvn]=cc1.lcl[idimvn, drop = F]
    l.lcl = l.lcl + (dim.lcl+1-idimvn)
  }
  
  grc.lcl = gss.zeros(row.lcl,col.lcl)
  j=1
  mm = 1
  c1=1
  while(!(j== row.lcl  )) {
    ndiag1 = 1%//%(nu2.lcl[mm, drop = F]%**%nu2.lcl[(mm+1):dim.lcl, drop = F])
    ndiag1 = gss.diagrv(diag(dim.lcl-mm),ndiag1)
    nhor1 = (-0.5)%**%(omegastar[mm, (mm+1):dim.lcl]%//%nuomega.lcl[mm,mm, drop = F])
    nhor2 = (-0.5)%**%(omegastar[mm, (mm+1):dim.lcl]%//%(t((nu1.lcl[(mm+1):dim.lcl, drop = F]))))
    grc.lcl[(j+1):(j+dim.lcl-mm),c1:(c1+dim.lcl-mm-1)]=ndiag1
    grc.lcl[j,c1:(c1+dim.lcl-1-mm)]=nhor1
    k=j
    l=1
    while(!(l==dim.lcl+1-mm)) {
      grc.lcl[k+dim.lcl+2-l-mm,l+c1-1]=nhor2[l, drop = F]
      k=k+dim.lcl+2-l-mm
      l=l+1
    }
    c1=c1+dim.lcl-mm
    j = j+dim.lcl+1-mm
    mm=mm+1
  }
  output = list(-grm.lcl,grc.lcl)
  return(output)
  #/* note that you need a negative in the first returned term
  #be careful with this  */
}

test.grad.corcov = function() {
  lit_omega = c(1, 2, 3)
  omegastar = diag(3)
  gradcorcov(c(2, 4, 6), lit_omega, omegastar)
}


gaomegab = function(x1,x2) { 
  l = nrow(x1)
  k = ncol(x1)
  m = nrow(x2)
  n = ncol(x2)
  g = t(x1)%x%x2    
  if(.omsymmetric) {
    temp1 = vecsymmetry(k) 
    g = temp1%**%g
    if(.omdiagonal) {
      tempselmatrix = vecdup(diag(k))            
      g = gss.selif(g,tempselmatrix)
    }
  }
  return(g)
}


multruncbivariate = function(mu,cov,trpoint) { 
  
  mu = matrix(mu)
  trpoint = matrix(trpoint)
  
  m = nrow(mu)
  muuntrunc = mu[1:2, drop = F]
  v11 = cov[1:2,1:2, drop = F]
  v11inv = solve(v11)
  output = bivariatenormaltrunc(muuntrunc,v11,trpoint)  
  mu1 = output[[1]]
  omega11 = output[[2]]
  v12 = cov[1:2,3:m, drop = F]
  v21 = t(v12)
  v22 = cov[3:m,3:m, drop = F]
  mu2 = (mu1)%|%(mu[3:m, drop = F]+((v21%**%v11inv)%**%(mu1-muuntrunc)))     
  omega12 = omega11%**%v11inv%**%v12
  omg = (omega11%~%omega12)%|%((t(omega12))%~%(v22-(v21%**%(v11inv-v11inv%**%omega11%**%v11inv)%**%v12)))  
  output = list(mu2,omg)
  return(output)
}  

#/*****************************************************************************************************************************************************
#                         SET 6: NUMERICAL GRADIENT FOR CDFMVNANALYTIC
#*****************************************************************************************************************************************************/

numgrad.cdfmvnanalytic = function(mu, cov, x, s=.Random.seed, grad.step=10^-8) {
  
  size = length(mu)
  base = cdfmvnanalytic(mu, cov, x, s)[[1]]
  
  mu.grad = mu*0
  for (i in 1:size) {
    mu2 = mu
    mu2[i] = mu2[i] + grad.step
    p.new = cdfmvnanalytic(mu2, cov, x, s)[[1]]
    mu.grad[i] = (p.new - base)/grad.step
  }
  
  x.grad = x*0
  for (i in 1:size) {
    x2 = x
    x2[i] = x2[i] + grad.step
    p.new = cdfmvnanalytic(mu, cov, x2, s)[[1]]
    x.grad[i] = (p.new - base)/grad.step
  }
  
  cov.grad = rep(0, times=size*(size+1)/2)
  count = 0
  for (i in 1:size) {
    for (j in i:size) {
      count = count + 1
      cov2 = cov
      cov2[i, j] = cov2[i, j] + grad.step
      if (i != j) cov2[j, i] = cov2[j, i] + grad.step
      p.new = cdfmvnanalytic(mu, cov2, x, s)[[1]]
      cov.grad[count] = (p.new - base)/grad.step
    }
  }
  return(list(base, matrix(mu.grad), matrix(cov.grad), matrix(x.grad)))
}

numgrad.cdfbvn = function(l1, l2, rho, grad.step = 10^-8) {
  base = gss.cdfbvn(l1, l2, rho)
  gl1 = (gss.cdfbvn(l1+grad.step, l2, rho) - base)/grad.step
  gl2 = (gss.cdfbvn(l1, l2+grad.step, rho) - base)/grad.step
  grho = (gss.cdfbvn(l1, l2, rho+grad.step) - base)/grad.step
  return(list(gl1, gl2, grho))
}

general.numgrad = function(func, ..., const.args=NULL, grad.step=10^-4) {
  inputs = list(...)
  outputs = list()
  base = do.call(func, c(inputs, const.args))
  for (i in 1:length(inputs)) {
    input = inputs[[i]]
    if (is.null(dim(input))) input = matrix(input)
    if (dim(input)[2]==1) {
      g = input*0
      for (j in 1:length(input)) {
        input2 = input
        input2[j] = input2[j] + grad.step
        inputs2 = inputs
        inputs2[[i]] = input2
        g[j] = (do.call(func, c(inputs2, const.args)) - base)/grad.step
      }
    } else {
      g = matrix(rep(0, nrow(input)*(nrow(input)+1)/2))
      count = 0
      for (j in 1:nrow(input)) {
        for (k in j:nrow(input)) {
          count = count+1
          input2 = input
          input2[k, j] = input2[k, j] + grad.step
          if (k != j) input2[j, k] = input2[j, k] + grad.step
          inputs2 = inputs
          inputs2[[i]] = input2
          g[count] = (do.call(func, c(inputs2, const.args)) - base)/grad.step
        }
      }
    }
    outputs[[i]] = g
  }
  return(outputs)
}

seedless.cdfmvnanalytic = function(mu, V, x, s) {
  return(cdfmvnanalytic(mu, V, x, s)[[1]])
}

# general.numgrad(seedless.cdfmvnanalytic, mu, V, x, const.args = list(s=1), grad.step = 0.0001)

# gradnoncdfn(2, 1.5, 3)
# general.numgrad(noncdfn, 2, 1.5, 3)
