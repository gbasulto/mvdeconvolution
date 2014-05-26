
rm(list = ls())

#' Randon sample of bivariate gamma
#' 
#' The correlated gammas where generated from a
#' bivariate normal whose variance matrix has 1 
#' in the diagonal and r otherwise, then each coordinate
#' was used to generate a gamma distribution.
#' @param n sample size.
#' @param shape1 Shape parameter for the first coordinate.
#' @param rate1 Rate parameter for the first coordinate.
#' @param shape2 Shape parameter for the second coordinate.
#' @param rate2 Rate parameter for the second coordinate.
#' @param r correlation of the normal entries.
#' @return A nx2 matrix with the samples.
#' @seealso dgam
#' @example
#' examples/bivgamma_ex.R
#' @export
rgam <- function(n, shape1, rate1, shape2, rate2, r)
  {
  Sig <- matrix(c(1, r, r, 1), 2, 2)
  Mu <- c(0, 0)
  norm <- rmvnorm(n, Mu, Sig)
  unif <- pnorm(norm)
  gam <- cbind(qgamma(unif[, 1], shape= shape1, rate= rate1), 
               qgamma(unif[, 2], shape= shape2, rate= rate2))
  return(gam)
}

#' Density of bivariate gamma
#' 
#' The correlated gammas where generated from a
#' bivariate normal whose variance matrix has 1 
#' in the diagonal and r otherwise, then each coordinate
#' was used to generate a gamma distribution.
#' @param val 2-columns matrix with the values to evaluate 
#' the density
#' @param shape1 Shape parameter for the first coordinate.
#' @param rate1 Rate parameter for the first coordinate.
#' @param shape2 Shape parameter for the second coordinate.
#' @param rate2 Rate parameter for the second coordinate.
#' @param r correlation of the normal entries.
#' @return A nx2 matrix with the samples.
#' @seealso rgam
#' @example
#' examples/bivgamma_ex.R
#' @export
dgam <- function(val, shape1, rate1, shape2, rate2, r){
  Sig <- matrix(c(1, r, r, 1), 2, 2)
  Mu <- c(0, 0)
  
  if(is.null(dim(val))) dim(val) <- c(1, 2)
  g1 <- val[, 1]
  g2 <- val[, 2]
  
  ##   Gamma cdfs
  FG1 <- pgamma(g1, shape= shape1, rate= rate1)
  FG2 <- pgamma(g2, shape= shape2, rate= rate2)
  ##   Gamma densities
  dG1 <- dgamma(g1, shape= shape1, rate= rate1)
  dG2 <- dgamma(g2, shape= shape2, rate= rate2)
  ##   Aux. quantities
  QFG1 <- qnorm(FG1)
  QFG2 <- qnorm(FG2)
  
  
  dvals <- dG1/dnorm(QFG1)*dG1/dnorm(QFG1)*
    dmvnorm(cbind(QFG1, QFG2), Mu, Sig)

  return(dvals)
}

## # Rotated gamma density
## # 
## # Function that computes the density of the rotated 
## # independent gammas (X,Y). The rotation is given in  
## # the matrix inv_rot which is the inverse
## # of the rotation matrix
## f_rotgamma <- function(xy, a1, b1, a2, b2, inv_rot)
##   {
##     rxy <- inv_rot %*% t(xy)
##     return(dgamma(rxy[1,], a1, b1)*dgamma(rxy[2, ], a2, b2))
## }
## 
## sample_f_rotgamma <- function(n, a1, b1, a2, b2, rot) {
##     ##Function to sample from the rotated gamma
##     ##rot is the rotation matrix
## 
##     s <- matrix(c(rgamma(n, a1, b1), rgamma(n, a2, b2)), ncol= 2)
##     return(t(rot%*%t(s)))
## }






