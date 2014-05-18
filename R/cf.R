#' Empirical characteristic function
#' 
#' It computes the d-dimensional empirical characteristic function
#' @param t d-variate vector where the ecf will be evaluated
#' @param samp nxd matrix with the n samples
#' @return A complex value
#' @example
#' examples/ecf_ex.R
#' @export
ecf <- function(t, samp)
{
  samp <- as.matrix(samp)
  prod <- samp %*% t
  return(mean(cos(prod)) + 1i*mean(sin(prod)))
}

#' Characteristic function of multivatiate normal kde
#' 
#' Characteristic function of the multivariate kernel density estimator,
#' assuming that 
#' @param t d-variate vector where the ecf will be evaluated
#' @param samp nxd matrix with the n samples
#' @param Sigma is the bandwidth matrix (dxd).
#' @return A complex value
#' @example
#' examples/kde.cf_ex.R
#' @export
kde.cf <- function(t, samp, Sigma)
{
  ##
  ## Description:
  ## It returns the characteristic function of the
  ## bivariate normal kernel for a single point.
  
  samp <- as.matrix(samp)
  modu <- exp(- .5*sum(t * (Sigma %*% t)))
  arg <- samp %*% t
  real <- modu*mean(cos(arg))
  imag <- modu*mean(sin(arg))
  
  return(real + 1i*imag)
}