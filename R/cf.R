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
  samp <- as.matrix(samp)
  modu <- exp(- .5*sum(t * (Sigma %*% t)))
  arg <- samp %*% t
  real <- modu*mean(cos(arg))
  imag <- modu*mean(sin(arg))
  
  return(real + 1i*imag)
}

#' FT of product kernels
#' 
#' It computes the Fourier transform of several kernels 
#' with support in -1 to 1.
#' @param t d-dimensional vector to evaluate the FT.
#' @param kernel Value from 1 to 5. 1, sinc kernel; 
#' 2, VP kernel;
#' 3, kernel whose FT is prop. to triweight kernel;
#' 4, kernel whose FT is prop. to Tricube kernel, and
#' 5, flat-top kernel.
#' @example
#' examples/ker_ex.R
#' @export
ker <- function(t, kernel)
{
  ## Dimension
  d <- length(t)
  ## Absolute value
  at <- abs(t)
  ## If any is greater than 1 in abs. value, it returns 0.
  if(any(at > 1)) return(0)
  
  phiK <- switch(kernel,
                 1,  # sinc kernel
                 prod(1 - at),  #VP kernel
                 prod(1 - at^2)^3, #Triweight
                 prod(1 - at^3)^3, #Tricube
                 prod((1 - 2*(at - 0.5))^(at >= 0.5)) #Flat Kernel
  ) 
  return(phiK)
}
