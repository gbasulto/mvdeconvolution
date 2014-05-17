#' Empirical characteristic function
#' 
#' It computes the d-dimensional empirical characteristic function
#' @param t d-variate vector where the ecf will be evaluated
#' @param samp nxd matrix with the n samples
#' @return A complex value
#' @examples
#' ## Univariate normal
#' samp <- rnorm(50)
#' ## qplot(samp, )
ecf <- function(t, samp)
{
  ## Computes ecf of a bivariate sample (univariate).
  ## It returns a complex number with the value
  ##
  ## Arguments:
  ## t: Argument of the ecf, a bivariate (univ.) vector.
  ## samp: Sample, a 2-columns matrix.
  ##
  ## Output:
  ## A complex (real) number.
  
  samp <- as.matrix(samp)
  prod <- samp %*% t
  return(mean(cos(prod)) + 1i*mean(sin(prod)))
}