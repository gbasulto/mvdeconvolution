#' Integrand for Fourier transform
#' 
#' This is the integrand in the bivariate deconvolution formula, except
#' for the exp(-itx) term.
#' @param t A bivariate vector.
#' @param EMethod Method used to approximate the characteristic function of 
#' the error. 1 or 2 if ecf of the char. fnc. of the kde is used, respectively.
#' @param error Matrix with 2 columns, either with the pure sample of the error, or with the 
#' differences obtained from the panel data structure.
#' @param sigE Bandwidth (2x2) matrix for the error. It is NULL by default 
#' and it is required only when EMethod = 2.
#' @param truncate NULL by default. If the approximated characteristic function
#' of the error is smaller than truncate, it will make the integrand zero.
#' This provides numerical stability.
#' @param h Bivariate vector with the bandwidth parameters.
#' @param samp 2-column matrix with the contaminated sample.
#' @param columns Number of columns of the panel data. If there is one column, 
#' it assumes that instead panel data, a pure sample of the error is provided.
#' @param Kkernel Kernel to be used. See 'ker' function.
#' @return A complex value.
#' @seealso kerdecon
integrand <- function(t, EMethod, error, sigE, truncate, h,
                      samp, columns, Kkernel)
{
  ##   phiE is modified according to if it is panel data structure, or it 
  ## there is a pure sample of the error.
  if(columns == 1)
  {
    phiE <- (switch(EMethod,
                    ecf(t= t, samp= error),
                    kde.cf(t= t, samp= error, Sigma= sigE)))
  } else
  {
    t_aux <- t/columns
    phiDelta <- (switch(EMethod,
                        ecf(t= t_aux, samp= error),
                        kde.cf(t= t_aux, samp= error,
                               Sigma= sigE))
    )
    phiE <- sqrt(max(0, Re(phiDelta)))^columns
  }
  
  ## If cf of the error is too small, it as  signs the value 0.
  if(!is.null(truncate)) if(Mod(phiE) < truncate) return(0)
  
  s <- h * t
  abs_s <- abs(s)
  
  phiK <- switch(Kkernel,
                 all(abs_s < 1),  # sinc kernel
                 prod(1 - abs_s)*all(abs_s < 1),  #VP kernel
                 prod(1 - abs_s^2)^3*all(abs_s < 1), #Triweight
                 prod(1 - abs_s^3)^3*all(abs_s < 1), #Tricube
                 prod( flat_Ker(abs_s) )*all(abs_s < 1) ) #Flat Kernel
  
  
  phiY <- ecf(t= t, samp= samp)
  
  return(phiY*phiK/phiE)
}

#' Bivariate kernel deconvolution estimator
#' 
#' It computes a bivariate kernel deconvolution estimator using the Fast
#' Fourier transform (FFT) on a grid.
#' @param resol An integer specifying the grid size on each coordinate. Power of 
#' two are faster.
#' @param samp 2-column matrix with the contaminated sample.
#' @param error Matrix with 2 columns, either with the pure sample of the error, or with the 
#' @param truncate NULL by default. If the approximated characteristic function
#' of the error is smaller than truncate, it will make the integrand zero.
#' This provides numerical stability.
#' @param h Bivariate vector with the bandwidth parameters.
#' @param EMethod Method used to approximate the characteristic function of 
#' the error. 1 or 2 if ecf of the char. fnc. of the kde is used, respectively.
#' differences obtained from the panel data structure. It can also accept "ecf" or "kde"
#' respectively.
#' @param columns Number of columns of the panel data. If there is one column, 
#' it assumes that instead panel data, a pure sample of the error is provided.
#' @param Kkernel Kernel to be used. See 'ker' function.
#' @param coord1Range A bivariate vector with the limits of first coordinate of
#' the grid to estimate the density.
#' @param coord2Range A bivariate vector with the limits of second coordinate of
#' the grid to estimate the density.
#' @param sigE Bandwidth (2x2) matrix for the error. It is NULL by default 
#' and it is required only when EMethod = 2.
#' @param pstve TRUE/FALSE determining whether returning zero for those 
#' points where the estimate in negative.
#' @return A list with the following elements:
#' @param x1 Vector of size resol with the first coordinate grid.
#' @param x2 Vector of size resol with the second coordinate grid.
#' @param z resol x resol grid with the estimated values.
#' @export
kerdecon <- function(resol, samp, error, truncate, h, EMethod,
                       columns, Kkernel,
                       coord1Range, coord2Range,
                       sigE = NULL, pstve = TRUE){  
  ## Change EMethod value to numeric, if necessary.
  if(!is.numeric(EMethod)) EMethod <- ifelse(EMethod == "ecf", 1, 2)
  
  f <- function(t1, t2) integrand(t= c(t1, t2), EMethod = EMethod,
                                  error = error,
                                  sigE = sigE,
                                  truncate = truncate,
                                  h= h, samp = samp,
                                  columns = columns,
                                  Kkernel)
  
  res <- FIntegral(f= f, n= 2, m= resol, a= - 1/h, b= 1/h,
           c= c(coord1Range[1], coord2Range[1]),
           d= c(coord1Range[2], coord2Range[2]), r = -1, s = -1)
  
  ##       Computing proportion of negative values and changing those
  ##       by zero.
  z <- Re(res$ft)
  ## cat(paste("Prop. of neg. vals. =", mean(z < 0), "\n"))
  if(pstve) z[z < 0] <- 0
  
  return(list(x1 = res$w[1, ],
              x2 = res$w[2, ],
              z = z))
}
