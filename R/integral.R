#' Multivariate Fourier integral
#' 
#' This function computes the a Fourier integral of a function with
#' support in a hyper-rectangle. Right now, this function only
#' computes univariate and bivariate continuous Fourier transform
#' based on the paper 'Fast computation of multidimensional Fourier
#' integrals' by Inverarity (2002). It is the formula (4.1).
#' @param f function from R^2 or R to which the FT will be applied.
#' @param n Dimension of the function f (1 or 2).
#' @param m Resolution of the integral
#' @param a nx1 vector. Lower integration limit.
#' @param b nx1 vector. Upper integration limit.
#' @param c nx1 vector. Lower limit of w.
#' @param d nx1 vector. Upper limit of w.
#' @param r Power in (4.1).
#' @param s Scale constant in (4.1).
#' @examples
#' ##  library(ggplot2)
#' ##
#' ##  Computing characteristic function of
#' ##  univariate normal on -1, 10.
#' ##
#' cf <- FIntegral(dnorm, n = 1, m = 2^12, -20, b = 20, c = -1, d = 10,
#'         r = 1, s = 1)
#' values <- data.frame(t = cf$w,
#'                      cf = c(Re(cf$ft),
#'                             sapply(cf$w, function(t) exp(-t^2/2))),
#'                      type = rep(c("approx", "real"), length(cf$w))
#'                      )
#' ggplot2:::qplot(t, cf, data = values, col = type, geom = "line")
#'
#' ## Real anf imag. parts of the characteristic function of an exponential
#' ## distribution, approximated with 128 points.
#' chf <- function(t) 1/(1 - 1i*t)
#' f <- function(x) dexp(x, 1)
#' cf <- FIntegral(f, n = 1, m = 2^7, a = 0, b = 5, c = -3, d = 10, r = 1, s = 1)
#' ##
#' realEst <- Re(cf$ft)
#' imagEst <- Im(cf$ft)
#' real <- Re(chf(cf$w))
#' imag <- Im(chf(cf$w))
#' m <- length(cf$w)
#' values <- data.frame(w = cf$w,
#'                     val = c(realEst, imagEst, real, imag),
#'                     Part = rep(rep(c("Real", "Imag"), 2), each = m),
#'                     Type = rep(c("FT", "Direct"), each = 2*m))
#' 
#' ggplot2:::qplot(w, val, data = values, geom = "line",
#'       facets = . ~ Part, colour = Type)
#'
#' ## Characteristic function of a bivariate normal distribution
#' chf <- function(t1, t2) exp(-(t1^2 + t2^2)/2)
#' f <- function(x, y) dnorm(x)*dnorm(y)
#' 
#' cf <- FIntegral(f, n = 2, m = 2^8, a = c(-6, -6), b = c(6, 6),
#'         c = c(-3, -3), d = c(3, 3), r = 1, s = 1)
#' 
#' persp(Re(cf$ft), col = "lightblue", phi = 15, theta = 30,
#'       shade = .3, border = NA)
#' @export
FIntegral <- function(f, n, m, a, b, c, d, r, s)
{
    ## Description:
    ## This function computes univariate and bivariate continuous
    ## Fourier tranform based on the paper by Inverarity (2002):
    ## "Fast computation of multidimensional Fourier integrals".
    ## It is the formula (4.1) on the paper.
    ##
    ## Arguments:
    ## f: Function from R^2 or R to C to which we will apply
    ##    the ft.
    ## n: Dimension of the function above.
    ## m: Resolution of the integral.
    ## a: nx1 vector. Lower integration limit.
    ## b: nx1 vector. Upper integration limit.
    ## d: nx1 vector. Lower limit of w.
    ## l: nx1 vector. Upper limit of w.
    ## r: Power in (4.1).
    ## s: Scale constant in (4.1).
    ##
    ## Output:
    ## w: vector or matrix with the values for which the cft was computed.
    ## ft: Continuous Fourier transform values at w.
    ## 

    ## This is an adjustment for the upper limit:
    d <- c + m*(d - c)/(m - 1)
    
    ## r = 1 is equivalent to the following:
    if(s != 1)
        {
            out <- FIintegral(f, n, m, a, b, s*c,
                     s*d, r, 1)
            w <- out$w/s
            return(list(w = w,
                        ft = abs(s)^(n/2)*out$ft))
        }
    

    if(n == 1)
        {
            ##  The next two lines are there because the code below
            ##  computes the integral with negative sign. By adding this
            ##  two lines, we compute the univ. formula 4.1
            c <- -c
            d <- -d
            
            bet <- (b - a)/m
            gam <- (d - c)/m
            del <- bet*gam/2
            J1 <- 0:(m - 1)
            J2 <- m:(2*m - 1)
            t <- a + bet*J1
            w <- c + gam*J1
            y <- c(f(t)*complex(argument = -J1*(bet*c + del*J1)),
                   rep(0, m))
            z <- complex(argument = del*(c(J1^2, (J2 - 2*m)^2)))
            val <- bet*complex(argument = -(a*w + del*J1^2))*
                fft(fft(y)*fft(z), inverse = T)/
                    (2*pi)^((1 - r)/2)/(2*m)

            ##  ... The same with this line.
            w <- -w
            
            return(list(w = w,
                        ft = val[J1 + 1]))
        }
    
    ##  nx1 vectors.
    bet <- (b - a)/m
    gam <- (d - c)/m
    del <- bet*gam/2
    a_hat <- a + bet/2
    
    ##  Aux. mx1 vectors
    J1 <- 0:(m - 1)
    J2 <- m:(2*m - 1)
    
    ## nxm matrices
    t <- sweep(bet %o% J1, 1, a_hat, "+")

    w <- sweep(gam %o% J1, 1, c, "+")
    
    ## nx2m matrix
    auxArg <- cbind(- del %o% (J1^2), - del %o% (J2 - 2*m)^2)
    z <- exp(1i*auxArg)
 ##   cat(dim(z))
 ##   cat("\n")

    
    ## Starting here, the program will work only for 2-dim. ft.

    ## m x m matrices
                                        # Exponential in 4.4, mxm
    auxArg <- outer(J1, J1,
                    function(j1, j2) j1*(bet[1]*c[1] + del[1]*j1)
                    + j2*(bet[2]*c[2] + del[2]*j2))
    aux1 <- exp(1i*auxArg)
                                        #     f(t) in 4.4, mxm
    aux2 <- apply(matrix(t[2, ]), 1, function(y)
                  apply(matrix(t[1, ]), 1, function(x) f(x, y)))

    ## 2m x 2m matrix
                                        #  y in 4.4, first filled out
                                        #  with zeros.
    y <- matrix(0, 2*m, 2*m)
    y[1:m, 1:m] <- aux1*aux2
    
    ##    Univariate dft, mx1 vectors.
    dft1 <- drop(fft(z[1, ]))
    dft2 <- drop(fft(z[2, ]))

    ##    Values to apply inverse dft in 4.8: 2m x 2m matrix.
    dft <- fft(y) * (dft1 %o% dft2)

    ##   mxm matrix of exponentials in 4.8
    aux1 <- drop(a_hat[1]*w[1, ] + del[1]*J1^2)
    aux2 <- drop(a_hat[2]*w[2, ] + del[2]*J1^2)
#
    expo <- complex(argument= (aux1 %o% aux2)) # mxm
    expo <- exp(1i*outer(aux1, aux2, '+')) # mxm
    fact <- prod(bet)*((2*pi)^(1 - r))^(-n/2) # real
    idft <- (fft(dft, inverse = T)/(2*m)^2)[1:m, 1:m]
    
    ##    FT
    val <- expo*fact*idft
    
    return(list(w = w,
                ft = val))
}
