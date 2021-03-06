% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/integral.R
\name{FIntegral}
\alias{FIntegral}
\title{Multivariate Fourier integral}
\usage{
FIntegral(f, n, m, a, b, c, d, r, s)
}
\arguments{
\item{f}{function from R^2 or R to which the FT will be applied.}

\item{n}{Dimension of the function f (1 or 2).}

\item{m}{Resolution of the integral}

\item{a}{nx1 vector. Lower integration limit.}

\item{b}{nx1 vector. Upper integration limit.}

\item{c}{nx1 vector. Lower limit of w.}

\item{d}{nx1 vector. Upper limit of w.}

\item{r}{Power in (4.1).}

\item{s}{Scale constant in (4.1).}
}
\description{
This function computes the a Fourier integral of a function with
support in a hyper-rectangle. Right now, this function only
computes univariate and bivariate continuous Fourier transform
based on the paper 'Fast computation of multidimensional Fourier
integrals' by Inverarity (2002). It is the formula (4.1).
}
\examples{
##  library(ggplot2)
##
##  Computing characteristic function of
##  univariate normal on -1, 10.
##
cf <- FIntegral(dnorm, n = 1, m = 2^12, -20, b = 20, c = -1, d = 10,
        r = 1, s = 1)
values <- data.frame(t = cf$w,
                     cf = c(Re(cf$ft),
                            sapply(cf$w, function(t) exp(-t^2/2))),
                     type = rep(c("approx", "real"), length(cf$w))
                     )
ggplot2:::qplot(t, cf, data = values, col = type, geom = "line")

## Real anf imag. parts of the characteristic function of an exponential
## distribution, approximated with 128 points.
chf <- function(t) 1/(1 - 1i*t)
f <- function(x) dexp(x, 1)
cf <- FIntegral(f, n = 1, m = 2^7, a = 0, b = 5, c = -3, d = 10, r = 1, s = 1)
##
realEst <- Re(cf$ft)
imagEst <- Im(cf$ft)
real <- Re(chf(cf$w))
imag <- Im(chf(cf$w))
m <- length(cf$w)
values <- data.frame(w = cf$w,
                    val = c(realEst, imagEst, real, imag),
                    Part = rep(rep(c("Real", "Imag"), 2), each = m),
                    Type = rep(c("FT", "Direct"), each = 2*m))

ggplot2:::qplot(w, val, data = values, geom = "line",
      facets = . ~ Part, colour = Type)

## Characteristic function of a bivariate normal distribution
chf <- function(t1, t2) exp(-(t1^2 + t2^2)/2)
f <- function(x, y) dnorm(x)*dnorm(y)

cf <- FIntegral(f, n = 2, m = 2^8, a = c(-6, -6), b = c(6, 6),
        c = c(-3, -3), d = c(3, 3), r = 1, s = 1)

persp(Re(cf$ft), col = "lightblue", phi = 15, theta = 30,
      shade = .3, border = NA)
}

