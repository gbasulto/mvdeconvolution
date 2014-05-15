
library(ggplot2)

source("integral.R")


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

qplot(w, val, data = values, geom = "line",
      facets = . ~ Part, colour = Type)



## Bivariate
chf <- function(t1, t2) exp(-(t1^2 + t2^2)/2)
f <- function(x, y) dnorm(x)*dnorm(y)

system.time({bb <- I(f, n = 2, m = 2^8, a = c(-6, -6), b = c(6, 6),
        c = c(-3, -3), d = c(3, 3), r = 1, s = 1)})

persp(Re(bb$ft), col = "lightblue", phi = 15, theta = 30,
      shade = .3, border = NA)
#     Comparing with original function.
vals <- outer(c(bb$w[1, ]), c(bb$w[2, ]), chf)
##persp(Mod(bb$ft - vals), col = "lightblue", phi = 15, theta = 30,
##      shade = .3, border = NA)
hist(Mod(bb$ft - vals))
hist(Re(bb$ft - vals))
hist(Im(bb$ft - vals))

