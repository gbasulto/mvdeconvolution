
library(mvdeconvolution)
library(reshape2)
library(ggplot2)

## Example 1
## Real and imag. parts of exponential dist.
samp <- rexp(60)
t <- seq(-3, 3, len = 150)
vals <- sapply(t, kde.cf, samp, .1)
real <- Re(vals)
imag <- Im(vals)
## ## Plot values. Uncomment code below
## df <- data.frame(t, real, imag)
## df.molten <- melt(df, id = "t")
## head(df.molten)
## qplot(t, value, data = df.molten, geom = 'line',
##       facets = . ~ variable)

## Example 2
## Real part of bivariate normal emp. ch. fnc.
samp <- cbind(rnorm(100), rnorm(100))
t1 <- seq(-3, 3, len = 150)
t2 <- t1
t <- expand.grid(t1, t2)
ecf.vals <- apply(t, 1, kde.cf, samp, diag(1, 2))
re.vals <- matrix(Re(ecf.vals), 150, 150)
persp(t1, t2, re.vals)
