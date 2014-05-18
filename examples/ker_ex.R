
library(mvdeconvolution)
library(reshape2)
library(ggplot2)

## Example for the univariate case
## Plot all the possible kernels
t <- seq(-1.2, 1.2, length.out = 150)
sinc <- sapply(t, ker, 1)
VP <- sapply(t, ker, 2)
triw <- sapply(t, ker, 3)
tric <- sapply(t, ker, 4)
flat <- sapply(t, ker, 5)

## ## Plot. Uncomment code below
## df <- data.frame(t, sinc, VP, triw, tric, flat)
## df.molten <- melt(df, id = "t", variable.name = "kernel")
## qplot(t, value, data = df.molten,
##       geom = "line", colour = kernel)

