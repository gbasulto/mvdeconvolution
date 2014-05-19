
library(mvdeconvolution)

###################################
##  Test relationship between input and output correlations
###################################
#
r1 <- seq(-1, 1, length.out= 20)
ff <- function(r) cor(rgam(2e4, .5, .5, 9, 2, r))[1, 2]
##ff <- function(r) rgam(1e5, .5, 5, .5, 5, r)$correl ## Interesting case
r2 <- sapply(r1, ff)
plot(r2 ~ r1, type = "l", main = "r2 (sample correl. gammas) vs. \n r1 (correl. normal)")
abline(0, 1, col = 2)
legend("topleft", legend = c("Correlations", "Identity"), 
       col = c(1:2), pch = 20)

###################################
## Observe how a biv. gamma looks like
###################################
##
n <- 1e3
shape1 <- 20
rate1 <-  2
shape2 <- 5
rate2 <-  1
r <- -.6

sam <- rgam(n, shape1, rate1, shape2, rate2, r)
colnames(sam) <- c("x", "y")

ran1 <- range(sam[, 1])
ran2 <- range(sam[, 2])
x <- seq(ran1[1], ran1[2], length.out= 100)
y <- seq(ran2[1], ran2[2], length.out= 90)
z <- apply(matrix(y), 1, function(a)
  dgam(cbind(x, a), shape1, rate1, shape2, rate2, r))


plot(y ~ x, data = sam, col = "gray",
     main = "Sample with true density overlaid")
contour(x, y, z, add = T)






