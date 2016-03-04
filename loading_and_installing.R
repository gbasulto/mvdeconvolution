
rm(list = ls())
library(devtools)
library(roxygen2)
has_devel()

install()

## Create vignettes
build_vignettes()

library(mvdeconvolution)
library(help = mvdeconvolution)

? FIntegral
? ecf
? kde.cf
? ker
? rgam
? dgam
? kerdecon


# remove.packages("mvdeconvolution")


## Install from Github
install_github("gbasulto/mvdeconvolution")
