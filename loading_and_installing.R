
library(devtools)
library(roxygen2)

has_devel()

document()

install()

## Create vignettes
build_vignettes()

## Check
##check()

## Generate compressed source to be sent to CRAN
## build()

library(mvdeconvolution)
library(help = mvdeconvolution)




mvdeconvolution:::FIntegral
FIntegral
## FIntegral()
? FIntegral
? ecf


remove.packages("mvdeconvolution")


## Install from Github
install_github("gbasulto/mvdeconvolution")
