
library(devtools)
library(roxygen2)
has_devel()

document()

install()

library(mvdeconvolution)
library(help = mvdeconvolution)

mvdeconvolution:::FIntegral
FIntegral
## FIntegral()
? FIntegral
## Create vignettes
build_vignettes()


remove.packages("mvdeconvolution")


