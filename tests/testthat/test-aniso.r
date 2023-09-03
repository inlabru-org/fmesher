library(devtools)
setwd("C:/Users/liaml/OneDrive/Documentos/")
load_all("fmesher")


# Source each file
lapply(files, source)

mesh <- fm_rcdt_2d_inla(globe = 1)
plot(mesh)

kappa <- function(xi) {
  return(sum(xi^2))
}

vec <- function(xi) {
  return(c(xi[1] + 1, xi[2] - 2))
}

nodes <- mesh$loc

kappa_values <- apply(nodes, 1, kappa)
vec_values <- t(apply(nodes, 1, vec))
aniso=list(kappa_values,vec_values)
a<-fm_fem_aniso(mesh,aniso)

