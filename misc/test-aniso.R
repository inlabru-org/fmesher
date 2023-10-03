library(devtools)
setwd("C:/Users/liaml/OneDrive/Documentos/")
load_all("fmesher")
library(ggplot2)
#PARAMETER VALUES NONSTATIONARY
kp=1
v1=0
v2=5
lambda=exp(sqrt(v1^2+v2^2))
stretch=sqrt(lambda)
#PARAMETER VALUES STATIONARY
nu=2-2/2
rh=sqrt(8*nu)/kp


#KAPPA FIELD
kappa <- function(x) {
  #return((x[1]-2)^2+(x[2]-2)^2) #Correlation length becomes small, lots of variation
  return(kp)
}

vec <- function(x) {
  #return(c(0, 0))
  #return(c(100, 0.1)) #Gives very small field, how come?
  return(c(x[1]-5, x[2]-5)) #This gives large correlation (slow variation) in x coord and faster in y which is correct.
                  #Unsure about why closer to zero though. Also variation in y direction should be much faster
}
#MESH FOR STATTIONAY FIELD
l=10
library(sf)
boundary_sf = st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10),c(0,0)))))
boundary = fm_as_segm(boundary_sf)

mesh1 <- fm_mesh_2d_inla(
  boundary = boundary,
  max.edge = c(0.2, 1))
nodes1 <- mesh1$loc
plot(mesh1)
# #MESH FOR NONSTATIONARY FIELD
# mesh2 <- fm_mesh_2d_inla(
#   boundary = fm_segm(rbind(c(0, 0), c(l*stretch, 0), c(l*stretch, l/stretch), c(0, l/stretch)), is.bnd = TRUE),
#   max.edge = c(stretch*0.5, 0.5))
# #plot(mesh2)
# nodes2 <- mesh2$loc

kappa_values <- apply(nodes1, 1, kappa)
vec_values <- t(apply(nodes1, 1, vec))
aniso=list(kappa_values,vec_values)
#SAMPLE OF MATERN WITH CORRELATION RANGE GIVEN BY KAPPA
samplebase= fm_matern_sample(mesh1,alpha = 2,rho = rh,sigma=1)
#SAMPLE OF ANISOTROPIC FIELD SHOULD BE EQUAL TO MATERN DEFORMED BY u(H^{-1/2}x)
 sampleaniso=fm_aniso_sample(mesh1,aniso)
samplebase= fm_matern_sample(mesh1,alpha = 2,rho = rh,sigma=1)
#SAMPLE OF ANISOTROPIC FIELD SHOULD BE EQUAL TO MATERN DEFORMED BY u(H^{-1/2}x)
# kappa_values2 <- apply(nodes1, 1, kappa)
# vec_values2 <- t(apply(nodes1, 1, vec))
# aniso2=list(kappa_values2,vec_values2)
# # sampleaniso2=fm_aniso_sample(mesh1,aniso2)


database<- data.frame(x = nodes1[,1],
                      y = nodes1[,2],
                      u = samplebase)
dataaniso <- data.frame(x = nodes1[,1],
                   y = nodes1[,2],
                   u = sampleaniso)
# dataaniso2 <- data.frame(x = nodes1[,1],
#                         y = nodes1[,2],
#                         u = sampleaniso2)
pxl = fm_pixels(mesh1, mask= boundary_sf)
pxl$uisotropic <- fm_evaluate(mesh1,loc = pxl,
                              field= database$u)
pxl$uaniso <- fm_evaluate(mesh1,loc = pxl,
                              field= dataaniso$u)


ggplot(pxl)+
  geom_tile(aes(geometry=geometry,fill =uisotropic),
            stat = "sf_coordinates", alpha=1) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"), limits=c(min(database$u), max(database$u))) +
  coord_equal()+
  xlab("X Coordinate") + ylab("Y Coordinate") +
  labs(fill = "Field u isotropic") +
  theme(legend.position = "bottom")

  ggplot(pxl)+
    geom_tile(aes(geometry=geometry,fill =uaniso),
              stat = "sf_coordinates", alpha=1) +
                scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF"), limits=c(min(dataaniso$u), max(dataaniso$u))) +
                coord_equal()+
                xlab("X Coordinate") + ylab("Y Coordinate") +
                labs(fill = "Field u anisotropic") +
                theme(legend.position = "bottom")





