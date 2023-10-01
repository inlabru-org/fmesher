library(devtools)
load_all("fmesher")
mesh <- fm_rcdt_2d(loc = rbind(c(0,0),c(1,0),c(1,1),c(0,1)),
                   tv = rbind(c(1,2,3),c(3,4,1)))
plot(mesh)
nodes1=mesh$loc

kp=1
v1=1
v2=0
vnorm=sqrt(v1^2+v2^2)
lambda=exp(vnorm)
stretch=sqrt(lambda)

mesh2 <- mesh
mesh2$loc <- mesh2$loc%*% diag(c(1/stretch,stretch,1))
#PARAMETER VALUES STATIONARY
nu=2-2/2
rh=sqrt(8*nu)/kp

#PARAMETERS OLD ANISO
gamma= sqrt(cosh(vnorm)-sqrt(cosh(vnorm)^2-1))
u1=sqrt(sinh(vnorm)+sinh(vnorm)*v1/vnorm)
u2 =sinh(vnorm)*v2/(vnorm*u1)
unorm=sqrt(u1^2+u2^2)
gamma_=sqrt(-unorm^2/2+sqrt(unorm^4/4+1))
hfugl=gamma_^2*diag(2)+t(t(c(u1,u2)))%*%t(c(u1,u2))
Hus= cosh(vnorm)*diag(2)+sinh(vnorm)/vnorm *matrix(c(v1,v2,v2,-v1),2,2)

print(rbind(hfugl,Hus))


#KAPPA FIELD
kappa <- function(x) {
  #return((x[1]+2)^2+(x[2]+2)^2) Correlation length becomes small, lots of variation
  return(kp)
}

vec <- function(x) {
  return(c(v1, v2)) #This gives large correlation (slow variation) in x coord and faster in y which is correct.
  #Unsure about why closer to zero though. Also variation in y direction should be much faster
}
kappa_values <- apply(nodes1, 1, kappa)
vec_values <- t(apply(nodes1, 1, vec))
aniso=list(kappa_values,vec_values)

a = fm_fem(mesh2)
b = fm_fem_aniso(mesh,aniso)
C = fm_fem(mesh,aniso= list(gamma=rep(gamma,mesh$n),vec=matrix(c(u1,u2,0),mesh$n,3,byrow=TRUE)))
print(cbind(a$c0,b$c0))
print(cbind(a$c1,b$c1))
print(cbind(a$g1,b$g1))
print(cbind(a$g2,b$g2))
print(rbind(a$g1,b$g1,C$g1aniso))
print((a$g1+C$g1aniso)/2)


