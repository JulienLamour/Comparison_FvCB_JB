### Simulation of electron transport response to irradiance with various curvature factor (Theta), including negative ones

J <- function(Q ,alpha , phi, theta, Jmax=100){
  if(theta==0){J=Q*alpha*phi*Jmax/(Q*alpha*phi+Jmax)}else{J=(Q*alpha*phi+Jmax-((Q*alpha*phi+Jmax)^2-4*theta*Q*alpha*phi*Jmax)^0.5)/(2*theta)}
}

Q=0:2000
alpha=0.85
phi=0.425
Jmax=100
jpeg(filename = 'Figures/Extended_data_figure_1.jpeg',width = 130,height = 130,units = 'mm',res=300)
plot(x=NULL,y=NULL,xlim=range(Q*1.05),ylim=c(0,1.3*Jmax),ylab=expression(italic(J)),xlab=expression(italic(Q)))
abline(c(0,alpha*phi),lwd=5)
abline(h = Jmax,lwd=5)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=1,Jmax=100),col="grey",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.9,Jmax=100),col="grey",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.7,Jmax=100),col="grey",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0,Jmax=100),col="grey",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=-2,Jmax=100),col="grey",lwd=2)
dev.off()