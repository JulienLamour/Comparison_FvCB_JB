### Simulation of electron transport response to irradiance with various curvature factor (Theta), including negative ones

J <- function(Q ,alpha , phi, theta, Jmax=100){
  if(theta==0){J=Q*alpha*phi*Jmax/(Q*alpha*phi+Jmax)}else{J=(Q*alpha*phi+Jmax-((Q*alpha*phi+Jmax)^2-4*theta*Q*alpha*phi*Jmax)^0.5)/(2*theta)}
}

J_smith <- function(Q ,alpha , phi, Jmax=100){
  Q*alpha*phi/sqrt(1+(Q*alpha*phi/Jmax)^2)
}

J_lin <- function(Q, alpha, phi){
  Q*alpha*phi
}

Q=0:2000
alpha=0.85
phi=0.425
Jmax=100

jpeg(filename = 'Figures/Suppementary_Figure_NoteS1.jpeg',width = 150,height = 130,units = 'mm',res=300)
par(mar=c(5.1,5.1,2.1,2.1))
plot(x=NULL,y=NULL,xlim=range(Q),ylim=c(0,1.1*Jmax),ylab=expression(Electron~transport~rate~italic(J)~(mu*mol~m^-2~s^-1)),xlab=expression(Photosynthetic~active~radiation~(mu*mol~m^-2~s^-1)))
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0,Jmax=100),col="deeppink",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.7,Jmax=100),col="deepskyblue",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.9,Jmax=100),col="deepskyblue",lwd=2)
lines(x=Q,y=J_smith(Q=Q,alpha=alpha,phi=phi,Jmax=100),col="darkseagreen",lwd=2) 
lines(x=Q,y=J_lin(Q=Q,alpha=alpha,phi=phi),col="grey",lwd=2) 
abline(c(0,alpha*phi),lwd=2,lty=3)
abline(h = Jmax,lwd=2,lty=3)

text(x = 1000, y = 75, labels = "Eqn S4", cex = 1, col = "deeppink", srt = 15)
text(x = 1000, y = 87, labels = expression(Eqn~S1*','~italic(theta)~'='~0.7), cex = 1, col = "deepskyblue", srt = 10)
text(x = 510, y = 95, labels = expression(Eqn~S1*','~italic(theta)~'='~0.9), cex = 1, col = "deepskyblue", srt = 25)
text(x = 680, y = 89, labels = expression(Eqn~S2), cex = 1, col = "darkseagreen", srt = 15)
text(x = 115, y = 60, labels = expression(Eqn~S3), cex = 1, col = "grey", srt = 80)
dev.off()

#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length = n + 1)
#  hcl(h = hues, l = 65, c = 100)[1:n]
#}
#
#ggplot(data=result,aes(x=Q,y=J,group=Theta))+
#  geom_line(color=gg_color_hue(2)[2])+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank())+
#  geom_abline(slope=alpha*phi,linetype = "dashed")+geom_hline(yintercept = Jmax,linetype="dashed")+
#  ylab(expression(Electron~transport~rate~italic(J)~(mu*mol~m^-2~s^-1)))+xlab(expression(Irradiance~italic(Q)~(mu*mol~m^-2~s^-1)))
  
