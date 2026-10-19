### Simulation of electron transport response to irradiance with various curvature factor (Theta), including negative ones

J <- function(Q ,alpha , phi, theta, Jmax=100){
  if(theta==0){J=Q*alpha*phi*Jmax/(Q*alpha*phi+Jmax)}else{J=(Q*alpha*phi+Jmax-((Q*alpha*phi+Jmax)^2-4*theta*Q*alpha*phi*Jmax)^0.5)/(2*theta)}
}

Q=0:2000
alpha=0.85
phi=0.425
Jmax=100

result=data.frame()
for (Theta in c(0,0.7,0.9,1)){
  result=rbind.data.frame(result,data.frame(Q=Q,J=J(Q=Q,alpha=alpha,phi=phi,theta=Theta,Jmax=Jmax),Theta = Theta))
}

jpeg(filename = 'Figures/Extended_data_figure_1.jpeg',width = 150,height = 130,units = 'mm',res=300)
par(mar=c(5.1,5.1,2.1,2.1))
plot(x=NULL,y=NULL,xlim=range(Q*1.05),ylim=c(0,1.3*Jmax),ylab=expression(Electron~transport~rate~italic(J)~(mu*mol~m^-2~s^-1)),xlab=expression(Photosynthetic~active~radiation~(mu*mol~m^-2~s^-1)))
abline(c(0,alpha*phi),lwd=5,lty=3)
abline(h = Jmax,lwd=5,lty=3)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=1,Jmax=100),col="deepskyblue",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.9,Jmax=100),col="deepskyblue",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0.7,Jmax=100),col="deepskyblue",lwd=2)
lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=0,Jmax=100),col="deepskyblue",lwd=2)
#lines(x=Q,y=J(Q=Q,alpha=alpha,phi=phi,theta=-2,Jmax=100),col="grey",lwd=2)
dev.off()

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggplot(data=result,aes(x=Q,y=J,group=Theta))+
  geom_line(color=gg_color_hue(2)[2])+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank())+
  geom_abline(slope=alpha*phi,linetype = "dashed")+geom_hline(yintercept = Jmax,linetype="dashed")+
  ylab(expression(Electron~transport~rate~italic(J)~(mu*mol~m^-2~s^-1)))+xlab(expression(Irradiance~italic(Q)~(mu*mol~m^-2~s^-1)))
  
