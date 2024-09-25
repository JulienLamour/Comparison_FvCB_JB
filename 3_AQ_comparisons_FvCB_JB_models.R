#########################################################################
###   In this code, I do simulations with the JB and FvCB models to   ###
###   analyse the response to light of both models.                   ###
###   First, I show that the convexity of a light curve (An ~ Q) is   ###
###   different to the convexity of J ~ Q. As a consequence,          ###
###   equation 5 can't be used to estimate Theta in the FvCB model.   ###
###   Second, I show that the convexity of the light curves produced  ###
###   by the JB model is much lower than that of the FvCB model.      ###
###   Mathematically, the JB model is a special case of the FvCB      ###
###   model for Theta = 0 if eta = 1. Here, I show that when we       ###
###   consider the effect of Ci on eta, the apparent convexity of     ###
###   the JB model is indeed close to zero                            ###
###   Finally, I show that if Theta is not 0, then Vqmax is different ###
###   to Jmax.                                                        ###
#########################################################################

library(ggplot2)
library(bbmle)
library(here)
path=here()
setwd(path)
source('R models/Farquhar_von_Caemmerer_Berry_photosynthesis_model.R')

### Simulation of an AQ curve using the FvCB model showing an Ac limitation at high light
### This simulations shows that the shape of the light curve is not only the result of
### the parameter Theta (curvature of electron transport) but also the result of the transition
### from Aj to Ac. This makes the use of Long equation unsuitable to estimate Theta in the FvCB model at ambient CO2.
### Many studies use this empirical equation to derive theta but this is incorrect, as discussed by 
### (Ziegler-Jöns & Selinger, 1987; Ögren & Evans, 1993).
jpeg(filename = 'Figures/Other_Effect_Ac_curvature.jpeg',width = 120,height = 120,units = 'mm',res=300)
par(mar=c(5,5,1,1))
VcmaxRef=50
Qin=0:2000
param=f.make.param(VcmaxRef=VcmaxRef,TpRef=9999,thetacj = 1,JmaxRef = VcmaxRef*1.67)
physio=f.A(PFD=Qin,Tleaf=27+273.15,Cs=400,Tair=25+273.15,RH=70,
           param=param)
plot(x=Qin,y=physio$A,type='l',col='darkgrey',lwd=4,ylim=c(-1.5,12),xlim=c(0,2000),xlab=expression(italic(Q)~mu*mol~m^-2~s^-1),ylab=expression(italic(A)[n]~mu*mol~m^-2~s^-1))
lines(x=Qin,y=rep(physio$Ac,length(physio$Aj)),type='l',col='darkblue',lwd=2,lty=2)
lines(x=Qin,y=physio$Aj,type='l',col='darkred',lwd=2,lty=2)
legend("bottomright",legend=c(expression(italic(A)[n]),expression(italic(A)[c]),expression(italic(A)[j])),lty=c(1,2,2),
       pch=c(NA,NA,NA),
       col=c("dark grey","dark blue","dark red"),bty="n",lwd=c(4,2,2),
       seg.len=2,cex=1,pt.cex=1)
physio$Qin=Qin
dev.off()


## Simulation of A-Q curves with the FvCB and JB models
## This simulation shows that the JB model predicts a low apparent convexity (Theta close to 0) 
## In these simulations, the photosynthesis rate, A at Q = 1800 is identical for the JB and FvCB models
## as a result of using a theoretical equation to find the correspondence between
## Vqmax and Jmax f.Jmax_to_Vqmax() (Eqn 11).
result=result_JB=Theta=data.frame()

## Environmental variables
PFD = seq(0,2500,1)
Tleaf = Tair = 25+273.15
Cs = 400
Tair = 25+273.15
RH = 70


for(JmaxRef in c(25,100,200,300)){
  for(g1 in 2:8){
   
    ## FvCB model
    source('R models/Farquhar_von_Caemmerer_Berry_photosynthesis_model.R')
    
    param=f.make.param(VcmaxRef=9999,TpRef=9999,JmaxRef = JmaxRef,thetacj = 1,thetaip = 1,g1=g1,RdRef=1)
    physio=f.A(PFD=PFD,Tleaf=Tleaf,Cs=Cs,Tair=Tair,RH=RH,param=param)
    physio$JmaxRef=JmaxRef
    physio$Qin=PFD
    physio$g1=g1
    physio=as.data.frame(physio)
    result=rbind.data.frame(result,physio)
    
    ## JB model
    source('R models/Johnson_Berry_photosynthesis_model.R')
    ## Ci1800 is Ci for Qin == 1800
    Ci1800=physio[physio$Qin==1800,"Ci"] 
    ## Ci0 is Ci for Qin == 0
    Ci0=physio[physio$Qin==0,"Ci"]
    param_JB=f.make.param_JB(VcmaxRef=9999,TpRef=9999,VqmaxRef = 9999,thetacj = 1,thetaip = 1,g1=g1,RdRef=1)
    eta0=f.eta(Ci = Ci0,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )

    ## Parameterization of alpha1 such as the quantum yield of FvCB and JB are standardized
    alpha1=param[['aQY']]/param_JB[['phi1max']]*eta0
    
    ## Calculation of VqmaxRef equivalent to JmaxRef.
    VqmaxRef=f.Jmax_to_Vqmax(Q = 1800,Jmax = JmaxRef,Ci = Ci1800,Beta = 1-alpha1) ## Eqn 9
    param_JB=f.make.param_JB(VcmaxRef=9999,TpRef=9999,VqmaxRef = VqmaxRef,thetacj = 1,thetaip = 1,g1=g1,RdRef=1,Beta=1-alpha1)
    physio_JB=f.A_JB(PFD=PFD,Tleaf=Tleaf,Tair=Tair,RH=RH,Cs=Cs,
                       param=param_JB)
    physio_JB$JmaxRef=JmaxRef
    physio_JB$Qin=PFD
    physio_JB$g1=g1
    physio_JB$Tleaf=Tleaf
    
    ## Calculation of the apparent Theta
    fit_Theta=f.fitting(measures=as.data.frame(physio_JB),Start=list(JmaxRef=70,Theta=0.4),param=f.make.param(RdRef = 1,VcmaxRef = 9999,TpRef=9999,thetacj = 1,thetaip = 1),type = 'AQ',modify.init = FALSE,do.plot = TRUE)
    Theta=rbind.data.frame(Theta,data.frame(Theta=fit_Theta[[2]]@coef[["Theta"]],g1=g1,JmaxRef=JmaxRef))
    result_JB=rbind.data.frame(result_JB,as.data.frame(physio_JB))
    }  
}

hist(Theta$Theta)

## Creating a figure with the AQ curves simulated by both models
result$model='FvCB'
result_JB$model='JB'
result=result[,colnames(result_JB)[colnames(result_JB)%in%colnames(result)]]
result_all=rbind.data.frame(result,result_JB[,colnames(result_JB)[colnames(result_JB)%in%colnames(result)]])

result_all=result_all[order(result_all$Qin),]
result_all$group=paste(result_all$JmaxRef,result_all$model)
a=ggplot(data=result_all[result_all$g1==4,],aes(x=Qin,y=A,color=model,group=group))+
  geom_line()+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank())+
  ylab(expression(italic(A)[n]~mu*mol~m^-2~s^-1)) + ylim(c(-2,50))+
  xlab(expression(italic(Q)~mu*mol~m^-2~s^-1)) +
  annotate('text', x = 1700, y = result_all[result_all$JmaxRef==25&result_all$Q==1700&result_all$model=="FvCB"&result_all$g1==4,"A"]+1.5,label = "italic(J)[max]~'='~25",parse = TRUE,size=3.5) +
  annotate('text', x = 1700, y = result_all[result_all$JmaxRef==100&result_all$Q==1700&result_all$model=="FvCB"&result_all$g1==4,"A"]+1.5,label = "italic(J)[max]~'='~100",parse = TRUE,size=3.5) +
  annotate('text', x = 1700, y = result_all[result_all$JmaxRef==200&result_all$Q==1700&result_all$model=="FvCB"&result_all$g1==4,"A"]+2,label = "italic(J)[max]~'='~200",parse = TRUE,size=3.5) +
  annotate('text', x = 1700, y = result_all[result_all$JmaxRef==300&result_all$Q==1700&result_all$model=="FvCB"&result_all$g1==4,"A"]+2.75,label = "italic(J)[max]~'='~300",parse = TRUE,size=3.5)

print(a)
jpeg(filename = 'Figures/Figure_1.jpeg',width = 95*1.3,height = 68*1.3,units = 'mm',res=300)
a
dev.off()

