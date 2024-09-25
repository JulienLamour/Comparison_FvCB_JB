################################################################################
###   The aim of this code is to compare the simulation of A-Ci curves with  ###
###   the FvCB and JB models.                                                ###
###   The part of the ACi curve limited by the Rubisco carboxylation         ###
###   rate (Ac) is identical in both models but the part limited by the      ###
###   electron transport rate (Aj) is expected to change due to the Ci       ###
###   effect on J in the JB model and the absence of Ci effect on J in the   ###
###   FvCB model. The difference is expected to be low as the effect of Ci   ###
###   on eta and therefore on J is low. We test this effect quantitatively   ###
###   in this code                                                           ###
################################################################################

library(cowplot)
library(ggplot2)
library(bbmle)
library(here) 
library(ggrepel)

setwd(here())
source('R models/Farquhar_von_Caemmerer_Berry_photosynthesis_model.R')
source('R models/Johnson_Berry_photosynthesis_model.R')

## Simulation of an A-Ci curve using the FvCB model

Tleaf = 25+273.16
Ci = seq(40,1400,10)
Qsat = 1800
Rd25 = 1.5
Ci_cross = 800 # Ci at which we calculate Jmax and Vqmax such as J_FvCB == J_JB


data_all=data.frame()
ls.Vcmax=c(10,50,120,200) ## We consider a wide range of Vcmax and Jmax values to encompass multiple situations
ls.ratio=c(1.3,1.67,2) ## We consider a wide range of Jmax:Vcmax ratio values to encompass multiple situations
for (Vcmax in ls.Vcmax){
  for (ratio in ls.ratio){
    param=f.make.param(VcmaxRef=Vcmax,TpRef=9999,JmaxRef = Vcmax*ratio,RdRef = Rd25)
    ## Simulation of an A-Ci curve with the FvCB model
    simu_FvCB=f.ACi(PFD=Qsat,Tleaf=Tleaf,Ci=Ci,param=param)
    data_FvCB=data.frame(Tleaf=Tleaf,
                    Ci=Ci,
                    Qin=Qsat,
                    Tair=Tleaf,
                    A=simu_FvCB$A,
                    Ac=simu_FvCB$Ac,
                    Aj=simu_FvCB$Aj,
                    Vcmax=Vcmax,
                    ratio=ratio,
                    model="FvCB")
    ## Two parameters need to be standardized so the FvCB and JB model have the same quantum yield and the same rate of electron transport at saturation irradiance.
    ## The first parameter is alpha1 (see Eqn 10 in the article). 
    ## The second parameter is Vqmax (see Eqn 11)
    ## Parameterization of alpha1 such as the quantum yield of FvCB and JB are standardized
    param_JB=f.make.param_JB()
    eta800=f.eta(Ci = Ci_cross,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )
    alpha1=param[['aQY']]/param_JB[['phi1max']]*eta800 ## Eqn 10
    
    ## We use a Vqmax that corresponds to Jmax in the FvCB model using a standardization equation such as J_FvCB (2000) = J_JB (2000)
    Vqmax_th=f.Jmax_to_Vqmax(Q = Qsat,Jmax = Vcmax*ratio,Ci = Ci_cross,Beta = 1-alpha1) ## Eqn 11 
    data_FvCB$Vqmax=Vqmax_th
    ## Simulation of a similar ACi curve using the JB model parameterized with standardized parameters
    param_JB=f.make.param_JB(VcmaxRef=Vcmax,TpRef=9999,VqmaxRef = Vqmax_th,RdRef = Rd25,Beta = 1-alpha1)
    simu_JB=f.ACi_JB(PFD=Qsat,Tleaf=Tleaf,Ci=Ci,param=param_JB)
    data_JB=data.frame(Tleaf=Tleaf,
                       Ci=Ci,
                       Qin=Qsat,
                       Tair=Tleaf,
                       A=simu_JB$A,
                       Ac=simu_JB$Ac,
                       Aj=simu_JB$Aj,
                       Vcmax=Vcmax,
                       ratio=ratio,
                       model="JB",
                       Vqmax=Vqmax_th)
    data=rbind.data.frame(data_FvCB,data_JB)
    data_all=rbind.data.frame(data_all,data)
  }
}

data_all$group=paste(data_all$Vcmax,data_all$ratio)
data_label=data_all[data_all$Ci==1400&data_all$model=="FvCB",]
data_label$label=paste(expression(italic(V)[cmax]~"="),data_label$Vcmax,expression(","~italic(J)[max]~"="),round(data_label$Vcmax*data_label$ratio),expression(","~italic(V)[qmax]~"="),round(data_label$Vqmax),sep="*")
## Creating a figure with both A-Ci curves
a=ggplot(data=data_all,aes(x=Ci,y=A,color=model,group=interaction(Vcmax,ratio,model)))+
  geom_line(aes(linetype=model),linewidth=0.7)+scale_linetype_manual(values = c(1,5))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank())+
  xlab(expression(italic(C)[i]~µmol~mol^-1)) +geom_text_repel(data=data_label,aes(x=Ci,y=A+2,label=label),parse=T,color="black",hjust=0,vjust=0,direction = "y",min.segment.length = 2,force=0.01,force_pull = 100)+ #geom_rect(xmin=300,xmax=1500,ymin=13,ymax=17,color="grey",alpha=0)+
  ylab(expression(italic(A)~mu*mol~m^-2~s^-1)) +
  xlim(c(0,2500)) #+ ylim(c(10,35))
a

## Creating a figure with the effect of Ci on eta
data_eta=data.frame(Ci=0:1500,eta=f.eta(Ci=0:1500))
b=ggplot(data=data_eta,aes(x=Ci,y=eta))+
  geom_line()+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab(expression(italic(C)[i]~µmol~mol^-1)) +
  ylab(expression(italic(eta))) +
  xlim(c(0,1500)) +
  scale_y_continuous(labels = function(x) format(x,digits=3),limits = c(1,9/8)) 


jpeg(filename = 'Figures/Figure_S1.jpeg',width = 170,height = 170,units = 'mm',res=600)
a+theme(legend.position = c(0.15, 0.875))
dev.off()


jpeg(filename = 'Figures/Figure_S2.jpeg',width = 150,height = 150,units = 'mm',res=600)
b
dev.off()