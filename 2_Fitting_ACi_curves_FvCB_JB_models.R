###############################################################################
###  The aim of this code is to fit measured A-Ci curves and compare the    ###
###  parameters estimated with the FvCB and JB models.                      ###
###  We also estimate the accuracy of Eqn 11 to estimate Vqmax from Jmax    ###
###############################################################################
library(ggplot2)
library(cowplot)
library(bbmle)
library(here)

path=here()
setwd(path)
source('R models/Farquhar_von_Caemmerer_Berry_photosynthesis_model.R')
source('R models/Johnson_Berry_photosynthesis_model.R')
source('R models/fit_ACi.R')
source('R models/Tools.R')

########################################
###   Importing the A-Ci datasets    ###
########################################
## For each dataset, we used a specific R code to import the authors data that we called "0_Import_transform_ACi_data.R"
## we did a QaQc of the data to remove bad curves or bad points in the file "1_QaQc_curated_ACi.R".
## Here we run all these codes to create a combined dataset that we then analyse

## 1) Run all the datasets import processes to be sure that everything is up to date

ls_Rfiles <- dir(path = path, recursive = TRUE)
ls_Rfiles <- ls_Rfiles[which(grepl(x=ls_Rfiles,pattern="0_Import_transform_ACi_data.R",ignore.case = TRUE))]
ls_Rfiles <- file.path(path,ls_Rfiles)
for(file in ls_Rfiles){
  print(file)
  source(file)
}

## 2) Run all the QAQC to be sure that everything is up to date

ls_Rfiles <- dir(path = path, recursive = TRUE)
ls_Rfiles <- ls_Rfiles[which(grepl(x=ls_Rfiles,pattern="1_QaQc_curated_ACi.R",ignore.case = TRUE))]
ls_Rfiles <- file.path(path,ls_Rfiles)
for(file in ls_Rfiles){
  print(file)
  source(file)
}

## 3) Combining the datasets

ACi_colnames=c("SampleID","Record","A","Ci","CO2s","Patm","Qin","RHs","Tleaf","QC","Dataset")
ls_Rdata <- dir(path = path, recursive = TRUE)
ls_Rdata <- ls_Rdata[which(grepl(x=ls_Rdata,pattern="1_QC_ACi_data.Rdata",ignore.case = TRUE))]
ls_Rdata <- file.path(path,ls_Rdata)
ACi=data.frame()
for(file in ls_Rdata){
  print(file)
  load(file,verbose=TRUE)
  ACi=rbind.data.frame(ACi,curated_data[,ACi_colnames])
}

ACi$Tleaf=ACi$Tleaf+273.16 ## Conversion of the leaf temperature in kelvin
ACi$SampleID_num=as.numeric(as.factor(paste(ACi$Dataset,ACi$SampleID))) ## Creating a SampleID numerical identifier
ACi=ACi[order(ACi$SampleID_num,ACi$Ci),] ## Sorting the points in the ACi curves so the Ci are in an increasing order. It helps with the plots
data_Qin=ACi[-which(duplicated(ACi$SampleID_num)),] ## Keeping the first point of each curves

##############################################################
###   Fitting the ACi curves with the FvCB and JB models   ###
##############################################################
setwd(path)
Do_Fitting = FALSE ## Change to TRUE if you want to fit the curves (the fitting takes around 30 minutes on my computer)

## Creation of a list of parameters for the FvCB model
## Note that I put 0 in the temperature dependence parameters so Vcmax, Jmax etc.. are estimated at Tleaf and not at 25°C
param=f.make.param(RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Creation of a list of parameters for the JB model (this first list is then modified)
param_JB=f.make.param_JB()
## We standardize the quantum yield so it is similar between models (Eqn 10)
eta800=f.eta(Ci = 800,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )
alpha1=param[['aQY']]/param_JB[['phi1max']]*eta800 ## Eqn 10
## Creation of a list of parameters for the JB model
param_JB=f.make.param_JB(Beta=1-alpha1,RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,VqmaxHa = 0,etacHd = 0,etacS = 0,etalHd = 0,etalS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Note that the functions f.fit_ACi creates pdf that show the fittings of each curves, with the Ac, Ac+Aj and Ac+Aj+Ap limitations.
## One pdf called "2_ACi_fitting_best_model.pdf" also showed the best model.
if(Do_Fitting){
  Bilan_FvCB=f.fit_ACi(measures=ACi,param = param)
  Bilan_JB=f.fit_ACi_JB(measures=ACi,param = param_JB )
  save(Bilan_FvCB,Bilan_JB,data_Qin,file=file.path(path,"Outputs","Fitted_ACi_curves.Rdata"))
} else {load('Outputs/Fitted_ACi_curves.Rdata',verbose=TRUE) }

Bilan=cbind.data.frame(Bilan_FvCB,Bilan_JB)

## Adding the light used for measuring the A-Ci curves
Bilan=merge(x=Bilan,y=data_Qin[,c('SampleID','SampleID_num','Qin','Dataset')],by.x='SampleID_num',by.y='SampleID_num')

## Comparing only the fittings that use the same limitations.
Bilan=Bilan[-which(Bilan$model!=Bilan$JB_model),] 


###################################
###   Comparing the parameters  ###
###################################

## Estimating Vqmax from Jmax using Eqn 11 parameterized with a default value for the light irradiance used for measuring the A-Ci curves
## and an average value of Ci
Bilan$Gstar=f.arrhenius(PRef = param[['GstarRef']],Ha = param[['GstarHa']],Tleaf = Bilan$Tleaf)
Bilan$Vqmax_pred_1300=f.Jmax_to_Vqmax(Q = 1300,Gstar = Bilan$Gstar,Ci=800,Jmax=Bilan$JmaxRef,Beta = 1-alpha1)

## Performance of the prediction
plot(Bilan$JB_VqmaxRef,Bilan$Vqmax_pred_1300)
abline(c(0,1))
reg_1=lm(Bilan$Vqmax_pred_1300~Bilan$JB_VqmaxRef)
summary(reg_1)
R2_1=summary(reg_1)$adj.r.squared
RMSE_1=summary(reg_1)$sigma

## Estimating Vqmax from Jmax using Eqn 11 parameterized with a default value for the light irradiance used for measuring the A-Ci curves
## and an average value of Ci
Bilan$Vqmax_pred=f.Jmax_to_Vqmax(Q = Bilan$Qin,Gstar = Bilan$Gstar,Ci=800,Beta=1-alpha1,Jmax=Bilan$JmaxRef)
plot(Bilan$JB_VqmaxRef,Bilan$Vqmax_pred)
abline(c(0,1))
reg_2=lm(Bilan$Vqmax_pred~Bilan$JB_VqmaxRef)
R2_2=summary(reg_2)$adj.r.squared
RMSE_2=summary(reg_2)$sigma
Bilan$Q=as.numeric(Bilan$Qin)

## Figures

size_p=2
shape_p=1

a=ggplot(data=Bilan,aes(x=JmaxRef,y=JB_VqmaxRef,color=Q)) +
  geom_abline(slope=1,intercept=0,color='grey')+ xlim(c(0,650))+ylim(c(0,650)) +
  geom_point(size=size_p,shape=shape_p)+scale_color_gradientn(colours = rainbow(5))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(Fitted~italic(V)[qmax]~mu*mol~m^-2~s^-1)) +
  xlab(expression(Fitted~italic(J)[max]~mu*mol~m^-2~s^-1))

b=ggplot(data=Bilan,aes(x=JB_VqmaxRef,y=Vqmax_pred,color=Q))+
  geom_abline(slope=1,intercept=0,color='grey')+scale_color_gradientn(colours = rainbow(5))+
  geom_point(size=size_p, shape=shape_p)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(Predicted~italic(V)[qmax]~mu*mol~m^-2~s^-1)) +xlim(c(0,650))+ylim(c(0,650))+
  xlab(expression(Fitted~italic(V)[qmax]~mu*mol~m^-2~s^-1)) +annotate('text',x=0,y=600,label=paste('R² =', round(R2_2,3)),hjust=0)+annotate('text',x=0,y=530,label=paste('RMSE =', round(RMSE_2,1)),hjust=0)+labs(color=expression(Q~mu*mol~m^-2~s^-1))

c=ggplot(data=Bilan,aes(x=JB_VqmaxRef,y=Vqmax_pred_1300,color=Q))+
  geom_abline(slope=1,intercept=0,color='grey')+xlim(c(0,650))+ylim(c(0,650))+
  geom_point(size=size_p, shape=shape_p)+scale_color_gradientn(colours = rainbow(5))+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(Predicted~italic(V)[qmax]~mu*mol~m^-2~s^-1)) +
  xlab(expression(Fitted~italic(V)[qmax]~mu*mol~m^-2~s^-1))+annotate('text',x=0,y=600,label=paste('R² =', round(R2_1,2)),hjust=0)+annotate('text',x=0,y=530,label=paste('RMSE =', round(RMSE_1,2)),hjust=0)


setwd(path)
leg=get_legend(b)
jpeg(filename = 'Figures/Figure_3.jpeg',width = 220,height = 75,units = 'mm',res=300)
plot_grid(a+theme(legend.position = "none"),b+theme(legend.position = "none"),c+theme(legend.position = "none"),leg,align='hv',ncol=4,labels = c('(a)','(b)','(c)'),rel_widths = c(0.29,0.29,0.29,0.13),label_x = -0.05)
dev.off()
