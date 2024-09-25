###############################################################################
###  The aim of this code is to fit empirical A-Q curves from various       ###
###  datasets using the most common parameterization used in terrestrial    ###
###  biosphere models. The actual comparison of the fitting is made         ###
###  in another file (5_Comparison_FvCB_JB_error.R)                         ###
###############################################################################

library(cowplot)
library(bbmle)
library(ggplot2)
library(here)
path=here()
setwd(path) 

source('R models/Farquhar_von_Caemmerer_Berry_photosynthesis_model.R')
source('R models/Johnson_Berry_photosynthesis_model.R')
source("R models/fit_AQ.R")


###################################
###   Importing the datasets    ###
###################################

## 1) Run all the datasets import to be sure that everything is up to date
path=here()
ls_Rfiles = dir(path = path, recursive = TRUE)
ls_Rfiles = ls_Rfiles[which(grepl(x = ls_Rfiles,pattern = "0_Import_transform_AQ_data.R",ignore.case = TRUE))]
ls_Rfiles = file.path(path,ls_Rfiles)
for(file in ls_Rfiles){
  print(file)
  source(file)
}

## 2) Combine the datasets
AQ_colnames = c("SampleID","Record","A","Ci","CO2s","Patm","Qin","RHs","Tleaf","gsw","Dataset","Species","PFT","SunShade")
ls_Rdata = dir(path = path, recursive = TRUE)
ls_Rdata = ls_Rdata[which(grepl(x = ls_Rdata,pattern="0_curated_data.Rdata",ignore.case = TRUE)&grepl(x = ls_Rdata,pattern="A-Q",ignore.case = TRUE))]
ls_Rdata = file.path(path,ls_Rdata)
AQ = data.frame()
for(file in ls_Rdata){
  print(file)
  load(file,verbose=TRUE)
  AQ = rbind.data.frame(AQ,curated_data[,AQ_colnames])
}

AQ$Tleaf = AQ$Tleaf+273.16 ## Conversion to kelvin
AQ$SampleID_num = as.numeric(as.factor(paste(AQ$Dataset,AQ$SampleID)))

# Creating a table that lists all the samples per dataset
Table_sample = AQ[-which(duplicated(AQ$SampleID_num)),c("SampleID","SampleID_num","Dataset","PFT","Species")]
AQ = AQ[order(AQ$SampleID_num,AQ$Qin),] ## Sorting the points in the AQ curves so the Qin are in an increasing order. It helps with the plots

#########################
###   Data exclusion  ###
#########################

AQ$QC="ok"
# Here we flag bad data
AQ[which(AQ$RHs>85),"QC"]="bad" # Removing high RH data following Busch et al. (2024) treshold
AQ[which(AQ$Ci<150|AQ$Ci>1000|is.na(AQ$Ci)),"QC"]="bad" # Removing extremely low or high Ci values
AQ[which(AQ$gsw<0.01),"QC"]="bad" # Removing low gsw values that bias Ci estimation
AQ=AQ[AQ$QC=="ok",]
n_AQ=tapply(AQ$A,AQ$SampleID_num,length)
AQ=AQ[!AQ$SampleID_num%in%names(n_AQ[n_AQ<6]),] # Removing curves made of too few points to estimate the photosynthetic parameters
min_Q=tapply(X = AQ$Qin,INDEX = AQ$SampleID_num,FUN = min)
AQ=AQ[!AQ$SampleID_num%in%names(min_Q[min_Q>150]),]# Removing curves where the first light level is too high to evaluate the curvature of the curve
N_curves_dataset=tapply(AQ$SampleID_num,AQ$Dataset,function(x){length(unique(x))}) ## Number of curves in each dataset
N_curves=sum(N_curves_dataset)


#############################################################
###   Fitting the AQ curves with the FvCB and JB models   ###
#############################################################
setwd(path)

## Note that I put 0 in the temperature dependence parameters so Vcmax, Jmax etc.. are estimated at Tleaf and not at 25°C
before=Sys.time()

# Here, we use the parameters from Caemmerer et al. 2009. 
# i.e: abso = 0.85,aQY = 0.425, and Theta= 0.7 
#Von Caemmerer, Susanne & Farquhar, Graham & Berry, Joseph. (2009). Biochemical Model of C3 Photosynthesis. 10.1007/978-1-4020-9237-4_9. 
setwd(file.path(path,'/Outputs/Fits_AQ_FvCB_aQY0425_Theta07'))
param_FvCB = f.make.param(RdHa = 0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)
Bilan_FvCB = f.fit_AQ(measures = AQ, param = param_FvCB)
Theta=param_FvCB[["Theta"]]
abso=param_FvCB[["abso"]]
aQY=param_FvCB[["aQY"]]
save(param_FvCB, Bilan_FvCB,Table_sample,Theta,abso,aQY,file="Fitted_AQ_FvCB_aQY0425_Theta07.Rdata")


## Here we use a lower apparent quantum yield and a higher theta, as in Medlyn et al. (2002)
# Medlyn BE, Dreyer E, Ellsworth D, Forstreuter M, Harley PC, Kirschbaum MUF, Le Roux X, Montpied P, Strassemeyer J, Walcroft A, et al. 2002. Temperature response of parameters of a biochemically based model of photosynthesis. II. A review of experimental data. Plant, Cell & Environment 25: 1167–1179.
setwd(file.path(path,'/Outputs/Fits_AQ_FvCB_aQY0372_Theta09'))
param_FvCB_Med = f.make.param(aQY = 0.093*4,Theta = 0.9,RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)
Bilan_FvCB_Med = f.fit_AQ(measures = AQ, param = param_FvCB_Med)
Theta=param_FvCB_Med[["Theta"]]
abso=param_FvCB_Med[["abso"]]
aQY=param_FvCB_Med[["aQY"]]
save(param_FvCB_Med, Bilan_FvCB_Med,Table_sample,Theta,abso,aQY,file="Fitted_AQ_FvCB_aQY0372_Theta09.Rdata")


## Fitting the data with the JB model with Parameterization of alpha1 such as the quantum yield of JB corresponds to that of FvCB Caemmerer
setwd(file.path(path,'/Outputs/Fits_AQ_JB_aQY0425'))
Ci0=mean(AQ[AQ$Qin<1,"Ci"])
param_JB=f.make.param_JB()
eta0=f.eta(Ci = Ci0,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )

alpha1=param_FvCB[['aQY']]/param_JB[['phi1max']]*eta0
param_JB_Beta = f.make.param_JB(Beta=1-alpha1,RdHa = 0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,VqmaxHa = 0,etacHd = 0,etacS = 0,etalHd = 0,etalS = 0,TpHa = 0,TpHd = 0,TpS = 0)
Bilan_JB_Beta = f.fit_AQ_JB(measures = AQ, param = param_JB_Beta)
save(param_JB_Beta, Bilan_JB_Beta,Table_sample,Ci0,eta0,alpha1,file="Fitted_AQ_JB_aQY0425.Rdata")

## Parameterization of alpha1 such as the quantum yield of JB corresponds to that of FvCB Medlyn
setwd(file.path(path,'/Outputs/Fits_AQ_JB_aQY0372'))
Ci0=mean(AQ[AQ$Qin<1,"Ci"])

param_JB=f.make.param_JB()
eta0=f.eta(Ci = Ci0,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )

alpha1=param_FvCB_Med[['aQY']]/param_JB[['phi1max']]*eta0
param_JB_Beta_Med = f.make.param_JB(Beta=1-alpha1,RdHa = 0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,VqmaxHa = 0,etacHd = 0,etacS = 0,etalHd = 0,etalS = 0,TpHa = 0,TpHd = 0,TpS = 0)
Bilan_JB_Beta_Med = f.fit_AQ_JB(measures = AQ, param = param_JB_Beta_Med)
save(param_JB_Beta_Med, Bilan_JB_Beta_Med,Table_sample,Ci0,eta0,alpha1,file="Fitted_AQ_JB_aQY0372.Rdata")

after=Sys.time()