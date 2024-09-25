###############################################################
###   The aim of this code is to compare and analyse the    ###
###   performance of the JB and FvCB models for fitting     ###
###   AQ curves. It uses the fittings made in the previous  ###  
###   step ("4_Fitting_AQ_curves_FvCB_JB_models.R")         ###
###############################################################

library(here)
library(ggplot2)
library(nlme)
path=here()
setwd(path)

## The folowing function comes from: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
## It can be used to calculate confidence and prediction intervals for mixed models
bolker_ci <- function(model, newdat, pred_int = FALSE, conf_level = 0.95) {
  if(class(model) != "lme") {
    stop("works for lme-models only")
  }
  z <- round(qnorm((1-conf_level)/2, lower.tail = FALSE), 2)
  newdat$predicted <- predict(model, newdat, level = 0)
  Designmat <- model.matrix(formula(model)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(model) %*% t(Designmat))
  newdat$se <- sqrt(predvar)
  newdat$conf.high <- newdat$predicted - z*newdat$se
  newdat$conf.low <- newdat$predicted + z*newdat$se
  if(pred_int == TRUE) {
    newdat$se2 <- sqrt(predvar+model$sigma^2)
    newdat$predint_l <- newdat$predicted - z*newdat$se2
    newdat$predint_h <- newdat$predicted + z*newdat$se2
  }
  newdat$x=newdat[,1]
  newdat
}

##############################################################################
###   Comparing the performance of the FvCB and JB models parameterized    ###
###   with various set of values for the quantum yield and curvature       ###
##############################################################################
compar_all=data.frame()
ls_Rdata=dir(path = "Outputs",pattern = "Fitted_AQ_",recursive = TRUE) # List all the Rdata files in the output folder
for (file in ls_Rdata){
  filenames=load(file.path(path,"Outputs",file)) #Load the Rdata file and save the name of the variables included in the Rdata
  Bilan=eval(as.name(filenames[2])) # Load the "Bilan" or "Bilan_JB" dataframe where the fitting results are stored into the variable Bilan. 
  Param=eval(as.name(filenames[1])) # Load the parameters used to fit the JB and FvCB model to the data into the  variable "Param"
  if(ncol(Bilan)==9){ # Check if the model is JB or FvCB based on the dimension of the Bilan data.frame
    eta0=eval(as.name(filenames[5])) # Load eta0
    alpha1=eval(as.name(filenames[6])) # Load alpha1
    Bilan$model="JB"
    Bilan$aQY=round(eta0^-1*alpha1*Param[["phi1max"]],2) ## Calculate aQY based on Equation 10 in the paper
    Bilan$Theta=NA
    Bilan$JmaxRef=NA} else {
    Bilan$model="FvCB"
    Bilan$aQY=round(Param[['aQY']],2) 
    Bilan$Theta=Param[['Theta']]  
    Bilan$VqmaxRef=NA
    }
  compar_all=rbind.data.frame(compar_all,Bilan[,c("sigma", "SampleID_num", "model", "aQY", "Theta")])
}

compar_all$model=paste(compar_all$model,"aQY",compar_all$aQY,"Theta",compar_all$Theta)

## Comparison of model RMSE
lme(sigma~0+model,random=~1|SampleID_num,data=compar_all)

## Among all the models, we specifically compare "Bilan_FvCB", i.e the FvCB model parameterized as in Caemerrer et al. 2009 (Theta = 0.7, aQY = 0.425)
## and Bilan_JB_Beta, the JB model parameterized such as the quantum yield is 0.425 (similar to Caemerrer et al. 2009)
Bilan2=Bilan_FvCB
Bilan2_JB=Bilan_JB_Beta
colnames(Bilan2_JB)=paste(colnames(Bilan2_JB),'JB',sep='_')

compar=cbind.data.frame(Bilan2,Bilan2_JB)
compar=merge(x=compar,y=Table_sample,by.x="SampleID_num")
compar$Delta_RMSE=compar$sigma-compar$sigma_JB ## Calculating the relative performance of the models

## Mixed model with the PFT as a fixed effect and the species as random effects
reg_performance=lme(Delta_RMSE~0+PFT,~1|Species,data=compar,na.action = "na.omit")
summary(reg_performance)
P_val=summary(reg_performance)$tTable
P_val=data.frame(PFT=unique(Table_sample$PFT),P_val=P_val[,5])

CI=bolker_ci(model=reg_performance,newdat = data.frame(PFT=unique(compar$PFT)),pred_int = TRUE,conf_level = 0.95)
size_p=0.4

jpeg(filename = 'Figures/Figure_2a.jpeg',width = 80,height = 80,units = 'mm',res=600)
ggplot(data=compar,aes(x=PFT,y=Delta_RMSE,color=PFT)) +geom_hline(yintercept = 0,color="grey")+geom_boxplot(data=CI,aes(ymin=predicted-1.96*se2,lower=predicted-1.96*se,middle=predicted,upper=predicted+1.96*se,ymax=predicted+1.96*se2,x=PFT,color=PFT),stat="identity",inherit.aes = FALSE)+ geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.5,size=size_p)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(expression(Delta*RMSE~mu*mol~m^-2~s^-1))+ylim(c(-1.01*max(compar$Delta_RMSE),1.01*max(compar$Delta_RMSE)))
dev.off()

## Mixed model with the Species as a random effect nested within the PFT random effect
reg_performance=lme(Delta_RMSE~1,~1|PFT/Species,data=compar,na.action = "na.omit")
summary_reg=summary(reg_performance)
b=ggplot(data=compar,aes(x=Delta_RMSE))+geom_histogram()+geom_vline(xintercept = 0,color="grey")+geom_vline(xintercept = reg_performance$coefficients$fixed,color="black")+geom_vline(xintercept = reg_performance$coefficients$fixed+c(-1,1)*summary_reg$tTable[2],color="black",linetype=2)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Number of A-Q curves")+xlab(expression(Delta~RMSE~mu*mol~m^-2~s^-1))

jpeg(filename = 'Figures/Figure_2b.jpeg',width = 80,height = 80,units = 'mm',res=600)
b
dev.off()

mean(compar$sigma) ## 
mean(compar$sigma_JB)
length(which(compar$Delta_RMSE>0))/length(compar$Delta_RMSE)

##########################################################################
###   Same code as above for the Medlyn et al. 2002 parameterization   ###
##########################################################################
Bilan2=Bilan_FvCB_Med
Bilan2_JB=Bilan_JB_Beta_Med
colnames(Bilan2_JB)=paste(colnames(Bilan2_JB),'JB',sep='_')

compar=cbind.data.frame(Bilan2,Bilan2_JB)
compar=merge(x=compar,y=Table_sample,by.x="SampleID_num")
compar$Delta_RMSE=compar$sigma-compar$sigma_JB ## Calculating the relative performance of the models

## Mixed model with the PFT as a fixed effect and the species as random effects
reg_performance=lme(Delta_RMSE~0+PFT,~1|Species,data=compar,na.action = "na.omit")
summary(reg_performance)
P_val=summary(reg_performance)$tTable
P_val=data.frame(PFT=unique(Table_sample$PFT),P_val=P_val[,5])

CI=bolker_ci(model=reg_performance,newdat = data.frame(PFT=unique(compar$PFT)),pred_int = TRUE,conf_level = 0.95)
size_p=0.4

jpeg(filename = 'Figures/Extended_data_Figure_2a.jpeg',width = 80,height = 80,units = 'mm',res=600)
ggplot(data=compar,aes(x=PFT,y=Delta_RMSE,color=PFT)) +geom_hline(yintercept = 0,color="grey")+geom_boxplot(data=CI,aes(ymin=predicted-1.96*se2,lower=predicted-1.96*se,middle=predicted,upper=predicted+1.96*se,ymax=predicted+1.96*se2,x=PFT,color=PFT),stat="identity",inherit.aes = FALSE)+ geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.5,size=size_p)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(expression(Delta*RMSE~mu*mol~m^-2~s^-1))+ylim(c(-1,1.01*max(compar$Delta_RMSE)))
dev.off()

## Mixed model with the Species as a random effect nested within the PFT random effect
reg_performance=lme(Delta_RMSE~1,~1|PFT/Species,data=compar,na.action = "na.omit")
summary_reg=summary(reg_performance)
b=ggplot(data=compar,aes(x=Delta_RMSE))+geom_histogram()+geom_vline(xintercept = 0,color="grey")+geom_vline(xintercept = reg_performance$coefficients$fixed,color="black")+geom_vline(xintercept = reg_performance$coefficients$fixed+c(-1,1)*summary_reg$tTable[2],color="black",linetype=2)+theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylab("Number of A-Q curves")+xlab(expression(Delta~RMSE~mu*mol~m^-2~s^-1))

jpeg(filename = 'Figures/Extended_data_Figure_2b.jpeg',width = 80,height = 80,units = 'mm',res=600)
b
dev.off()

mean(compar$sigma)
mean(compar$sigma_JB)
length(which(compar$Delta_RMSE>0))/length(compar$Delta_RMSE)


