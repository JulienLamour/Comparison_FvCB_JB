###############################################################
###   The aim of this code is to compare and analyse the    ###
###   performance of the JB and FvCB models for fitting     ###
###   AQ curves. It uses the fittings made in the previous  ###  
###   step ("4_Fitting_AQ_curves_FvCB_JB_models.R")         ###
###############################################################

library(here)
library(ggplot2)
library(cowplot)
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

## The code below loads the outputs of the fittings that were produced in the precedent step (4_Fitting_AQ_curves_FvCB_JB_models.R).
## The outputs are combined in a dataframe called compar_all
ls_Rdata=dir(path = "Outputs",pattern = "Fitted_AQ_",recursive = TRUE) # List all the Rdata files in the output folder
for (file in ls_Rdata){
  filenames=load(file.path(path,"Outputs",file)) #Load the Rdata file and save the name of the variables included in the Rdata
}

## Among all the models, we specifically compare "Bilan_FvCB", i.e the FvCB model parameterized as in Caemerrer et al. 2009 (Theta = 0.7, aQY = 0.425)
## and Bilan_JB_Beta, the JB model parameterized such as the quantum yield is 0.425 (similar to Caemerrer et al. 2009)
Bilan2=Bilan_FvCB[[1]]
Bilan2_JB=Bilan_JB_Beta[[1]]
colnames(Bilan2_JB)=paste(colnames(Bilan2_JB),'JB',sep='_')

compar=cbind.data.frame(Bilan2,Bilan2_JB)
compar=merge(x=compar,y=Table_sample,by.x="SampleID_num")
compar$Delta_RMSE=compar$sigma-compar$sigma_JB ## Calculating the relative performance of the models

reg = lm(VqmaxRef_JB ~ 0+ I(JmaxRef^2)+JmaxRef, data = compar)

equation_label = sprintf("y == %.3f * x^2 %s%.2f * x ~~~ (R^2 == %.2f)", 
                         coef(reg)[1], ifelse(coef(reg)[2] >= 0, "+", "-"), abs(coef(reg)[2]),summary(reg)$r.squared)

reg_plot = ggplot(compar[paste(compar$model,"JB",sep="_") == compar$model_JB,], aes(x = JmaxRef, y = VqmaxRef_JB)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  geom_point(shape = 21, fill = "#0072B2", color = "black", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red",formula = y ~ 0+ I(x^2)+x) +
  annotate("text", x = 10, y = 700, label = equation_label, parse = TRUE, hjust = 0) +
  xlim(c(0, 350)) + ylim(c(0, 700)) +
  labs(
    x = expression(Cb[6]*f[max] ~ (mu * mol ~ m^-2 ~ s^-1)),
    y = expression(J[max] ~ (mu * mol ~ m^-2 ~ s^-1))
  ) +
  theme_bw()


reg_plot

## Mixed model with the PFT as a fixed effect and the species as random effects
reg_performance=lme(Delta_RMSE~0+PFT,~1|Species,data=compar,na.action = "na.omit")
summary(reg_performance)
P_val=summary(reg_performance)$tTable
P_val=data.frame(PFT=unique(Table_sample$PFT),P_val=P_val[,5])
Mean_val=summary(reg_performance)$tTable
Mean_val=data.frame(PFT=unique(Table_sample$PFT),Mean_val=Mean_val[,1])

CI=bolker_ci(model=reg_performance,newdat = data.frame(PFT=unique(compar$PFT)),pred_int = TRUE,conf_level = 0.95)
size_p=0.4

asterisk<-function(P_val){
  if(P_val>=0.05){lab=""}
  if(P_val<0.05){lab="*"}
  if(P_val<0.01){lab="**"}
  if(P_val<0.001){lab="***"}
  return(lab)
}

a=ggplot(data=compar,aes(x=PFT,y=Delta_RMSE,color=PFT)) +geom_hline(yintercept = 0,color="grey")+geom_boxplot(data=CI,aes(ymin=predicted-1.96*se2,lower=predicted-1.96*se,middle=predicted,upper=predicted+1.96*se,ymax=predicted+1.96*se2,x=PFT,color=PFT),stat="identity",inherit.aes = FALSE,width=0.7)+ geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.5,size=size_p)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",plot.caption = element_text(hjust = 0.5, vjust = 1, size = 10), plot.caption.position = "panel" )+
  labs(x="",y=expression(Delta*RMSE~mu*mol~m^-2~s^-1),caption="(RMSE FvCB - RMSE JB)")+ylim(c(-1.01*max(compar$Delta_RMSE),1.01*max(compar$Delta_RMSE)))+coord_flip() + 
  annotate("text",label=asterisk(P_val[1,2]),color="black",x = 1+0.4,y=Mean_val[1,2])+ 
  annotate("text",label=asterisk(P_val[2,2]),color="black",x = 2+0.4,y=Mean_val[2,2])+ 
  annotate("text",label=asterisk(P_val[3,2]),color="black",x = 3+0.4,y=Mean_val[3,2])+ 
  annotate("text",label=asterisk(P_val[4,2]),color="black",x = 4+0.4,y=Mean_val[4,2])+ 
  annotate("text",label=asterisk(P_val[5,2]),color="black",x = 5+0.4,y=Mean_val[5,2])+ 
  annotate("text",label=asterisk(P_val[6,2]),color="black",x = 6+0.4,y=Mean_val[6,2])+
  annotate("text",label=asterisk(P_val[7,2]),color="black",x = 7+0.4,y=Mean_val[7,2])

## Mixed model with the Species as a random effect nested within the PFT random effect
reg_performance=lme(Delta_RMSE~1,~1|Species,data=compar,na.action = "na.omit")
summary_reg=summary(reg_performance)

b=ggplot(data=compar,aes(x= Delta_RMSE, fill = cut(Delta_RMSE, breaks = c(-Inf, 0, Inf))))+
  geom_histogram(boundary = 0)+ scale_fill_manual(
    values = c("(-Inf,0]" = "#007680", "(0, Inf]" = "#d57800")) +
  geom_vline(xintercept = 0,color="grey")+geom_vline(xintercept = reg_performance$coefficients$fixed,color="black")+geom_vline(xintercept = reg_performance$coefficients$fixed+c(-1,1)*summary_reg$tTable[2],color="black",linetype=2)+
  labs(x=expression(Delta*RMSE~mu*mol~m^-2~s^-1),caption="(RMSE FvCB - RMSE JB)",y="Number of A-Q curves") +
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.caption = element_text(hjust = 0.5, vjust = 1, size = 10), plot.caption.position = "panel",legend.position = "none")+
  annotate("text",label="17%",color="white",x = -0.25,y=10)+
  annotate("text",label="83%",color="white",x = 0.45,y=35)+xlim(c(-1.05,1.05))+
  annotate("text",label="FvCB is better",color="#007680",x=-1,y=150,hjust=0)+
  annotate("text",label="JB is better",color="#d57800",x=-1,y=125,hjust=0)
print(b)


jpeg(filename = 'Figures/Figure_3.jpeg',width = 160,height = 80,units = 'mm',res=600)
plot_grid(b,a,labels = c("A","B"))
dev.off()

mean(compar$sigma) ## 
mean(compar$sigma_JB)
length(which(compar$Delta_RMSE>0))/length(compar$Delta_RMSE)
table(compar$model)/sum(table(compar$model))
table(compar$model_JB)/sum(table(compar$model_JB))
  
##########################################################################
###   Same code as above for the Medlyn et al. 2002 parameterization   ###
##########################################################################
Bilan2=Bilan_FvCB_Med[[1]]
Bilan2_JB=Bilan_JB_Beta_Med[[1]]
colnames(Bilan2_JB)=paste(colnames(Bilan2_JB),'JB',sep='_')

compar=cbind.data.frame(Bilan2,Bilan2_JB)
compar=merge(x=compar,y=Table_sample,by.x="SampleID_num")
compar$Delta_RMSE=compar$sigma-compar$sigma_JB ## Calculating the relative performance of the models

## Mixed model with the PFT as a fixed effect and the species as random effects
reg_performance=lme(Delta_RMSE~0+PFT,~1|Species,data=compar,na.action = "na.omit")
summary(reg_performance)
P_val=summary(reg_performance)$tTable
P_val=data.frame(PFT=unique(Table_sample$PFT),P_val=P_val[,5])
Mean_val=summary(reg_performance)$tTable
Mean_val=data.frame(PFT=unique(Table_sample$PFT),Mean_val=Mean_val[,1])

CI=bolker_ci(model=reg_performance,newdat = data.frame(PFT=unique(compar$PFT)),pred_int = TRUE,conf_level = 0.95)
size_p=0.4

a=ggplot(data=compar,aes(x=PFT,y=Delta_RMSE,color=PFT)) +geom_hline(yintercept = 0,color="grey")+geom_boxplot(data=CI,aes(ymin=predicted-1.96*se2,lower=predicted-1.96*se,middle=predicted,upper=predicted+1.96*se,ymax=predicted+1.96*se2,x=PFT,color=PFT),stat="identity",inherit.aes = FALSE,width=0.7)+ geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.5,size=size_p)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",plot.caption = element_text(hjust = 0.5, vjust = 1, size = 10), plot.caption.position = "panel" )+
  labs(x="",y=expression(Delta*RMSE~mu*mol~m^-2~s^-1),caption="(RMSE FvCB - RMSE JB)")+ylim(c(-1.01*max(compar$Delta_RMSE),1.01*max(compar$Delta_RMSE)))+coord_flip() + 
  annotate("text",label=asterisk(P_val[1,2]),color="black",x = 1+0.4,y=Mean_val[1,2])+ 
  annotate("text",label=asterisk(P_val[2,2]),color="black",x = 2+0.4,y=Mean_val[2,2])+ 
  annotate("text",label=asterisk(P_val[3,2]),color="black",x = 3+0.4,y=Mean_val[3,2])+ 
  annotate("text",label=asterisk(P_val[4,2]),color="black",x = 4+0.4,y=Mean_val[4,2])+ 
  annotate("text",label=asterisk(P_val[5,2]),color="black",x = 5+0.4,y=Mean_val[5,2])+ 
  annotate("text",label=asterisk(P_val[6,2]),color="black",x = 6+0.4,y=Mean_val[6,2])+
  annotate("text",label=asterisk(P_val[7,2]),color="black",x = 7+0.4,y=Mean_val[7,2])

## Mixed model with the Species as a random effect nested within the PFT random effect
reg_performance=lme(Delta_RMSE~1,~1|Species,data=compar,na.action = "na.omit")
summary_reg=summary(reg_performance)

b=ggplot(data=compar,aes(x= Delta_RMSE, fill = cut(Delta_RMSE, breaks = c(-Inf, 0, Inf))))+
  geom_histogram(boundary=0)+ scale_fill_manual(
    values = c("(-Inf,0]" = "#007680", "(0, Inf]" = "#d57800")) +
  geom_vline(xintercept = 0,color="grey")+geom_vline(xintercept = reg_performance$coefficients$fixed,color="black")+geom_vline(xintercept = reg_performance$coefficients$fixed+c(-1,1)*summary_reg$tTable[2],color="black",linetype=2)+
  labs(x=expression(Delta*RMSE~mu*mol~m^-2~s^-1),caption="(RMSE FvCB - RMSE JB)",y="Number of A-Q curves") +
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.caption = element_text(hjust = 0.5, vjust = 1, size = 10), plot.caption.position = "panel",legend.position = "none")+
  annotate("text",label="11%",color="white",x = -0.25,y=5)+
  annotate("text",label="89%",color="white",x = 0.75,y=50)+
  annotate("text",label="FvCB is better",color="#007680",x=-1,y=150,hjust=0)+
  annotate("text",label="JB is better",color="#d57800",x=-1,y=125,hjust=0)
print(b)

mean(compar$sigma)
mean(compar$sigma_JB)
length(which(compar$Delta_RMSE>0))/length(compar$Delta_RMSE)

jpeg(filename = 'Figures/Figure_S2.jpeg',width = 160,height = 80,units = 'mm',res=600)
plot_grid(b,a,labels = c("A","B"))
dev.off()

### What quantum yield should we use?
#Bilan042=Bilan_JB_Beta[[1]]
#Bilan037=Bilan_JB_Beta_Med[[1]]
#colnames(Bilan037)=paste(colnames(Bilan037),'037',sep='_')
#
#compar=cbind.data.frame(Bilan042,Bilan037)
#compar$Delta_RMSE = compar$sigma-compar$sigma_037
#compar=merge(x=compar,y=Table_sample,by.x="SampleID_num")
#reg_performance=lme(Delta_RMSE~1,~1|Species,data=compar,na.action = "na.omit")
#summary_reg=summary(reg_performance)

#######################################################################
## Comparison of the FvCB model set with Theta = 0 and the JB model ###
#######################################################################
## We have shown that the JB model and the FvCB model parameterized with Theta = 0 are mathematically similar.
## They are not exactly identical mathematicaly because of eta which depends on Ci
## Here, we test wether the FvCB paramaterized with Theta = 0 result in similar fits to the JB model (They do fit curves similarly).
## Comparison of the FvCB model parameterized with theta = 0 and the JB model, for aQY = 0.425
Bilan_JB_Beta2 = Bilan_JB_Beta[[1]]
colnames(Bilan_JB_Beta2)=paste(colnames(Bilan_JB_Beta2),'JB',sep='_')
compar=cbind.data.frame(Bilan_FvCB_0[[1]],Bilan_JB_Beta2)
reg = lm(sigma_JB~sigma,compar)
coef(reg)
summary(reg)

a=ggplot(data=compar,aes(x=sigma,y=sigma_JB)) +
  geom_abline(slope=1,intercept=0,color='grey')+ 
  geom_point(size = 0.5)+xlim(c(0,3.1))+ylim(c(0,3.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(RMSE~JB~mu*mol~m^-2~s^-1)) +
  xlab(expression(RMSE~FvCB~(theta~'='~0)~mu*mol~m^-2~s^-1)) +
  annotate("text", x = 0, y = 3, 
           hjust = 0, vjust = 1,
           label = paste0(
             "Slope = 1.0", "\n",
             "Intercept = 0.0", "\n",
             "R² > 0.99"))
a

## Comparison of the FvCB model parameterized with theta = 0 and the JB model, for aQY = 0.372
colnames(Bilan_JB_Beta_Med[[1]])=paste(colnames(Bilan_JB_Beta_Med[[1]]),'JB',sep='_')
compar=cbind.data.frame(Bilan_FvCB_Med_0[[1]],Bilan_JB_Beta_Med[[1]])
reg = lm(sigma_JB~sigma,compar)
summary(reg)

b=ggplot(data=compar,aes(x=sigma,y=sigma_JB)) +
  geom_abline(slope=1,intercept=0,color='grey')+ 
  geom_point(size=0.5)+xlim(c(0,3.1))+ylim(c(0,3.1))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab(expression(RMSE~JB~mu*mol~m^-2~s^-1)) +
  xlab(expression(RMSE~FvCB~(theta~'='~0)~mu*mol~m^-2~s^-1)) +
  annotate("text", x = 0, y = 3, 
           hjust = 0, vjust = 1,
           label = paste0(
             "Slope = 1.0", "\n",
             "Intercept = 0.0", "\n",
             "R² > 0.99"))
b

jpeg(filename = 'Figures/Figure_S3.jpeg',width = 160,height = 80,units = 'mm',res=600)
plot_grid(a,b,labels = c("A","B"))
dev.off()



###################################################################################
###   Comparison of the fittings with that of the FvCB model with a free Theta  ###
###################################################################################

## Verification : sigma of the model with a free theta should be >= than sigma of the model with a fixed theta
xlim_sigma = c(min(Bilan_FvCB_Free[[1]]$sigma,Bilan_FvCB[[1]]$sigma),max(Bilan_FvCB_Free[[1]]$sigma,Bilan_FvCB[[1]]$sigma))
plot(Bilan_FvCB_Free[[1]]$sigma, Bilan_FvCB[[1]]$sigma, xlim=xlim_sigma,ylim=xlim_sigma)
abline(c(0,1),col="red") ## This is correct

## Comparing AIC
xlim_AIC = c(min(Bilan_FvCB_Free[[1]]$AIC,Bilan_FvCB[[1]]$AIC),max(Bilan_FvCB_Free[[1]]$AIC,Bilan_FvCB[[1]]$AIC))
plot(Bilan_FvCB_Free[[1]]$AIC, Bilan_FvCB[[1]]$AIC, xlim=xlim_AIC,ylim=xlim_AIC)
abline(c(0,1),col="red") 
hist(Bilan_FvCB_Free[[1]]$AIC- Bilan_FvCB[[1]]$AIC)
compar_AIC = Bilan_FvCB_Free[[1]]
compar_AIC$AIC_fixed = Bilan_FvCB[[1]]$AIC
compar_AIC$AIC_fixed_0 = Bilan_FvCB_0[[1]]$AIC
compar_AIC = merge(x = compar_AIC, y = Table_sample, by= "SampleID_num")
res_JB =Bilan_JB_Beta[[1]][,c("SampleID_num","AIC")]
colnames(res_JB) = c("SampleID_num","AIC_JB")
compar_AIC = merge(x = compar_AIC, y = res_JB, by = "SampleID_num")
compar_AIC$Delta_AIC_plus_2=compar_AIC$AIC+2-compar_AIC$AIC_fixed ## AIC Theta Free + 2 - AIC Theta = 0.7
compar_AIC$Delta_AIC_plus_2_0=compar_AIC$AIC+2-compar_AIC$AIC_fixed_0 ## AIC Theta Free + 2 - AIC Theta = 0
compar_AIC$Delta_AIC_plus_2_JB=compar_AIC$AIC+2-compar_AIC$AIC_JB ## AIC Theta Free + 2 - AIC JB = 0

summary(lme(Delta_AIC_plus_2~1,~1|Species,data=compar_AIC,na.action = "na.omit"))
summary(lme(Delta_AIC_plus_2_0~1,~1|Species,data=compar_AIC,na.action = "na.omit"))
summary(lme(Delta_AIC_plus_2_JB~1,~1|Species,data=compar_AIC,na.action = "na.omit"))

length(which(Bilan_FvCB_Free[[1]]$AIC<(Bilan_FvCB[[1]]$AIC-2)))/length(Bilan_FvCB_Free[[1]]$AIC)
length(which(Bilan_FvCB_Free[[1]]$AIC<(Bilan_FvCB_0[[1]]$AIC-2)))/length(Bilan_FvCB_Free[[1]]$AIC)
length(which(Bilan_FvCB_Free[[1]]$AIC<(Bilan_JB_Beta[[1]]$AIC-2)))/length(Bilan_FvCB_Free[[1]]$AIC)

## Estimating Jmax is more problematic when both when Theta is free and Vcmax is estimated togather with jmax
# Define a helper function to safely calculate correlation from a single sample
calc_cor <- function(sample_item, p1 = "Theta", p2 = "JmaxRef") {
  
  # 1. Extract the mle2 S4 object (2nd element in each sample list)
  mle_obj <- sample_item[[2]]
  
  # 2. Extract the vcov matrix
  v_mat <- tryCatch(mle_obj@vcov, error = function(e) NULL)
  
  # 3. Check if vcov is NULL, not a matrix, or contains NAs
  if (is.null(v_mat) || !is.matrix(v_mat) || any(is.na(v_mat)) || any(is.nan(v_mat))) {
    return(NA_real_)
  }
  
  # 4. Check if both parameters exist in the matrix dimnames
  if (!all(c(p1, p2) %in% colnames(v_mat))) {
    return(NA_real_)
  }
  
  # 5. Extract variances and covariance
  var_p1  <- v_mat[p1, p1]
  var_p2  <- v_mat[p2, p2]
  cov_p12 <- v_mat[p1, p2]
  
  # Check for positive non-zero variances
  if (var_p1 <= 0 || var_p2 <= 0) {
    return(NA_real_)
  }
  
  # 6. Compute correlation: Cov(X,Y) / sqrt( Var(X) * Var(Y) )
  cor_val <- cov_p12 / sqrt(var_p1 * var_p2)
  return(cor_val)
}

cor_JmaxRef_VcmaxRef = sapply(Bilan_FvCB_Free[[3]],FUN = function(x) calc_cor(x, p1 = "VcmaxRef", p2 = "JmaxRef"))
cor_JmaxRef_Theta = sapply(Bilan_FvCB_Free[[3]],FUN = function(x) calc_cor(x,p1 = "Theta", p2 = "JmaxRef"))
cor_Theta_VcmaxRef = sapply(Bilan_FvCB_Free[[3]],FUN = function(x) calc_cor(x,p1 = "Theta", p2 = "VcmaxRef"))
n_bad = length(which(is.na(cor_JmaxRef_Theta)|abs(cor_JmaxRef_Theta)>=0.95|is.na(cor_JmaxRef_VcmaxRef)|abs(cor_JmaxRef_VcmaxRef)>=0.95|is.na(cor_Theta_VcmaxRef)|abs(cor_Theta_VcmaxRef)>=0.95))

n_bad = length(which(is.na(cor_JmaxRef_Theta)|abs(cor_JmaxRef_Theta)>=0.95))
n_bad / length(cor_JmaxRef_Theta) ## In 26 % of the cases, estimating JmaxRef is unreliable (high correlation). (Aj and Ac are considered)

n_bad = length(which(is.na(Bilan_FvCB_Free[[1]]$cor_JmaxRef_Theta)|abs(Bilan_FvCB_Free[[1]]$cor_JmaxRef_Theta)>=0.95))
n_bad / length(Bilan_FvCB_Free[[1]]$cor_JmaxRef_Theta) ## In 22% of the cases, estimating JmaxRef is unreliable when considering the best fit (i.e including Aj or Aj + Ac)

#############################
#### Statistics on Theta  ###
#############################

data_Theta = Bilan_FvCB_Free[[1]]
data_Theta[is.na(data_Theta$cor_JmaxRef_Theta)|abs(data_Theta$cor_JmaxRef_Theta)>=0.95,c("Theta","JmaxRef")] = NA
data_Theta$Trans_Theta = log(1 - data_Theta$Theta)
hist(data_Theta$Trans_Theta)
hist(data_Theta$Theta)
data_Theta = merge(x=data_Theta,y=Table_sample,by.x="SampleID_num", by.y = "SampleID_num")
data_Theta = data_Theta[order(as.numeric(data_Theta$SampleID_num)),]
nrow(data_Theta[!is.na(data_Theta$Theta),])

# Mixed model with the PFT as a fixed effect and the species as random effects
reg_Trans_Theta=lme(Trans_Theta~0+PFT,~1|Species,data=data_Theta,na.action = "na.omit")
summary(reg_Trans_Theta)

P_val=summary(reg_Trans_Theta)$tTable
P_val=data.frame(PFT=unique(Table_sample$PFT),P_val=P_val[,5])
Mean_val=summary(reg_Trans_Theta)$tTable
Mean_val=data.frame(PFT=unique(Table_sample$PFT),Mean_val=Mean_val[,1])

CI=bolker_ci(model=reg_Trans_Theta,newdat = data.frame(PFT=unique(data_Theta$PFT)),pred_int = TRUE,conf_level = 0.95)
CI$predicted_untrans = 1 - exp(CI$predicted)
CI$low_predint_untrans = 1  - exp(CI$predint_l)
CI$high_predint_untrans = 1 - exp(CI$predint_h)
CI$low_confint_untrans = 1  - exp(CI$conf.low)
CI$high_confint_untrans = 1 - exp(CI$conf.high)

size_p = 0.4
a=ggplot(data=data_Theta,aes(x=PFT,y=Theta,color=PFT)) +geom_hline(yintercept = 0,color="grey")+
  scale_y_continuous(breaks = sort(unique(c(pretty(data_Theta$Theta), 1))))+
  geom_boxplot(data=CI,aes(ymin=low_predint_untrans,lower=low_confint_untrans,middle=predicted_untrans,upper=high_confint_untrans,ymax=high_predint_untrans,x=PFT,color=PFT),stat="identity",inherit.aes = FALSE,width=0.7)+
  geom_jitter(shape=16, position=position_jitter(0.2),alpha=0.5,size=size_p)+
  theme_bw()+theme(legend.position = "none")+
  labs(x="",y=expression(theta))+coord_flip() +
  annotate("text",label=asterisk(P_val[1,2]),color="black",x = 1+0.4,y=CI[CI$x == "Crop","predicted_untrans"])+ 
  annotate("text",label=asterisk(P_val[2,2]),color="black",x = 2+0.4,y=CI[CI$x == "Gymnosperm","predicted_untrans"])+ 
  annotate("text",label=asterisk(P_val[3,2]),color="black",x = 3+0.4,y=CI[CI$x == "Mediterranean","predicted_untrans"])+ 
  annotate("text",label=asterisk(P_val[4,2]),color="black",x = 4+0.4,y=CI[CI$x == "Subtropical","predicted_untrans"])+
  annotate("text",label=asterisk(P_val[5,2]),color="black",x = 5+0.4,y=CI[CI$x == "Temperate","predicted_untrans"])+ 
  annotate("text",label=asterisk(P_val[6,2]),color="black",x = 6+0.4,y=CI[CI$x == "Tropical","predicted_untrans"])+ 
  annotate("text",label=asterisk(P_val[7,2]),color="black",x = 7+0.4,y=CI[CI$x == "Tundra","predicted_untrans"])
a

reg_performance=lme(Trans_Theta~1,~1|Species,data=data_Theta,na.action = "na.omit")
summary_reg=summary(reg_performance)
CI=bolker_ci(model=reg_performance,newdat = data.frame(PFT=1),pred_int = TRUE,conf_level = 0.95)
CI$predicted_untrans = 1 - exp(CI$predicted)
CI$low_predint_untrans = 1  - exp(CI$predint_l)
CI$high_predint_untrans = 1 - exp(CI$predint_h)
CI$low_confint_untrans = 1  - exp(CI$conf.low)
CI$high_confint_untrans = 1 - exp(CI$conf.high)

b=ggplot(data=data_Theta,aes(x= Theta))+
  geom_histogram(fill = "#3B6D8C", color = "white", linewidth = 0.1, breaks = seq(-12, 1, by = 0.2)) +
  scale_x_continuous(breaks = sort(unique(c(pretty(data_Theta$Theta), 1)))) +
  geom_vline(xintercept = 0,color="grey")+geom_vline(xintercept = CI$predicted_untrans,color="black")+
  geom_vline(xintercept = as.vector(CI[1,c("low_confint_untrans")]),color="black",linetype=2)+
  geom_vline(xintercept = as.vector(CI[1,c("high_confint_untrans")]),color="black",linetype=2)+
  labs(x=expression(theta),y="Number of A-Q curves") +
  theme_bw() 

print(b)

jpeg(filename = 'Figures/Figure_4.jpeg',width = 160,height = 80,units = 'mm',res=600)
plot_grid(b,a,labels = c("A","B"),rel_widths = c(0.4,0.6))
dev.off()


## Effect of Theta free on Jmax

library(ggplot2)

compar = data.frame(Jmax_Free = data_Theta$JmaxRef, Jmax_07 = Bilan_FvCB[[1]]$JmaxRef)
reg = lm(Jmax_07 ~ Jmax_Free, data = compar)

equation_label = sprintf("y == %.2f * x %s %.2f ~~~ (R^2 == %.2f)", 
              coef(reg)[2], ifelse(coef(reg)[1] >= 0, "+", "-"), abs(coef(reg)[1]), summary(reg)$r.squared)

reg_plot = ggplot(compar, aes(x = Jmax_Free, y = Jmax_07)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  geom_point(shape = 21, fill = "#0072B2", color = "black", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = 10, y = 330, label = equation_label, parse = TRUE, hjust = 0) +
  coord_fixed(xlim = c(0, 350), ylim = c(0, 350)) +
  labs(
    x = expression(J[max] ~ (mu * mol ~ m^-2 ~ s^-1) * "," ~ theta ~ Free),
    y = expression(J[max] ~ (mu * mol ~ m^-2 ~ s^-1) * "," ~ theta == 0.7)
  ) +
  theme_bw()
jpeg(filename = 'Figures/Figure_S4.jpeg',width = 80,height = 80,units = 'mm',res=600)
print(reg_plot)
dev.off()
nrow(data_Theta[!is.na(data_Theta$JmaxRef),])
