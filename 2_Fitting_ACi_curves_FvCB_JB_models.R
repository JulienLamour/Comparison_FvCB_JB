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


Do_ACi_import = FALSE ## Change to TRUE if you want to import the A-Ci curves (Quick)
Do_Fitting = FALSE ## Change to TRUE if you want to fit the curves (the fitting takes around 30 minutes on my computer)


########################################
###   Importing the A-Ci datasets    ###
########################################
## For each dataset, we used a specific R code to import the authors data that we called "0_Import_transform_ACi_data.R"
## we did a QaQc of the data to remove bad curves or bad points in the file "1_QaQc_curated_ACi.R".
## Here we run all these codes to create a combined dataset that we then analyse



if(Do_ACi_import){
  ## 1) Run all the datasets import processes to be sure that everything is up to date
  
  ls_Rfiles = dir(path = path, recursive = TRUE)
  ls_Rfiles = ls_Rfiles[which(grepl(x=ls_Rfiles,pattern="0_Import_transform_ACi_data.R",ignore.case = TRUE))]
  ls_Rfiles = file.path(path,ls_Rfiles)
  for(file in ls_Rfiles){
    print(file)
    source(file)
  }
  
  ## 2) Run all the QAQC to be sure that everything is up to date
  ls_Rfiles = dir(path = path, recursive = TRUE)
  ls_Rfiles = ls_Rfiles[which(grepl(x=ls_Rfiles,pattern="1_QaQc_curated_ACi.R",ignore.case = TRUE))]
  ls_Rfiles = file.path(path,ls_Rfiles)
  for(file in ls_Rfiles){
    print(file)
    source(file)
  }
  
  ## 3) Combining the datasets
  
  ACi_colnames=c("SampleID","Record","A","Ci","CO2s","Patm","Qin","RHs","Tleaf","QC","Dataset")
  ls_Rdata = dir(path = path, recursive = TRUE)
  ls_Rdata = ls_Rdata[which(grepl(x=ls_Rdata,pattern="1_QC_ACi_data.Rdata",ignore.case = TRUE))]
  ls_Rdata = file.path(path,ls_Rdata)
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
  
  save(ACi,data_Qin,file = file.path(path,'/Outputs/Curated_ACi_database.Rdata'))
  
} else (load(file = file.path(path,'/Outputs/Curated_ACi_database.Rdata')))

##############################################################
###   Fitting the ACi curves with the FvCB and JB models   ###
##############################################################
setwd(path)
## Here, we fit the FvCB model to the data using the parameters from Caemmerer et al. 2009. 
## i.e: abso = 0.85,aQY = 0.425, and Theta= 0.7 
## Von Caemmerer, Susanne & Farquhar, Graham & Berry, Joseph. (2009). Biochemical Model of C3 Photosynthesis. 10.1007/978-1-4020-9237-4_9. 

## Creation of a list of parameters for the FvCB model
## Note that I put 0 in the temperature dependence parameters so Vcmax, Jmax etc.. are estimated at Tleaf and not at 25°C
param=f.make.param(RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Creation of a list of parameters for the JB model (this first list is then modified)
param_JB=f.make.param_JB()
## We standardize the quantum yield so it is similar between models
eta800=f.eta(Ci = 800,Gstar = param_JB[["GstarRef"]],etal = param_JB[["etalRef"]],etac = param_JB[["etacRef"]] )
alpha1=param[['aQY']]/param_JB[['phi1max']]*eta800 ## Eqn 10
## Creation of a list of parameters for the JB model
param_JB=f.make.param_JB(Beta=1-alpha1,RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,VqmaxHa = 0,etacHd = 0,etacS = 0,etalHd = 0,etalS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Note that the functions f.fit_ACi creates pdf that show the fittings of each curves, with the Ac, Ac+Aj and Ac+Aj+Ap limitations.
## One pdf called "2_ACi_fitting_best_model.pdf" also showed the best model.
if(Do_Fitting){
  #ACi = ACi[ACi$SampleID_num%in% unique(ACi$SampleID_num)[sample.int(size = 30, n = length(unique(ACi$SampleID_num)))],]
  Bilan_FvCB=f.fit_ACi(measures=ACi,param = param)
  Bilan_JB=f.fit_ACi_JB(measures=ACi,param = param_JB )
  save(Bilan_FvCB,Bilan_JB,data_Qin,file=file.path(path,"Outputs","Fitted_ACi_curves.Rdata"))
} else {load('Outputs/Fitted_ACi_curves.Rdata',verbose=TRUE) }

Bilan=cbind.data.frame(Bilan_FvCB,Bilan_JB) ## (NB: Bilan is a french word which I wrongly assumed also meant something in English)

## Adding the light used for measuring the A-Ci curves
Bilan=merge(x=Bilan,y=data_Qin[,c('SampleID','SampleID_num','Qin','Dataset')],by.x='SampleID_num',by.y='SampleID_num')

## Comparing only the fittings that use the same limitations.
Bilan=Bilan[-which(Bilan$model!=Bilan$JB_model),] 



## Here we fit the FvCB model to the data using Medlyn et al. (2002) set of parameters for defining the quantum yield and Theta
## i.e: abso = 0.85,aQY = 0.372, and Theta= 0.9
## Medlyn BE, Dreyer E, Ellsworth D, Forstreuter M, Harley PC, Kirschbaum MUF, Le Roux X, Montpied P, Strassemeyer J, Walcroft A, et al. 2002. Temperature response of parameters of a biochemically based model of photosynthesis. II. A review of experimental data. Plant, Cell & Environment 25: 1167–1179.
## Creation of a list of parameters for the FvCB model
## Note that I put 0 in the temperature dependence parameters so Vcmax, Jmax etc.. are estimated at Tleaf and not at 25°C
param_Med=f.make.param(aQY = 0.093*4,Theta = 0.9,RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Creation of a list of parameters for the JB model (this first list is then modified)
param_JB_Med=f.make.param_JB()
## We standardize the quantum yield so it is similar between models
eta800_Med=f.eta(Ci = 800,Gstar = param_JB_Med[["GstarRef"]],etal = param_JB_Med[["etalRef"]],etac = param_JB_Med[["etacRef"]] )
alpha1_Med=param_Med[['aQY']]/param_JB_Med[['phi1max']]*eta800_Med

## Creation of a list of parameters for the JB model
param_JB_Med=f.make.param_JB(Beta=1-alpha1_Med,RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,VqmaxHa = 0,etacHd = 0,etacS = 0,etalHd = 0,etalS = 0,TpHa = 0,TpHd = 0,TpS = 0)

## Note that the functions f.fit_ACi creates pdf that show the fittings of each curves, with the Ac, Ac+Aj and Ac+Aj+Ap limitations.
## One pdf called "2_ACi_fitting_best_model.pdf" also showed the best model.
if(Do_Fitting){
  Bilan_FvCB_Med=f.fit_ACi(measures=ACi,param = param_Med)
  Bilan_JB_Med=f.fit_ACi_JB(measures=ACi,param = param_JB_Med )
  save(Bilan_FvCB_Med,Bilan_JB_Med,data_Qin,file=file.path(path,"Outputs","Fitted_ACi_curves_Med.Rdata"))
} else {load('Outputs/Fitted_ACi_curves_Med.Rdata',verbose=TRUE) }

Bilan_Med=cbind.data.frame(Bilan_FvCB_Med,Bilan_JB_Med) ## (NB: Bilan is a french word which I wrongly assumed also meant something in English)

## Adding the light used for measuring the A-Ci curves
Bilan_Med=merge(x=Bilan_Med,y=data_Qin[,c('SampleID','SampleID_num','Qin','Dataset')],by.x='SampleID_num',by.y='SampleID_num')

## Comparing only the fittings that use the same limitations.
Bilan_Med=Bilan_Med[-which(Bilan_Med$model!=Bilan_Med$JB_model),] 


###################################
###   Comparing the parameters  ###
###################################

## Estimating Vqmax from Jmax using Eqn 11 parameterized with a default value for the light irradiance used for measuring the A-Ci curves
## and an average value of Ci
Bilan$Gstar=f.arrhenius(PRef = param[['GstarRef']],Ha = param[['GstarHa']],Tleaf = Bilan$Tleaf)
Bilan$Vqmax_pred_1800=f.Jmax_to_Vqmax(Q = 1800,Gstar = Bilan$Gstar,Ci=800,Jmax=Bilan$JmaxRef,Beta = 1-alpha1)

## Performance of the prediction
plot(Bilan$JB_VqmaxRef,Bilan$Vqmax_pred_1800)
abline(c(0,1))

## Estimating Vqmax from Jmax using Eqn 11 parameterized with the light irradiance used for measuring the A-Ci curves
## and an average value of Ci
Bilan$Vqmax_pred=f.Jmax_to_Vqmax(Q = Bilan$Qin,Gstar = Bilan$Gstar,Ci=800,Beta=1-alpha1,Jmax=Bilan$JmaxRef)
plot(Bilan$JB_VqmaxRef,Bilan$Vqmax_pred)
abline(c(0,1))
Bilan$Q=as.numeric(Bilan$Qin)

## Figures
library(patchwork)
get_fit_stats = function(model) {
  res = residuals(model)
  y_obs = model$model[[1]] # Extract observed response values
  rmse = sqrt(mean(res^2))
  nrmse = (rmse / mean(y_obs)) * 100
  list(rmse = rmse, nrmse = nrmse)
}

theme_shared = theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 10)
  )

q_color_scale = scale_color_viridis_c(
  name = expression(Q~(mu*mol~m^-2~s^-1)),
  option = "D"
)

reg_a = lm(JB_VqmaxRef ~ 0 + I(JmaxRef^2) + JmaxRef, data = Bilan)
stats_a = get_fit_stats(reg_a)

eq_a = sprintf(
  "atop(y == %.3f * x^2 %s %.2f * x, RMSE == %.1f ~ '(%.1f%%)')",
  coef(reg_a)[1],
  ifelse(coef(reg_a)[2] >= 0, "+", "-"),
  abs(coef(reg_a)[2]),
  stats_a$rmse,
  stats_a$nrmse
)

p_a = ggplot(Bilan, aes(x = JmaxRef, y = JB_VqmaxRef, color = Q)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.7, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ 0 + I(x^2) + x, color = "red", se = FALSE,linewidth = 0.7) +
  annotate("text", x = 100, y = 50, label = eq_a, parse = TRUE, hjust = 0, size = 3.2) +
  q_color_scale + coord_cartesian(xlim = c(0, 650), ylim = c(0, 650)) +
  labs(x = expression(Fitted~italic(J)[max]~(mu*mol~m^-2~s^-1)),
    y = expression(Fitted~italic(Cb)[6]*italic(f)[max]~(mu*mol~m^-2~s^-1))) + theme_shared

reg_b = lm(Vqmax_pred ~ JB_VqmaxRef, data = Bilan)
stats_b = get_fit_stats(reg_b)

eq_b = sprintf(
  "atop(y == %.2f * x %s %.2f, RMSE == %.1f ~ '(%.1f%%)')",
  coef(reg_b)[2],
  ifelse(coef(reg_b)[1] >= 0, "+", "-"),
  abs(coef(reg_b)[1]),
  stats_b$rmse,
  stats_b$nrmse
)

p_b = ggplot(Bilan, aes(x = JB_VqmaxRef, y = Vqmax_pred, color = Q)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.7, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE, linewidth = 0.7) +
  annotate("text", x = 20, y = 600, label = eq_b, parse = TRUE, hjust = 0, size = 3.2) +
  q_color_scale + coord_cartesian(xlim = c(0, 650), ylim = c(0, 650)) +
  labs(x = expression(Fitted~italic(Cb)[6]*italic(f)[max]~(mu*mol~m^-2~s^-1)),
    y = expression(Predicted~italic(Cb)[6]*italic(f)[max]~(mu*mol~m^-2~s^-1))) +theme_shared

reg_c = lm(Vqmax_pred_1800 ~ JB_VqmaxRef, data = Bilan)
stats_c = get_fit_stats(reg_c)

eq_c = sprintf(
  "atop(y == %.2f * x %s %.2f, RMSE == %.1f ~ '(%.1f%%)')",
  coef(reg_c)[2],
  ifelse(coef(reg_c)[1] >= 0, "+", "-"),
  abs(coef(reg_c)[1]),
  stats_c$rmse,
  stats_c$nrmse
)

p_c = ggplot(Bilan, aes(x = JB_VqmaxRef, y = Vqmax_pred_1800, color = Q)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 1.7, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE, linewidth = 0.7) +
  annotate("text", x = 20, y = 600, label = eq_c, parse = TRUE, hjust = 0, size = 3.2) +
  q_color_scale + coord_cartesian(xlim = c(0, 650), ylim = c(0, 650)) +
  labs( x = expression(Fitted~italic(Cb)[6]*italic(f)[max]~(mu*mol~m^-2~s^-1)),
    y = expression(Predicted~italic(Cb)[6]*italic(f)[max]~(mu*mol~m^-2~s^-1))) + theme_shared

final_figure = p_a + p_b + p_c +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")

print(final_figure)
setwd(path)
jpeg(filename = 'Figures/Figure_5.jpeg',width = 170,height = 150,units = 'mm',res=300)
print(final_figure)
dev.off()


## Ci transition

f.ci.treshold<-function(PFD,Tleaf,param){
  
  Kc=f.arrhenius(param[['KcRef']],param[['KcHa']],Tleaf)
  Ko=f.arrhenius(param[['KoRef']],param[['KoHa']],Tleaf)
  Gstar=f.arrhenius(param[['GstarRef']],param[['GstarHa']],Tleaf)
  
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],param[['VcmaxHa']],param[['VcmaxHd']],param[['VcmaxS']],Tleaf)
  Jmax=f.modified.arrhenius(PRef=param[['JmaxRef']],param[['JmaxHa']],param[['JmaxHd']],param[['JmaxS']],Tleaf)
  
  I2=PFD*param[['abso']]*(param[['aQY']])
  if(param[['Theta']]==0){J=I2*Jmax/(I2+Jmax)}else{J=(I2+Jmax-((I2+Jmax)^2-4*(param[['Theta']])*I2*Jmax)^0.5)/(2*(param[['Theta']]))}
  
  ci_t=(J*Kc*(1+param[['O2']]/Ko)-8*Gstar*Vcmax)/(4*Vcmax-J)
  return(ci_t)
}

for (i in 1:nrow(Bilan)){
  param=f.make.param(RdHa=0,RdHd = 0,RdS = 0,VcmaxHa = 0,VcmaxHd = 0,VcmaxS = 0,JmaxHa = 0,JmaxHd = 0,JmaxS = 0,TpHa = 0,TpHd = 0,TpS = 0)
  param[["VcmaxRef"]] = Bilan[i,"VcmaxRef"]
  param[["JmaxRef"]] = Bilan[i,"JmaxRef"]
  param[["RdRef"]] = Bilan[i,"RdRef"]
  Bilan[i, "Ci_transition"] = f.ci.treshold(PFD = Bilan[i,"Qin"],Tleaf = Bilan[i,"Tleaf"], param = param)
}

hist(Bilan$Ci_transition, na.rm = TRUE)
