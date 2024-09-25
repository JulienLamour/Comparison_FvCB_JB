#############################################################################
###   2 functions (f.fit_AQ and f.fit_AQ_JB) are coded in this file.      ###
###   They fit the FvCB and the JB to A-Q data using the Ac, or Ac + Aj   ###
###   photosynthetic limitations. The most parsimonious model is chosen   ###
###   based on the AIC criterion.                                         ###
###   3 pdf are produced by the functions. They represent the A-Q data    ###
###   as well as the best photosynthetic model fit. The pdf               ###
###   "4_AQ_fitting_best_model.pdf" shows the best model based on the     ###
###   AIC. These functions return a dataframe with the photosynthetic     ###
###   parameters for each A-Q curve                                       ###
#############################################################################


#' @param measures AQ curve dataframe with at least the columns A, Ci, Qin and Tleaf (in K)
#' @param param See f.make.param for details. Determine the parameters used for the fitting
#' @param VcmaxRef Value of VcmaxRef used to initialize the fitting procedure
#' @param JmaxRef Value of JmaxRef used to initialize the fitting procedure
#' @param RdRef Value of RdRef used to initialize the fitting procedure
#' @param Fit.Theta Should Theta be also fitted? TRUE / FALSE
#' @param Theta Value the curvature factor for the electron transport equation used to initialize the fitting procedure. Do not chose values above 0.5 since the fitting procedure automatically chose start values in the range up to 2*given value.
#' @param criterion Criterion for the comparisons of the models (AIC or BIC)
#'
#' @return Return a dataframe (Bilan) with for each curve in row the estimated parameters as well as the model used (Ac and Ac and Aj)
#' @export
#'
#' @examples
f.fit_AQ<-function(measures,param,VcmaxRef=50, JmaxRef=80, RdRef = 2, Theta = 0.45,Fit.Theta=FALSE,criterion='AIC'){
  param[['TpRef']]=9999
  param[['VcmaxRef']]=9999
  
  ## Fitting of the AQ curve using Aj limitation and a fixed value of Theta
  pdf(file = '4_AQ_fitting_Aj.pdf')
  result_Aj=by(data = measures,INDICES = list(measures$SampleID_num),
               FUN = function(x){f.fitting(measures = x,Start = list(JmaxRef=JmaxRef,RdRef=RdRef),
                                           param=param,id.name = 'SampleID_num',type = 'AQ')})
  dev.off()
  
  
  ## Fitting of the AQ curve using Aj and Ac limitations and a fixed value of Theta
  
  pdf(file = '4_AQ_fitting_Aj_Ac.pdf')
  result_Aj_Ac=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){
                    Start = list(JmaxRef=JmaxRef,RdRef=RdRef)
                    modify.init=TRUE
                    if(!is.null(result_Aj[[as.character(unique(x$SampleID_num))]][[2]])){
                      without_Ac_Start=result_Aj[[as.character(unique(x$SampleID_num))]][[2]]@coef
                      for(value in names(Start)){Start[[value]]=without_Ac_Start[value]}
                      # modify.init=FALSE
                    }
                    Start[['VcmaxRef']]=Start[['JmaxRef']]/1.67
                    
                    f.fitting(measures = x,Start = Start,
                              param=param,id.name = 'SampleID_num',modify.init =modify.init, type='AQ')
                  })
  dev.off()
  
  ## Fitting of the AQ curve using Ac Aj limitation and a free Theta
  if (Fit.Theta){
    pdf(file = '4_AQ_fitting_Aj_Ac_Theta.pdf')
    result_Aj_Ac_Theta=by(data = measures,INDICES = list(measures$SampleID_num),
                          FUN = function(x){
                            Start = list(JmaxRef=JmaxRef,RdRef=RdRef,VcmaxRef=VcmaxRef)
                            modify.init=TRUE
                            if(!is.null(result_Aj_Ac[[as.character(unique(x$SampleID_num))]][[2]])){
                              without_Theta_Start=result_Aj_Ac[[as.character(unique(x$SampleID_num))]][[2]]@coef
                              for(value in names(Start)){Start[[value]]=without_Theta_Start[value]}
                              modify.init=FALSE
                            }
                            Start[['Theta']]=0.45
                            
                            f.fitting(measures = x,Start = Start,
                                      param=param,id.name = 'SampleID_num',modify.init =modify.init, type='AQ')
                          })
    dev.off()
    
    
    ###Fitting of the curve using Aj limitation and a free Theta
    pdf(file = '4_AQ_fitting_Aj_Theta.pdf')
    result_Aj_Theta=by(data = measures,INDICES = list(measures$SampleID_num),
                       FUN = function(x){
                         Start = list(JmaxRef=JmaxRef,RdRef=RdRef)
                         modify.init=TRUE
                         if(!is.null(result_Aj[[as.character(unique(x$SampleID_num))]][[2]])){
                           without_Theta_Start=result_Aj[[as.character(unique(x$SampleID_num))]][[2]]@coef
                           for(value in names(Start)){Start[[value]]=without_Theta_Start[value]}
                           modify.init=FALSE
                         }
                         Start[['Theta']]=0.45
                         
                         f.fitting(measures = x,Start = Start,
                                   param=param,id.name = 'SampleID_num',modify.init =modify.init, type='AQ')
                       })
    dev.off()
    
  }
  
  
  ### Extracting the fitting metrics for each model and each curve
  
  res_nlf_Aj=as.data.frame(t(sapply(result_Aj,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
  res_nlf_Aj$SampleID_num=row.names(res_nlf_Aj)
  
  res_nlf_Aj_Ac=as.data.frame(t(sapply(result_Aj_Ac,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
  res_nlf_Aj_Ac$SampleID_num=row.names(res_nlf_Aj_Ac)
  
  if(Fit.Theta){  
    res_nlf_Aj_Ac_Theta=as.data.frame(t(sapply(result_Aj_Ac_Theta,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
    res_nlf_Aj_Ac_Theta$SampleID_num=row.names(res_nlf_Aj_Ac_Theta)
    
    res_nlf_Aj_Theta=as.data.frame(t(sapply(result_Aj_Theta,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
    res_nlf_Aj_Theta$SampleID_num=row.names(res_nlf_Aj_Theta)
  }
  
  res_nlf_Aj$VcmaxRef=NA
  res_nlf_Aj$Theta=param[['Theta']]
  res_nlf_Aj_Ac$Theta=param[['Theta']]
  
  res_nlf_Aj$model='Aj'
  res_nlf_Aj_Ac$model='Aj_Ac'
  
  
  res_nlf_Aj_Ac=res_nlf_Aj_Ac[,colnames(res_nlf_Aj)]
  
  if(Fit.Theta){
    res_nlf_Aj_Theta$VcmaxRef=NA
    res_nlf_Aj_Ac_Theta$model='Aj_Ac_Theta'
    res_nlf_Aj_Theta$model='Aj_Theta'
    res_nlf_Aj_Theta=res_nlf_Aj_Theta[,colnames(res_nlf_Aj)]
    res_nlf_Aj_Ac_Theta=res_nlf_Aj_Ac_Theta[,colnames(res_nlf_Aj)]
  }
  
  ## Finding the best model (Ac or Ac_Aj or Ac_Aj_Ap according to the BIC criterion)
  
  Bilan=res_nlf_Aj
  if (criterion == 'BIC'){
    Bilan[which(res_nlf_Aj_Ac$BIC<Bilan$BIC),]=res_nlf_Aj_Ac[which(res_nlf_Aj_Ac$BIC<Bilan$BIC),]
    if(Fit.Theta){
      Bilan[which(res_nlf_Aj_Ac_Theta$BIC<Bilan$BIC),]=res_nlf_Aj_Ac_Theta[which(res_nlf_Aj_Ac_Theta$BIC<Bilan$BIC),]
      Bilan[which(res_nlf_Aj_Theta$BIC<Bilan$BIC),]=res_nlf_Aj_Theta[which(res_nlf_Aj_Theta$BIC<Bilan$BIC),]
    }
  } else if (criterion == 'AIC'){
    Bilan[which(res_nlf_Aj_Ac$AIC<Bilan$AIC),]=res_nlf_Aj_Ac[which(res_nlf_Aj_Ac$AIC<Bilan$AIC),]
    if(Fit.Theta){
      Bilan[which(res_nlf_Aj_Ac_Theta$AIC<Bilan$AIC),]=res_nlf_Aj_Ac_Theta[which(res_nlf_Aj_Ac_Theta$AIC<Bilan$AIC),]
      Bilan[which(res_nlf_Aj_Theta$AIC<Bilan$AIC),]=res_nlf_Aj_Theta[which(res_nlf_Aj_Theta$AIC<Bilan$AIC),]
    }
  } else {print('Wrong criterion, chose between BIC and AIC')}
  
  colnames(Bilan)=c("sigma","JmaxRef","RdRef","BIC","AIC","Tleaf","SampleID_num","VcmaxRef","Theta","model") 
  
  ### Creating a pdf with the best model for each curve
  
  pdf(file = '4_AQ_fitting_best_model.pdf')
  by(data = measures,INDICES = list(measures$SampleID_num),
     FUN = function(x){
       param_model=param
       param_leg=Bilan[Bilan$SampleID_num==unique(x$SampleID_num),c('VcmaxRef','JmaxRef','RdRef','Theta','sigma')]
       param_leg[is.na(param_leg)]=9999
       param_model[['VcmaxRef']]=param_leg[['VcmaxRef']];param_model[['JmaxRef']]=param_leg[['JmaxRef']];param_model[['RdRef']]=param_leg[['RdRef']];param_model[['Theta']]=param_leg[['Theta']]
       f.plot(measures=x,list_legend = as.list(param_leg),param = param_model,name =unique(x$SampleID_num),type = 'AQ',xlim=c(0,max(x$Qin)))
     }
  )
  dev.off()
  
  return(Bilan)
}



#' @param measures AQ curve dataframe with at least the columns A, Ci, Qin and Tleaf (in K)
#' @param param See f.make.param_JB for details. Determine the parameters used for the fitting
#' @param VcmaxRef Value of VcmaxRef used to initialize the fitting procedure
#' @param JmaxRef Value of JmaxRef used to initialize the fitting procedure
#' @param RdRef Value of RdRef used to initialize the fitting procedure
#' @param criterion Criterion for comparing the models (BIC or AIC)
#'
#' @return Return a dataframe (Bilan) with for each curve in row the estimated parameters as well as the model used (Ac, Ac and Aj, Ac and Aj and Ap)
#' @export
#'
#' @examples
f.fit_AQ_JB<-function(measures,param,VcmaxRef=50, VqmaxRef=90, RdRef = 2,criterion='AIC'){
  param[['VcmaxRef']]=9999
  
  ## Fitting of the AQ curve using Aj limitation
  pdf(file = '4_AQ_fitting_Aj_JB.pdf')
  result_Aj=by(data = measures,INDICES = list(measures$SampleID_num),
               FUN = function(x){f.fitting_JB(measures = x,Start = list(VqmaxRef=VqmaxRef,RdRef=RdRef),
                                              param=param,id.name = 'SampleID_num',type = 'AQ')})
  dev.off()
  
  
  ## Fitting of the AQ curve using Aj and Ac limitations 
  
  pdf(file = '4_AQ_fitting_Aj_Ac_JB.pdf')
  result_Aj_Ac=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){
                    Start = list(VqmaxRef=VqmaxRef,RdRef=RdRef)
                    modify.init=TRUE
                    if(!is.null(result_Aj[[as.character(unique(x$SampleID_num))]][[2]])){
                      without_Ac_Start=result_Aj[[as.character(unique(x$SampleID_num))]][[2]]@coef
                      for(value in names(Start)){Start[[value]]=without_Ac_Start[value]}
                      #modify.init=FALSE
                    }
                    Start[['VcmaxRef']]=Start[['VqmaxRef']]/1.67
                    
                    f.fitting_JB(measures = x,Start = Start,
                                 param=param,id.name = 'SampleID_num',modify.init =modify.init, type='AQ')
                  })
  dev.off()
  
  
  ## Extracting the fitting metrics for each model and each curve
  
  res_nlf_Aj=as.data.frame(t(sapply(result_Aj,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
  res_nlf_Aj$SampleID_num=row.names(res_nlf_Aj)
  
  res_nlf_Aj_Ac=as.data.frame(t(sapply(result_Aj_Ac,FUN = function(x) c(x[[2]]@coef,BIC=BIC(x[[2]]),AIC=AIC(x[[2]]),Tleaf=x[[3]]['Tleaf']))))
  res_nlf_Aj_Ac$SampleID_num=row.names(res_nlf_Aj_Ac)
  
  res_nlf_Aj$VcmaxRef=NA
  res_nlf_Aj$model='Aj_JB'
  res_nlf_Aj_Ac$model='Aj_Ac_JB'
  
  res_nlf_Aj=res_nlf_Aj[,colnames(res_nlf_Aj_Ac)]
  ## Finding the best model (Ac or Ac_Aj or Ac_Aj_Ap according to the BIC or AIC criterion)
  
  Bilan=res_nlf_Aj
  if (criterion=='BIC'){
    Bilan[which(res_nlf_Aj_Ac$BIC<Bilan$BIC),]=res_nlf_Aj_Ac[which(res_nlf_Aj_Ac$BIC<Bilan$BIC),]
    
  } else if (criterion=='AIC') {
    Bilan[which(res_nlf_Aj_Ac$AIC<Bilan$AIC),]=res_nlf_Aj_Ac[which(res_nlf_Aj_Ac$AIC<Bilan$AIC),]
  } else {print('Wrong criterion, chose between BIC and AIC')}
  
  colnames(Bilan)=c('sigma','VcmaxRef','VqmaxRef','RdRef','BIC','AIC','Tleaf','SampleID_num','model')
  
  ### Creating a pdf with the best model for each curve
  
  pdf(file = '4_AQ_fitting_best_model_JB.pdf')
  by(data = measures,INDICES = list(measures$SampleID_num),
     FUN = function(x){
       param_model=param
       param_leg=Bilan[Bilan$SampleID_num==unique(x$SampleID_num),c('VcmaxRef','VqmaxRef','RdRef','sigma')]
       param_leg[is.na(param_leg)]=9999
       param_model[['VcmaxRef']]=param_leg[['VcmaxRef']];param_model[['VqmaxRef']]=param_leg[['VqmaxRef']];param_model[['RdRef']]=param_leg[['RdRef']]
       f.plot_JB(measures=x,list_legend = as.list(param_leg),param = param_model,name =unique(x$SampleID_num),type = 'AQ')
     }
  )
  dev.off()
  
  return(Bilan)
}



