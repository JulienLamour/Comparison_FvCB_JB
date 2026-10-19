#############################################################################
###   2 functions (f.fit_ACi and f.fit_ACi_JB) are coded in this file.    ###
###   They fit the FvCB and the JB to A-Ci data using the Ac, Ac + Aj,    ###
###   or Ac + Aj + Ap photosynthetic limitations. The most parsimonious   ###
###   model is chosen based on the AIC criterion.                         ###
###   4 pdf are produced by the functions. They represent the A-Ci data   ###
###   as well as the best photosynthetic model fit. The pdf               ###
###   "2_ACi_fitting_best_model.pdf" shows the best model based on the    ###
###   AIC. These functions return a dataframe with the photosynthetic     ###
###   parameters for each A-Ci curve                                      ###
#############################################################################


#' @param measures ACi curve dataframe with at least the columns A, Ci, Qin and Tleaf (in K)
#' @param param See f.make.param for details. Determine the parameters used for the fitting
#' @param VcmaxRef Value of VcmaxRef used to initialize the fitting procedure
#' @param JmaxRef Value of JmaxRef used to initialize the fitting procedure
#' @param RdRef Value of RdRef used to initialize the fitting procedure
#' @param TpRef Value of TpRef used to initialize the fitting procedure
#'
#' @return Return a dataframe (Bilan) with for each curve in row the estimated parameters as well as the model used (Ac, Ac and Aj, Ac and Aj and Ap)
#' @export
#'
#' @examples
f.fit_ACi<-function(measures,param,VcmaxRef=50, JmaxRef=100, RdRef = 2, TpRef= 5){
  ## Basic checks
  if(any(measures$Ci < 0)) {stop("Ci contains negative values")}
  if(any(is.na(c(measures$A,measures$Ci,measures$Tleaf,measures$Qin,measures$SampleID_num)))) {stop("Ci, A, Tleaf or Qin have NA values")}
  if(any(measures$Tleaf < 273)) {stop("Check Tleaf values. Should be in Kelvin")}
  if(is.null(measures$SampleID_num)|is.null(measures$Ci)|is.null(measures$Tleaf)|is.null(measures$Qin)|is.null(measures$A)){stop("Your dataframe measures does not contain the colums SampleID_num or A or Tleaf or Ci")}
  n_points_curves=tapply(X=measures$SampleID_num,INDEX = measures$SampleID_num,FUN = function(x) length(x))
  min_Ci_curve=tapply(X=measures$Ci,INDEX = measures$SampleID_num,FUN = min)
  if(any(n_points_curves<3)){stop("Some of your A-Ci curves have less than 3 points")}
  if(any(min_Ci_curve>270)){stop("Some of your A-Ci curves have a minimum Ci above 270 ppm meaning that Vcmax cant be estimated")}
  
  param[['TpRef']]=9999
  
  ## Fitting of the ACi curve using Ac and Aj limitation
  pdf(file = 'Outputs/2_ACi_fitting_Ac_Aj.pdf')
  result_Ac_Aj=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){print(paste('SampleID_num:',unique(x$SampleID_num)))
                    f.fitting(measures = x,Start = list(JmaxRef=JmaxRef,RdRef=RdRef,VcmaxRef=VcmaxRef),
                              param=param,id.name = 'SampleID_num')})
  dev.off()
  
  param[['JmaxRef']]=9999
  ## Fitting of the ACi curve using only the Ac limitation
  pdf(file = 'Outputs/2_ACi_fitting_Ac.pdf')
  result_Ac=by(data = measures,INDICES = list(measures$SampleID_num),
               FUN = function(x){print(paste('SampleID_num:',unique(x$SampleID_num)))
                 f.fitting(measures = x,Start = list(RdRef=RdRef,VcmaxRef=VcmaxRef),
                           param=param,id.name = 'SampleID_num')})
  dev.off()
  
  ##Fitting of the curve using the Ac, Aj and Ap limitations
  pdf(file = 'Outputs/2_ACi_fitting_Ac_Aj_Ap.pdf')
  result_Ac_Aj_Ap=by(data = measures,INDICES = list(measures$SampleID_num),
                     FUN = function(x){
                       print(paste('SampleID_num:',unique(x$SampleID_num)))
                       Start = list(JmaxRef=JmaxRef,RdRef=RdRef,VcmaxRef=VcmaxRef)
                       modify.init=TRUE
                       if(!is.null(result_Ac_Aj[[as.character(unique(x$SampleID_num))]][[2]])){
                         without_Tp_Start=result_Ac_Aj[[as.character(unique(x$SampleID_num))]][[2]]@coef
                         modify.init=FALSE
                         for(value in names(Start)){Start[[value]]=without_Tp_Start[value]}
                       }
                       Start[['TpRef']]=Start[['VcmaxRef']]/10
                       
                       f.fitting(measures = x,Start = Start,
                                 param=param,id.name = 'SampleID_num',modify.init =modify.init )})
  dev.off()
  
  ##Fitting of the curve using the Ac and Ap limitations (without Aj)
  pdf(file = 'Outputs/2_ACi_fitting_Ac_Ap.pdf')
  result_Ac_Ap=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){
                    print(paste('SampleID_num:',unique(x$SampleID_num)))
                    Start = list(RdRef=RdRef,VcmaxRef=VcmaxRef)
                    modify.init=TRUE
                    if(!is.null(result_Ac_Aj_Ap[[as.character(unique(x$SampleID_num))]][[2]])){
                      without_Aj_Start=result_Ac_Aj_Ap[[as.character(unique(x$SampleID_num))]][[2]]@coef
                      modify.init=FALSE
                      for(value in names(Start)){Start[[value]]=without_Aj_Start[value]}
                    }
                    Start[['TpRef']]=Start[['VcmaxRef']]/8
                    
                    f.fitting(measures = x,Start = Start,
                              param=param,id.name = 'SampleID_num',modify.init =modify.init )})
  dev.off()
  
  ## Extracting the fitting metrics for each model and each curve
  res_nlf_Ac_Aj=as.data.frame(t(sapply(result_Ac_Aj,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else { coefs=rep(NA,4)
                                  std_dev=rep(NA,4)
                                  AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Aj$SampleID_num=row.names(res_nlf_Ac_Aj)
  
  res_nlf_Ac=as.data.frame(t(sapply(result_Ac,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,3)
                                std_dev=rep(NA,3)
                                AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac$SampleID_num=row.names(res_nlf_Ac)
  
  res_nlf_Ac_Aj_Ap=as.data.frame(t(sapply(result_Ac_Aj_Ap,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,5)
                                std_dev=rep(NA,5)
                                AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Aj_Ap$SampleID_num=row.names(res_nlf_Ac_Aj_Ap)
  
  res_nlf_Ac_Ap=as.data.frame(t(sapply(result_Ac_Ap,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,4)
      std_dev=rep(NA,4)
      AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Ap$SampleID_num=row.names(res_nlf_Ac_Ap)
  
  res_nlf_Ac$JmaxRef=NA
  res_nlf_Ac$TpRef=NA
  res_nlf_Ac$StdError_JmaxRef=NA
  res_nlf_Ac$StdError_TpRef=NA
  res_nlf_Ac_Aj$TpRef=NA
  res_nlf_Ac_Aj$StdError_TpRef=NA
  res_nlf_Ac_Aj$model='Ac_Aj'
  res_nlf_Ac$model='Ac'
  res_nlf_Ac_Aj_Ap$model='Ac_Aj_Ap'
  
  res_nlf_Ac_Ap$JmaxRef=NA
  res_nlf_Ac_Ap$StdError_JmaxRef=NA
  res_nlf_Ac_Ap$model='Ac_Ap'
  
  res_nlf_Ac_Aj=res_nlf_Ac_Aj[,colnames(res_nlf_Ac_Aj_Ap)]
  res_nlf_Ac=res_nlf_Ac[,colnames(res_nlf_Ac_Aj_Ap)]
  res_nlf_Ac_Ap=res_nlf_Ac_Ap[,colnames(res_nlf_Ac_Aj_Ap)]
  
  ## Finding the best model (Ac or Ac_Aj or Ac_Aj_Ap or Ac_Ap according to the AIC criterion)
  
  Bilan=res_nlf_Ac_Aj
  Bilan[which(res_nlf_Ac$AIC<res_nlf_Ac_Aj$AIC),]=res_nlf_Ac[which(res_nlf_Ac$AIC<res_nlf_Ac_Aj$AIC),]
  Bilan[which(res_nlf_Ac_Aj_Ap$AIC<Bilan$AIC),]=res_nlf_Ac_Aj_Ap[which(res_nlf_Ac_Aj_Ap$AIC<Bilan$AIC),]
  Bilan[which(res_nlf_Ac_Ap$AIC<Bilan$AIC),]=res_nlf_Ac_Ap[which(res_nlf_Ac_Ap$AIC<Bilan$AIC),]
  colnames(Bilan)=c("sigma","JmaxRef","VcmaxRef","TpRef","RdRef","StdError_sigma","StdError_JmaxRef","StdError_VcmaxRef","StdError_TpRef","StdError_RdRef","AIC","Tleaf","SampleID_num","model") 
  Bilan$Vcmax_method="A-Ci curve"
  
  ### Creating a pdf with the best model for each curve
  
  pdf(file = 'Outputs/2_ACi_fitting_best_model.pdf')
  A_limitations=by(data = measures,INDICES = list(measures$SampleID_num),
     FUN = function(x){
       param_leg=Bilan[Bilan$SampleID_num==unique(x$SampleID_num),c('VcmaxRef','JmaxRef','RdRef','TpRef')]
       param_leg[is.na(param_leg)]=9999
       param['VcmaxRef']=param_leg['VcmaxRef']
       param['JmaxRef']=param_leg['JmaxRef']
       param['RdRef']=param_leg['RdRef']
       param['TpRef']=param_leg['TpRef']
       f.plot(measures=x,list_legend = as.list(param_leg),param = param,name =unique(x$SampleID_num))
       simu=f.ACi(PFD=x$Qin,Tleaf = x$Tleaf,Ci = x$Ci,param = param)
       ## Number of points limited by Ac, Aj and Ap
       n_Ac=sum(simu$Ac<simu$Aj&simu$Ac<simu$Ap, na.rm = TRUE)
       n_Aj=sum(simu$Aj<simu$Ac&simu$Aj<simu$Ap, na.rm = TRUE)
       n_Ap=sum(simu$Ap<simu$Ac&simu$Ap<simu$Aj, na.rm = TRUE)
       n_lim=c(SampleID_num=unique(x$SampleID_num),n_Ac=n_Ac,n_Aj=n_Aj,n_Ap=n_Ap)
       return(n_lim)
     }
  )
  dev.off()
  
  A_limitations=do.call(rbind.data.frame, A_limitations)
  colnames(A_limitations)=c("SampleID_num","n_Ac","n_Aj","n_Ap")
  Bilan=cbind.data.frame(Bilan,A_limitations[,c("n_Ac","n_Aj","n_Ap")])
  ## I remove estimates of Ac, Aj or Ap that are based on too few points. 
  Bilan[Bilan$n_Ac<=2,c("VcmaxRef","StdError_VcmaxRef")]=NA
  Bilan[Bilan$n_Aj<2,c("JmaxRef","StdError_JmaxRef")]=NA
  Bilan[Bilan$n_Ap<2,c("TpRef","StdError_TpRef")]=NA
  return(Bilan)
}



#' @param measures ACi curve dataframe with at least the columns A, Ci, Qin and Tleaf (in K)
#' @param param See f.make.param for details. Determine the parameters used for the fitting
#' @param VcmaxRef Value of VcmaxRef used to initialize the fitting procedure
#' @param VqmaxRef Value of VqmaxRef used to initialize the fitting procedure
#' @param RdRef Value of RdRef used to initialize the fitting procedure
#' @param TpRef Value of TpRef used to initialize the fitting procedure
#'
#' @return Return a dataframe (Bilan) with for each curve in row the estimated parameters as well as the model used (Ac, Ac and Aj, Ac and Aj and Ap)
#' @export
#'
#' @examples
f.fit_ACi_JB<-function(measures,param,VcmaxRef=50, VqmaxRef=100, RdRef = 2, TpRef= 5){
  ## Basic checks
  if(any(measures$Ci < 0)) {stop("Ci contains negative values")}
  if(any(is.na(c(measures$A,measures$Ci,measures$Tleaf,measures$Qin,measures$SampleID_num)))) {stop("Ci, A, Tleaf or Qin have NA values")}
  if(any(measures$Tleaf < 273)) {stop("Check Tleaf values. Should be in Kelvin")}
  if(is.null(measures$SampleID_num)|is.null(measures$Ci)|is.null(measures$Tleaf)|is.null(measures$Qin)|is.null(measures$A)){stop("Your dataframe measures does not contain the colums SampleID_num or A or Tleaf or Ci")}
  n_points_curves=tapply(X=measures$SampleID_num,INDEX = measures$SampleID_num,FUN = function(x) length(x))
  min_Ci_curve=tapply(X=measures$Ci,INDEX = measures$SampleID_num,FUN = min)
  if(any(n_points_curves<3)){stop("Some of your A-Ci curves have less than 3 points")}
  if(any(min_Ci_curve>270)){stop("Some of your A-Ci curves have a minimum Ci above 270 ppm meaning that Vcmax cant be estimated")}
  
  param[['TpRef']]=9999
  
  ## Fitting of the ACi curve using Ac and Aj limitation
  pdf(file = 'Outputs/2_JB_ACi_fitting_Ac_Aj.pdf')
  result_Ac_Aj=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){print(paste('SampleID_num:',unique(x$SampleID_num)))
                    f.fitting_JB(measures = x,Start = list(VqmaxRef=VqmaxRef,RdRef=RdRef,VcmaxRef=VcmaxRef),
                              param=param,id.name = 'SampleID_num')})
  dev.off()
  
  param[['VqmaxRef']]=9999
  ## Fitting of the ACi curve using only the Ac limitation
  pdf(file = 'Outputs/2_JB_ACi_fitting_Ac.pdf')
  result_Ac=by(data = measures,INDICES = list(measures$SampleID_num),
               FUN = function(x){print(paste('SampleID_num:',unique(x$SampleID_num)))
                 f.fitting_JB(measures = x,Start = list(RdRef=RdRef,VcmaxRef=VcmaxRef),
                           param=param,id.name = 'SampleID_num')})
  dev.off()
  
  ##Fitting of the curve using the Ac, Aj and Ap limitations
  pdf(file = 'Outputs/2_JB_ACi_fitting_Ac_Aj_Ap.pdf')
  result_Ac_Aj_Ap=by(data = measures,INDICES = list(measures$SampleID_num),
                     FUN = function(x){
                       print(paste('SampleID_num:',unique(x$SampleID_num)))
                       Start = list(VqmaxRef=VqmaxRef,RdRef=RdRef,VcmaxRef=VcmaxRef)
                       modify.init=TRUE
                       if(!is.null(result_Ac_Aj[[as.character(unique(x$SampleID_num))]][[2]])){
                         without_Tp_Start=result_Ac_Aj[[as.character(unique(x$SampleID_num))]][[2]]@coef
                         modify.init=FALSE
                         for(value in names(Start)){Start[[value]]=without_Tp_Start[value]}
                       }
                       Start[['TpRef']]=Start[['VcmaxRef']]/10
                       
                       f.fitting_JB(measures = x,Start = Start,
                                 param=param,id.name = 'SampleID_num',modify.init =modify.init )})
  dev.off()
  
  ##Fitting of the curve using the Ac and Ap limitations (without Aj)
  pdf(file = 'Outputs/2_JB_ACi_fitting_Ac_Aj_Ap.pdf')
  result_Ac_Ap=by(data = measures,INDICES = list(measures$SampleID_num),
                  FUN = function(x){
                    print(paste('SampleID_num:',unique(x$SampleID_num)))
                    Start = list(RdRef=RdRef,VcmaxRef=VcmaxRef)
                    modify.init=TRUE
                    if(!is.null(result_Ac_Aj_Ap[[as.character(unique(x$SampleID_num))]][[2]])){
                      without_Aj_Start=result_Ac_Aj_Ap[[as.character(unique(x$SampleID_num))]][[2]]@coef
                      modify.init=FALSE
                      for(value in names(Start)){Start[[value]]=without_Aj_Start[value]}
                    }
                    Start[['TpRef']]=Start[['VcmaxRef']]/8
                    
                    f.fitting_JB(measures = x,Start = Start,
                              param=param,id.name = 'SampleID_num',modify.init =modify.init )})
  dev.off()
  
  
  ## Extracting the fitting metrics for each model and each curve
  res_nlf_Ac_Aj=as.data.frame(t(sapply(result_Ac_Aj,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else { coefs=rep(NA,4)
      std_dev=rep(NA,4)
      AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Aj$SampleID_num=row.names(res_nlf_Ac_Aj)
  
  res_nlf_Ac=as.data.frame(t(sapply(result_Ac,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,3)
      std_dev=rep(NA,3)
      AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac$SampleID_num=row.names(res_nlf_Ac)
  
  res_nlf_Ac_Aj_Ap=as.data.frame(t(sapply(result_Ac_Aj_Ap,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,5)
      std_dev=rep(NA,5)
      AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Aj_Ap$SampleID_num=row.names(res_nlf_Ac_Aj_Ap)
  
  res_nlf_Ac_Ap=as.data.frame(t(sapply(result_Ac_Ap,FUN = function(x){
    if(!is.null(x[[2]])){
      coefs=x[[2]]@coef
      std_dev=sqrt(diag(x[[2]]@vcov))
      names(std_dev)=paste('StdError',names(std_dev),sep='_')
      AICcurve=AIC(x[[2]])}else {coefs=rep(NA,4)
      std_dev=rep(NA,4)
      AICcurve=NA}
    return(c(coefs,std_dev,AIC=AICcurve,Tleaf=x[[3]]['Tleaf']))}
  )))
  res_nlf_Ac_Ap$SampleID_num=row.names(res_nlf_Ac_Ap)
  
  res_nlf_Ac$VqmaxRef=NA
  res_nlf_Ac$TpRef=NA
  res_nlf_Ac$StdError_VqmaxRef=NA
  res_nlf_Ac$StdError_TpRef=NA
  res_nlf_Ac_Aj$TpRef=NA
  res_nlf_Ac_Aj$StdError_TpRef=NA
  res_nlf_Ac_Aj$model='Ac_Aj'
  res_nlf_Ac$model='Ac'
  res_nlf_Ac_Aj_Ap$model='Ac_Aj_Ap'
  
  res_nlf_Ac_Ap$VqmaxRef=NA
  res_nlf_Ac_Ap$StdError_VqmaxRef=NA
  res_nlf_Ac_Ap$model='Ac_Ap'
  
  res_nlf_Ac_Aj=res_nlf_Ac_Aj[,colnames(res_nlf_Ac_Aj_Ap)]
  res_nlf_Ac=res_nlf_Ac[,colnames(res_nlf_Ac_Aj_Ap)]
  res_nlf_Ac_Ap=res_nlf_Ac_Ap[,colnames(res_nlf_Ac_Aj_Ap)]
  
  ## Finding the best model (Ac or Ac_Aj or Ac_Aj_Ap or Ac_Ap according to the AIC criterion)
  
  Bilan=res_nlf_Ac_Aj
  Bilan[which(res_nlf_Ac$AIC<res_nlf_Ac_Aj$AIC),]=res_nlf_Ac[which(res_nlf_Ac$AIC<res_nlf_Ac_Aj$AIC),]
  Bilan[which(res_nlf_Ac_Aj_Ap$AIC<Bilan$AIC),]=res_nlf_Ac_Aj_Ap[which(res_nlf_Ac_Aj_Ap$AIC<Bilan$AIC),]
  Bilan[which(res_nlf_Ac_Ap$AIC<Bilan$AIC),]=res_nlf_Ac_Ap[which(res_nlf_Ac_Ap$AIC<Bilan$AIC),]
  colnames(Bilan)=c("sigma","VcmaxRef","VqmaxRef","TpRef","RdRef","StdError_sigma","StdError_VcmaxRef","StdError_VqmaxRef","StdError_TpRef","StdError_RdRef","AIC","Tleaf","SampleID_num","model") 
  Bilan$Vcmax_method="A-Ci curve"
  
  ### Creating a pdf with the best model for each curve
  
  pdf(file = 'Outputs/2_JB_ACi_fitting_best_model.pdf')
  A_limitations=by(data = measures,INDICES = list(measures$SampleID_num),
     FUN = function(x){
       param_leg=Bilan[Bilan$SampleID_num==unique(x$SampleID_num),c('VcmaxRef','VqmaxRef','RdRef','TpRef')]
       param_leg[is.na(param_leg)]=9999
       param['VcmaxRef']=param_leg['VcmaxRef']
       param['VqmaxRef']=param_leg['VqmaxRef']
       param['RdRef']=param_leg['RdRef']
       param['TpRef']=param_leg['TpRef']
       f.plot_JB(measures=x,list_legend = as.list(param_leg),param = param,name =unique(x$SampleID_num))
       simu=f.ACi_JB(PFD=x$Qin,Tleaf = x$Tleaf,Ci = x$Ci,param = param)
       ## Number of points limited by Ac, Aj and Ap
       n_Ac=sum(simu$Ac<simu$Aj&simu$Ac<simu$Ap, na.rm = TRUE)
       n_Aj=sum(simu$Aj<simu$Ac&simu$Aj<simu$Ap, na.rm = TRUE)
       n_Ap=sum(simu$Ap<simu$Ac&simu$Ap<simu$Aj, na.rm = TRUE)
       n_lim=c(SampleID_num=unique(x$SampleID_num),n_Ac=n_Ac,n_Aj=n_Aj,n_Ap=n_Ap)
       return(n_lim)
     }
  )
  dev.off()
  A_limitations=do.call(rbind.data.frame, A_limitations)
  colnames(A_limitations)=c("SampleID_num","n_Ac","n_Aj","n_Ap")
  Bilan=cbind.data.frame(Bilan,A_limitations[,c("n_Ac","n_Aj","n_Ap")])
  ## I remove estimates of Ac, Aj or Ap that are based on too few points. 
  Bilan[Bilan$n_Ac<=2,c("VcmaxRef","StdError_VcmaxRef")]=NA
  Bilan[Bilan$n_Aj<2,c("VqmaxRef","StdError_VqmaxRef")]=NA
  Bilan[Bilan$n_Ap<2,c("TpRef","StdError_TpRef")]=NA
  colnames(Bilan)=paste("JB",colnames(Bilan),sep="_")
  return(Bilan)
}




