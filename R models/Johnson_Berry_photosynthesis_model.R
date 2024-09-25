################################################################################
###   This code contains various functions to simulate leaf gas exchange     ###
###   with the FvCB photosynthesis model modified by Johnson and Berry.      ###
################################################################################

#' @title Function to generate a list of parameters for the Johnson and Berry (2021) photosynthesis model
#' @param VcmaxRef Maximum rate of Rubisco for carboxylation micromol.m-2.s-1
#' @param VcmaxHa Energy of activation for Vcmax in J.mol-1
#' @param VcmaxHd Energy of desactivation for Vcmax in J.mol-1
#' @param VcmaxS Entropy term for Vcmax in J.mol-1.K-1
#' @param RdRef Respiration value at the reference temperature
#' @param RdHa Energie of activation for Rd in J.mol-1
#' @param RdHd Energy of desactivation for Rd in J.mol-1
#' @param RdS Entropy term for Rd in J.mol-1.K-1
#' @param KcRef Michaelis-Menten constant of Rubisco for CO2 at the reference temperature in micromol.mol-1
#' @param KcHa Energy of activation for Kc in J.mol-1
#' @param KoRef ichaelis-Menten constant of Rubisco for CO2 at the reference temperature in milimol.mol-1
#' @param KoHa Energy of activation for Ko in J.mol-1
#' @param GstarRef CO2 compensation point in absence of respiration in micromol.mol-1
#' @param GstarHa Enthalpie of activation for Gstar in J.mol-1
#' @param VqmaxRef Maximum activity of cytochrom b6f in micro mol electron-1 m-2 s-1
#' @param VqmaxHa Energy of activation for Vqmax in J.mol-1
#' @param VqmaxHd Energy of desactivation for Vqmax in J.mol-1. Note that Johnson et al. 2021 advice to use the desactivation terms on etal and etac
#' @param VqmaxS Entropy term for Vqmax in J.mol-1.K-1.  Note that Johnson et al. 2021 advice to use the desactivation terms on etal and etac
#' @param abso Absorptance of the leaf (0 - 1, unitless)
#' @param Beta PSII fraction of leaf absorptance
#' @param etacRef Coupling efficiency for CEF1 at the reference temperature (unitless)
#' @param etacHd Energy of desactivation for etac in J.mol-1
#' @param etacS Entropy term for etac in J.mol-1
#' @param etalRef Coupling efficiency for LEF at the reference temperature (unitless)
#' @param etalHd Energy of desactivation for etal in J.mol-1
#' @param etalS Entropy term for etal in J.mol-1
#' @param phi1max Maximum photochemical yield of PS1, calculated using eq 30a, 41 c and Table 1
#' @param thetacj Smoothing factor between Ac and Aj
#' @param thetaip Smoothing factor between Ai (smoothed minimum between Ac and Aj) and Ap
#' @param model.diff Type of model to simulate the transport between the leaf surface and the sub stomatal cavity. Chose between 'Fick', 'vCF' and 'MSWF'
#' @param gcw Cuticular conductance. THis parameter is only used when model.diff corresponds to 'MSWF'. See Lamour et al. 2021 and Marquez et al. 2021
#' @param model.gs Type of conductance model (USO, USO_simpl or BWB).
#' @param g0 Constant of the conductance model representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter of the conductance model
#' @param Patm Atmospheric pressure in kPa
#' @return List of parameters that can be used with the JB model of photosynthesis
#' @references 
#' #' @references Bernacchi, C.J., Singsaas, E.L., Pimentel, C., Portis Jr, A.R. and Long, S.P. (2001), Improved temperature response functions for models of Rubisco limited photosynthesis. Plant, Cell & Environment, 24: 253-259. doi:10.1111/j.1365-3040.2001.00668.x
#' Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x
#' Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D. (1995). Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies. Plant, Cell & Environment, 18(10), 1183-1200.
#' Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.
#' M?rquez, D.A., Stuart-Williams, H. & Farquhar, G.D. An improved theory for calculating leaf gas exchange more precisely accounting for small fluxes. Nat. Plants 7, 317-326 (2021). https://doi.org/10.1038/s41477-021-00861-w
#' Lamour, J., Davidson, K.J., Ely, K.S., Li, Q., Serbin, S.P. and Rogers, A. (2021), New calculations for photosynthesis measurement systems: what's the impact for physiologists and modelers?. New Phytol. https://doi.org/10.1111/nph.17762
#' @export
#'
#' @examples
f.make.param_JB<-function(R=NA,TRef=NA,VcmaxRef=NA,
                       VcmaxHa	=NA,
                       VcmaxHd	=NA,
                       VcmaxS	=NA,
                       VqmaxRef=NA,
                       VqmaxHa=	NA,
                       VqmaxHd=NA,
                       VqmaxS=NA,
                       TpRef=NA,
                       TpHa=NA,
                       TpHd=NA,
                       TpS=NA,
                       abso=NA,
                       Beta=NA,
                       GstarRef=NA, 
                       GstarHa=NA,
                       etacRef=NA,
                       etacHd=NA,
                       etacS=NA,
                       etalRef=NA,
                       etalHd=NA,
                       etalS=NA,
                       phi1max=NA,
                       KcRef=	NA,
                       KcHa=	NA,
                       KcQ10=NA,
                       KoRef=	NA,
                       KoHa=	NA,
                       RdRef=NA,
                       RdHa=	NA,
                       RdHd= NA,
                       RdS=NA,
                       O2=NA,
                       thetacj=NA,
                       thetaip=NA,
                       model.diff=NA,
                       gcw=NA,
                       model.gs=NA,
                       g0=NA,
                       g1=NA,
                       Patm=NA
                       ){
  param=list(R=8.314,TRef=298.16,VcmaxRef=100,VcmaxHa	=65330,VcmaxHd	=149250,VcmaxS	=485,
             VqmaxRef=350,VqmaxHa=	43540,VqmaxHd=0, VqmaxS=0, TpRef=9999,TpHa=53100,TpHd=150650,TpS=490,
             abso=0.85,Beta=0.52,GstarRef=42.75, GstarHa	=37830,etacRef=1,etacHd=152040,etacS=495,etalRef=0.75,etalHd=152040,etalS=495,phi1max=0.96,
             KcRef=	404.9,KcHa=	79430,KoRef=	278.4,KoHa=	36380,
             RdRef=1.5,RdHa=46390,RdHd=150650,RdS=490,O2=210,thetacj=0.999,thetaip=0.999,model.diff=0,gcw= 10*10^-3,model.gs=0,g0=0.01,g1=4.1,Patm=101)
  if(!is.na(model.diff)&model.diff=='vCF'){model.diff=1}else if(!is.na(model.diff)&model.diff=='MSWF'){model.diff=2}else if (!is.na(model.diff)&!model.diff%in%c('Fick','vCF','MSWF')){print('Unknown model.diff, should be Fick, vCF (mind the lower case) or MSWF')}
  if(!is.na(model.gs)&model.gs=="USO"){model.gs=0}else if(!is.na(model.gs)&model.gs=="USO_simpl"){model.gs=1}else if(!is.na(model.gs)&model.gs=="BWB"){model.gs=2}else if(!is.na(model.gs)&model.gs=="Nonlinear"){model.gs=3}else if(!is.na(model.gs)&!model.gs%in%c("USO","USO_simpl","BWB","Nonlinear")){print("Unknown model.gs")}

  param_fun=list(R=R,TRef=TRef,VcmaxRef=VcmaxRef,VcmaxHa=VcmaxHa,VcmaxHd	=VcmaxHd,VcmaxS	=VcmaxS,
             VqmaxRef=VqmaxRef,VqmaxHa=	VqmaxHa,VqmaxHd=VqmaxHd, VqmaxS=VqmaxS,TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
             abso=abso,Beta=Beta,GstarRef=GstarRef, GstarHa	=GstarHa,etacRef=etacRef,etacHd=etacHd,etacS=etacS,etalRef=etalRef,etalHd=etalHd,etalS=etalS,phi1max=phi1max,
             KcRef=	KcRef,KcHa=	KcHa,KoRef=	KoRef,KoHa=KoHa,
             RdRef=RdRef,RdHa=RdHa,RdHd=RdHd,RdS=RdS,O2=O2,thetacj=thetacj,thetaip=thetaip,model.diff=model.diff,gcw= gcw,model.gs=model.gs,g0=g0,g1=g1,Patm=Patm)
  modified=which(lapply(X=param_fun,FUN = is.na)==FALSE)
  if(length(modified)>0){
    for(i in 1: length(modified)){param[names(modified[i])]=param_fun[modified[i]]}
  }
  param_fun[names(param)]=param
  return(param_fun)
}

f.desactivation<-function(PRef,Hd,s,Tleaf,TRef=298.16,R=8.314){
  P=PRef*(1+exp((s*TRef-Hd)/(R*TRef)))/(1+exp((s*Tleaf-Hd)/(R*Tleaf)))
  return(P)
}

f.ACi_JB<-function(PFD,Tleaf,Ci,param){
  Kc=f.arrhenius(PRef = param[['KcRef']],Ha = param[['KcHa']],Tleaf=Tleaf)
  Ko=f.arrhenius(PRef = param[['KoRef']],Ha = param[['KoHa']],Tleaf=Tleaf)
  Gstar=f.arrhenius(PRef = param[['GstarRef']],Ha = param[['GstarHa']],Tleaf = Tleaf)
  
  ## Mitochondrial respiration
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],Ha = param[['RdHa']],Hd = param[['RdHd']],s=param[['RdS']],Tleaf = Tleaf)
  
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],Ha = param[['VcmaxHa']],Hd = param[['VcmaxHd']],s = param[['VcmaxS']],Tleaf = Tleaf)
  etal=f.desactivation(PRef=param[['etalRef']],Hd = param[['etalHd']],s = param[['etalS']],Tleaf=Tleaf)
  etac=f.desactivation(PRef=param[['etacRef']],Hd = param[['etacHd']],s = param[['etacS']],Tleaf=Tleaf)
  
  Vqmax=f.modified.arrhenius(PRef=param[['VqmaxRef']],Ha = param[['VqmaxHa']],Hd = param[['VqmaxHd']],s=param[['VqmaxS']],Tleaf = Tleaf)
  
  eta=1-etal/etac+(3+7*Gstar/Ci)/((4+8*Gstar/Ci)*etac) ## Eq 15 c
  
  Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
  
  alpha1=1-param[['Beta']]

  JP700j=Vqmax*PFD/(Vqmax/(param[['abso']]*alpha1*param[['phi1max']])+PFD) ## Eq 30
  JP680j=JP700j/eta ## Eq 30
  
  ## Calculation of PFDsat
  #static case
  #K=Vcmax*(4*Ci+8*Gstar)/(Ci+Kc*(1+param[['O2']]/Ko))
  #PFDsat=(K*Vqmax/(param[['abso']]*alpha1*param[['phi1max']]))/(eta^-1*Vqmax-K)
  
  #Dynamic case
  #PFDsat=-(sqrt(((Abs^2*Kp2^2*Ku2^2+(2*Abs^2*Kf+2*Abs^2*Kd)*Kp2^2*Ku2+(Abs^2*Kf^2+2*Abs^2*Kd*Kf+Abs^2*Kd^2)*Kp2^2)*Q^2*eta^2+(((2*Abs*Kf+2*Abs*Kd)*Kp2-2*Abs*Kp2^2)*Ku2^2+(4*Abs*Kf^2+8*Abs*Kd*Kf+4*Abs*Kd^2)*Kp2*Ku2+(2*Abs*Kf^2+4*Abs*Kd*Kf+2*Abs*Kd^2)*Kp2^2+(2*Abs*Kf^3+6*Abs*Kd*Kf^2+6*Abs*Kd^2*Kf+2*Abs*Kd^3)*Kp2)*Q*Vqmax*eta+((Kp2^2+(2*Kf+2*Kd)*Kp2+Kf^2+2*Kd*Kf+Kd^2)*Ku2^2+((2*Kf+2*Kd)*Kp2^2+(4*Kf^2+8*Kd*Kf+4*Kd^2)*Kp2+2*Kf^3+6*Kd*Kf^2+6*Kd^2*Kf+2*Kd^3)*Ku2+(Kf^2+2*Kd*Kf+Kd^2)*Kp2^2+(2*Kf^3+6*Kd*Kf^2+6*Kd^2*Kf+2*Kd^3)*Kp2+Kf^4+4*Kd*Kf^3+6*Kd^2*Kf^2+4*Kd^3*Kf+Kd^4)*Vqmax^2)*phi1P_max^2+((2*Abs*Kp2^2*Ku2^2+(4*Abs*Kf+4*Abs*Kd)*Kp2^2*Ku2+(2*Abs*Kf^2+4*Abs*Kd*Kf+2*Abs*Kd^2)*Kp2^2)*Q*Vqmax*eta^2+((2*Kp2^2+(2*Kf+2*Kd)*Kp2)*Ku2^2+((4*Kf+4*Kd)*Kp2^2+(4*Kf^2+8*Kd*Kf+4*Kd^2)*Kp2)*Ku2+(2*Kf^2+4*Kd*Kf+2*Kd^2)*Kp2^2+(2*Kf^3+6*Kd*Kf^2+6*Kd^2*Kf+2*Kd^3)*Kp2)*Vqmax^2*eta)*phi1P_max+(Kp2^2*Ku2^2+(2*Kf+2*Kd)*Kp2^2*Ku2+(Kf^2+2*Kd*Kf+Kd^2)*Kp2^2)*Vqmax^2*eta^2)+K*(Kd*(-2*Ku2*phi1P_max^2-4*Kf*phi1P_max^2+Kp2*(-2*phi1P_max^2-2*eta*phi1P_max))+Kf*(Kp2*(-2*phi1P_max^2-2*eta*phi1P_max)-2*Ku2*phi1P_max^2)-2*Kf^2*phi1P_max^2-2*Kd^2*phi1P_max^2-2*Kp2*Ku2*eta*phi1P_max)+Kd*(-Ku2*Vqmax*phi1P_max-2*Kf*Vqmax*phi1P_max+Kp2*Vqmax*(-phi1P_max-eta))+Kf*(Kp2*Vqmax*(-phi1P_max-eta)-Ku2*Vqmax*phi1P_max)-Kf^2*Vqmax*phi1P_max-Kd^2*Vqmax*phi1P_max+Kp2*Ku2*Vqmax*(-phi1P_max-eta))/(Abs*(Kp2*Ku2*eta*phi1P_max+Kf*Kp2*eta*phi1P_max+Kd*Kp2*eta*phi1P_max))
  
  Aj=JP680j/(4+8*Gstar/Ci)*(1-Gstar/Ci) ## Eq 31 without Rd
  Ac=Vcmax*Ci/(Ci+Kc*(1+param[['O2']]/Ko))*(1-Gstar/Ci) ## Eq 32 without Rd
  
  JP680c=Ac*(4+8*Gstar/Ci)/(1-Gstar/Ci) ##Eq 33a and 32
  JP700c=JP680c/eta
  
  JP700=pmin(JP700c,JP700j)## Johnson and Berry use a smoothing function, here I keep simple for now
  JP680=pmin(JP680c,JP680j)
  
  Ap=3*Tp
  
  ##
  Ai=f.smooth(A1 = Ac,A2 = Aj,theta=param[['thetacj']])*as.numeric(Ci>Gstar)+f.smooth(A1 = Ac,A2 = Aj,theta=param[['thetacj']],root = 2)*as.numeric(Ci<=Gstar)
  Ag=f.smooth(A1=Ai,A2=Ap,theta=param[['thetaip']])
  An=Ag-Rd
  
  result=data.frame(A=An,Ag=Ag,Aj=Aj-Rd,Ac=Ac-Rd,Ap=Ap-Rd)
  return(result)
}


#' @title Fitting function for photosynthesis data (light curve or ACi curve)
#' @description Function to fit model to data. The parameters to fit have to be described in the list Start.
#' All the other parameters of the f.ACi functions have to be in param. If the parameters from Start are repeated in param, the later one will be ignored.
#' This function uses two methods to fit the data. First by minimizing the residual sum-of-squares of the residuals and then by maximizing the likelihood function. The first method is more robust but the second one allows to calculate the confident interval of the parameters.
#' @param measures Data frame of measures obtained from gas exchange analyser with at least the columns A, Ci, Qin and Tleaf (in K). If RHs, Tair, Patm, VPDleaf are also present, there mean will be added in the output, but those columns are not needed to estimate the parameters
#' @param id.name Name of the colums in measures with the identifier for the curve.
#' @param Start List of parameters to fit with their initial values.
#' @param param See f.make.param() for details.
#' @param modify.init TRUE or FALSE, allows to modify the Start values before fitting the data
#' @param do.plot TRUE or FALSE, plot data and fitted curves ?
#' @return Return a list with 3 components, 1 the result of the optim function which is used to estimate the parameters, 2 the output of the function bbmle, 3 the mean variable of the environment during the measurement
#' @export
#'
#' @examples ##Simulation of a CO2 curve
#' data=data.frame(Tleaf=rep(300,20),
#' Ci=seq(40,1500,75),Qin=rep(2000,20),Tair=300,RHs=70,VPDleaf=2,Patm=101,A=f.ACi_JB(PFD=2000,Tleaf=300,Ci=seq(40,1500,75),
#' param=f.make.param_JB())$A+rnorm(n = 20,mean = 0,sd = 0.5))
#'
#' f.fitting(measures=data,id.name=NULL,Start=list(JmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param())
f.fitting_JB<-function(measures,id.name=NULL,Start=list(VqmaxRef=90,VcmaxRef=70,RdRef=1),param=f.make.param_JB(),modify.init=TRUE,do.plot=TRUE,type='ACi'){
  Fixed=param[!names(param)%in%names(Start)]
  if(modify.init){
    if('VqmaxRef'%in%names(Start)){Start[['VqmaxRef']]=f.arrhenius.inv(P = 6*1.6*(max(measures$A,na.rm=TRUE)+1),Ha = param[['VqmaxHa']],Tleaf = mean(measures$Tleaf,na.rm=TRUE),TRef = param[['TRef']],R = param[['R']])}
    if('VqmaxRef'%in%names(Start)&'VcmaxRef'%in%names(Start)){Start[['VcmaxRef']]=Start[['VqmaxRef']]/2.5}
    grille=expand.grid(lapply(X = Start,FUN = function(x){x*c(0.2,1,2)}))
    grille.list=apply(X=grille,MARGIN = 1,FUN=as.list)
    value=9999999
    l_param=0
    for(l in 1:nrow(grille)){
      MoindresCarres=optim(par=grille.list[[l]],fn=f.SumSq_JB,data=measures,Fixed=Fixed)
      if(!is.null(MoindresCarres)&MoindresCarres$value<value){value=MoindresCarres$value;l_param=l}
    }
    Start=grille.list[[l_param]]
  }
  
  if(is.null(id.name)){name=''}else{name=unique(measures[,id.name])}
  MoindresCarres<-Estimation2<-NULL
  try({
    MoindresCarres<-optim(par=Start,fn=f.SumSq_JB,data=measures,Fixed=Fixed)
    print(MoindresCarres)
    print(paste('sd',sqrt(MoindresCarres$value/NROW(measures))))
    Start$sigma=sqrt(MoindresCarres$value/NROW(measures))
    for(l.name in names(MoindresCarres$par)){Start[l.name]=MoindresCarres$par[[l.name]]}
    for(l.name in names(MoindresCarres$par)){param[l.name]=MoindresCarres$par[[l.name]]}
    #if(do.plot){f.plot(measures=measures,name=name,param =param,list_legend = Start,type=type)}
  })
  
  try({
    Estimation2=mle2(minuslogl = f.MinusLogL_JB,start = Start,fixed = Fixed,data = list(data=measures),method = "Nelder-Mead")
    print(summary(Estimation2))
    #conf=confint(Estimation2)
    #print(conf)
    for(i in names(Estimation2@coef[names(Estimation2@coef)%in%names(param)])){param[i]=Estimation2@coef[i]}
    if(do.plot){f.plot_JB(measures=measures,name=name,param =param,list_legend = as.list(Estimation2@coef),type=type)}
  })
  Envir=NA
  try({
    Envir=c(Tair=mean(measures$Tair,na.rm=TRUE),Tleaf=mean(measures$Tleaf,na.rm=TRUE),RHs=mean(measures$RHs,na.rm=TRUE),VPDleaf=mean(measures$VPDleaf,na.rm=TRUE),Qin=mean(measures$Qin,na.rm=TRUE),Patm=mean(measures$Patm,na.rm=TRUE))
  })
  return(list(MoindresCarres,Estimation2,Envir))
}


f.MinusLogL_JB<-function(data,sigma,R=NA,TRef=NA,VcmaxRef=NA,
                         VcmaxHa	=NA,
                         VcmaxHd	=NA,
                         VcmaxS	=NA,
                         VqmaxRef=NA,
                         VqmaxHa=	NA,
                         VqmaxHd= NA,
                         VqmaxS= NA,
                         TpRef=NA,
                         TpHa=NA,
                         TpHd=NA,
                         TpS=NA,
                         abso=NA,
                         Beta=NA,
                         GstarRef=NA, 
                         GstarHa=NA,
                         etacRef=NA,
                         etacHd=NA,
                         etacS=NA,
                         etalRef=NA,
                         etalHd=NA,
                         etalS=NA,
                         phi1max=NA,
                         KcRef=	NA,
                         KcHa=	NA,
                         KcQ10=NA,
                         KoRef=	NA,
                         KoHa=	NA,
                         RdRef=NA,
                         RdHa=	NA,
                         RdHd= NA,
                         RdS=NA,
                         O2=NA,
                         thetacj=NA,
                         thetaip=NA,
                         model.diff=NA,
                         gcw= NA,
                         model.gs=NA,
                         g0=NA,
                         g1=NA,
                         Patm=NA
                      ){
  VcmaxRef=abs(VcmaxRef)
  VqmaxRef=abs(VqmaxRef)
  param=list(R=R,TRef=TRef,VcmaxRef=VcmaxRef,VcmaxHa=VcmaxHa,VcmaxHd	=VcmaxHd,VcmaxS	=VcmaxS,
             VqmaxRef=VqmaxRef,VqmaxHa=	VqmaxHa,VqmaxHd=VqmaxHd, VqmaxS=VqmaxS,TpRef=TpRef,TpHa=TpHa,TpHd=TpHd,TpS=TpS,
             abso=abso,Beta=Beta,GstarRef=GstarRef, GstarHa	=GstarHa,etacRef=etacRef,etacHd=etacHd,etacS=etacS,etalRef=etalRef,etalHd=etalHd,etalS=etalS,phi1max=phi1max,
             KcRef=	KcRef,KcHa=	KcHa,KoRef=	KoRef,KoHa=KoHa,
             RdRef=RdRef,RdHa=RdHa,RdHd=RdHd,RdS=RdS,O2=O2,thetacj=thetacj,thetaip=thetaip,model.diff=model.diff,gcw= gcw,model.gs=model.gs,g0=g0,g1=g1,Patm=Patm)
  A_pred=f.ACi_JB(Ci=data$Ci,PFD=data$Qin,Tleaf=data$Tleaf,param=param)
  
  y<-dnorm(x=data$A,mean=A_pred$A,sd=(sigma),log=TRUE)
  return(-sum(y))
}

f.SumSq_JB<-function(Fixed,data,Start){
  if(!is.na(Start["VcmaxRef"])){Start["VcmaxRef"]=abs(Start["VcmaxRef"])}
  if(!is.na(Start["VqmaxRef"])){Start["VqmaxRef"]=abs(Start["VqmaxRef"])}
  param=c(Fixed,Start)
  y<-data$A-f.ACi_JB(Ci=data$Ci,PFD=data$Qin,Tleaf=data$Tleaf,param=param)$A
  return(sum(y^2))
}


f.plot_JB<-function(measures=NULL,list_legend,param,name='',type='ACi'){
  # Plot all data points
  if(type=='ACi'){x=measures$Ci
  xlab=expression(italic(C)[i]~ppm)}
  if(type%in%c('Aq','AQ')){x=measures$Qin
  xlab=expression(italic(Q)['in']~mu*mol~m^-2~s^-1)}
  if(!type%in%c('ACi','AQ','Aq')){print('type should be ACi or AQ')}
  plot(x=x,y=measures$A, main=name, xlab=xlab, ylab=expression(italic(A)~mu*mol~m^-2~s^-1),ylim=c(min(measures$A,na.rm = TRUE),1.15*max(measures$A,na.rm = TRUE)))
  if(!is.null(list_legend)){
    list_legend=list_legend[order(names(list_legend))]
    legend("bottomright",legend=mapply(FUN = function(x, i){paste(i,'=', round(x,2))}, list_legend, names(list_legend)),bty="n",cex=1)
  }
  legend("topleft",legend=c(expression(italic(A)[c]),expression(italic(A)[j]),expression(italic(A)[p]),expression(italic(A)),"Obs"),lty=c(2,2,2,1,0),
         pch=c(NA,NA,NA,NA,21),
         col=c("dark blue","dark red","dark green","dark grey","black"),bty="n",lwd=c(2,2,2,1,1),
         seg.len=2,cex=1,pt.cex=1)
  result=f.ACi_JB(Ci=measures$Ci,Tleaf=measures$Tleaf,PFD=measures$Qin,param=param)
  lines(x=x,y=result$A,col="dark grey",lwd=1)
  lines(x=x,y=result$Ac,lwd=2,col="dark blue",lty=2)
  lines(x=x,y=result$Aj,lwd=2,col="dark red",lty=2)
  lines(x=x,y=result$Ap,lwd=2,col="dark green",lty=2)
  box(lwd=1)
}


f.J_FvCB<-function(Q,alpha=0.85,Phi=0.425,Jmax=100,Theta=0.7){
  J=(Q*alpha*Phi+Jmax-sqrt((Q*alpha*Phi+Jmax)^2-4*Theta*Q*alpha*Phi*Jmax))/(2*Theta)
  return(J)
}

f.J_JB<-function(Q,Beta=0.5,phi1max=0.96,Vqmax=100,Ci=270,etal=0.75,etac=1,Gstar=42.75,alpha=0.85){
  alpha1=1-Beta
  eta=1-etal/etac+(3+7*Gstar/Ci)/((4+8*Gstar/Ci)*etac)
  J=Vqmax*Q/(Vqmax/(alpha*alpha1*phi1max)+Q)/eta ## Eq 30
  return(J)
}

f.Jmax_to_Vqmax<-function(Q=2000,alpha=0.85,Phi=0.425,Jmax=100,Theta=0.7,Beta=0.52,phi1max=0.96,Ci=800,etal=0.75,etac=1,Gstar=42.75){
  alpha1=1-Beta
  J=f.J_FvCB(Q=Q,alpha=alpha,Phi=Phi,Jmax=Jmax,Theta=Theta)
  eta_inv=1/(1-etal/etac+(3*Ci+7*Gstar)/((4*Ci+8*Gstar)*etac))
  Vqmax=Q*J/(Q*eta_inv-J/(alpha*alpha1*phi1max))
  return(Vqmax)
}

## Dependence of eta on Ci 
f.eta<-function(Ci,Gstar=42.75,etal=0.75,etac=1){
 eta=(1-etal/etac+(3*Ci+7*Gstar)/((4*Ci+8*Gstar)*etac))
  return(eta)
}


#param=f.make.param_JB()
#PFD=0:2000
#Cs=400
#Tleaf=25+273.16
#Tair=Tleaf
#RH=70
#physio=f.A_JB(PFD=0:2000,Cs = 400,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = param)
#plot(physio$A)
#abline(h=physio$Ac,col='blue')
#lines(physio$Aj,col='red')
#
#physio2=f.A_JB(PFD=2000,Cs = 10:2000,Tleaf = 25+273.16,Tair = 25+273.16,RH = 70,param = param)
#plot(physio2$A,ylim=c(-5,60))
#lines(physio2$Ac,col='blue')
#lines(physio2$Aj,col='red')

f.A_JB<-function(PFD,Cs,Tleaf,Tair,RH,param){
 
  Kc=f.arrhenius(PRef = param[['KcRef']],Ha = param[['KcHa']],Tleaf=Tleaf)
  Ko=f.arrhenius(PRef = param[['KoRef']],Ha = param[['KoHa']],Tleaf=Tleaf)
  Gstar=f.arrhenius(PRef = param[['GstarRef']],Ha = param[['GstarHa']],Tleaf = Tleaf)
  
  ## Mitochondrial respiration
  Rd=f.modified.arrhenius(PRef=param[['RdRef']],Ha = param[['RdHa']],Hd = param[['RdHd']],s=param[['RdS']],Tleaf = Tleaf)
  
  Vcmax=f.modified.arrhenius(PRef=param[['VcmaxRef']],Ha = param[['VcmaxHa']],Hd = param[['VcmaxHd']],s = param[['VcmaxS']],Tleaf = Tleaf)
  etal=f.desactivation(PRef=param[['etalRef']],Hd = param[['etalHd']],s = param[['etalS']],Tleaf=Tleaf)
  etac=f.desactivation(PRef=param[['etacRef']],Hd = param[['etacHd']],s = param[['etacS']],Tleaf=Tleaf)
  
  Vqmax=f.modified.arrhenius(PRef=param[['VqmaxRef']],Ha = param[['VqmaxHa']],Hd = param[['VqmaxHd']],s=param[['VqmaxS']],Tleaf = Tleaf)
  Tp=f.modified.arrhenius(PRef=param[['TpRef']],param[['TpHa']],param[['TpHd']],param[['TpS']],Tleaf)
  
  alpha1=1-param[['Beta']]
  
  JP700j=Vqmax*PFD/(Vqmax/(param[['abso']]*alpha1*param[['phi1max']])+PFD) ## Eq 30
  
  wi=0.61365*exp(17.502*(Tleaf-273.16)/(240.97+(Tleaf-273.16)))/(param[['Patm']])
  ws=0.61365*exp(17.502*(Tair-273.16)/(240.97+(Tair-273.16)))/(param[['Patm']])*RH/100
  
  ### VPDleaf for USO model in kPa
  ds=(wi-ws)*param[['Patm']]
  
  if(param[['model.diff']]==0){k=0;l=0;gcw=0;q=param[['g0']]}
  if(param[['model.diff']]==1){k=(wi-ws)/(2-(wi+ws));l=0;gcw=0;q=param[['g0']]} 
  if(param[['model.diff']]==2){k=(wi-ws)/(2-(wi+ws));l=gcw/20;q=param[['g0']]-gcw}
  
  ## Solution for the coupled equations representing the photosynthesis, the stomatal conductance and the co2 diffusion from the surface of the leaf to the substomatal cavity
  # A solution exists for A in the form A=((Ci-Gstar)x)/((Ci+y))-Rd. For the JB model Aj limitation we have: 
  x=(Vqmax*PFD*etac)/((Vqmax/(alpha1*param[['phi1max']])+PFD)*(4*etac-4*etal+3))
  y=(8*Gstar*etac-8*Gstar*etal+7*Gstar)/(4*etac-4*etal+3)
  
  Cij=f.solv(x=x,y=y,Cs=Cs,Rd=Rd,Gstar=Gstar,q=q,g1=param[['g1']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
  
  ## With the USO model, gsw can be negative under the light compensation point. Below, we use the constraint gsw>=g0 to keep positive gsw
  if(param[["model.gs"]]==0){g1=-1}else(g1=0)
  Cij_g0=f.solv(x=x,y=2*Gstar,Cs=Cs,Rd=Rd,Gstar=Gstar,q=q,g1=g1,ds=ds,RH=RH,k=k,l=l,model=param[["model.gs"]])
  Cij=pmax(Cij,Cij_g0)
  Cij[(Cij-Cs)>=0]=Cij_g0[(Cij-Cs)>=0]
  
  eta=1-etal/etac+(3+7*Gstar/Cij)/((4+8*Gstar/Cij)*etac) ## Eq 15 c
  JP680j=JP700j/eta ## Eq 30
  Aj=JP680j/(4+8*Gstar/Cij)*(1-Gstar/Cij) ## Eq 31 without Rd
  
  Cic=f.solv(x=Vcmax,y=Kc*(1+param[['O2']]/Ko),Cs=Cs,Rd=Rd,Gstar=Gstar,q=q,g1=param[['g1']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
  Ac=Vcmax*Cic/(Cic+Kc*(1+param[['O2']]/Ko))*(1-Gstar/Cic) ## Eq 32 without Rd
  
  Ap=3*Tp
  Cip=f.solv(x=3*Tp,y=-Gstar,Cs=Cs,Rd=Rd,Gstar=Gstar,q=q-gcw,g1=param[['g1']],ds=ds,RH=RH,k=k,l=l,model=param[['model.gs']])
  
  
  Ci=Cij
  if(!is.null(which(Ac<Aj))&length(Cic)==length(Cij)){Ci[which(Ac<Aj)]=Cic[which(Ac<Aj)]}
  if(!is.null(which(Ac<Aj))&length(Cic)!=length(Cij)){Ci[which(Ac<Aj)]=Cic}
  A=pmin(Ac,Aj)
  if(!is.null(which(Ap<A))){Ci[which(Ap<A)]=Cip[which(Ap<A)]}
  Ai=f.smooth(A1 = Ac,A2 = Aj,theta=param[['thetacj']])
  Ag=f.smooth(A1=Ai,A2=Ap,theta=param[['thetaip']])
  An=Ag-Rd
  
  gs=f.gs(A=An,Cs=Cs,ds=ds*1000,Rd=Rd,Gstar=Gstar,RH=RH,g0=param[['g0']],g1=param[['g1']],model =param[['model.gs']])
    
    Ec=gcw*(wi-ws)
    Es=gs*(wi-ws)/(1-(wi+ws)/2)
    ET=Ec+Es
    output=list(A=An,Ac=Ac,Aj=Aj,Ap=Ap,Ag=Ag,gs=gs,Ci=Ci,ds=ds,ET=ET,Ec=Ec,Es=Es)
    return(output)
}


## Apply only for USO models and BWB model of conductance
f.solv<-function(x,y,Cs,Rd,Gstar,g1,ds,RH,k,l,q,model){
  if(model=="USO"|model==0){
    m=1.6*(1+g1/(ds)^0.5)*1/Cs
    if(g1==-1){m=0} ## I use this trick to calculate Ci at g0, see f.A function.
  }else if(model=="USO_simpl"|model==1){
    m=1.6*(g1/(ds)^0.5)*1/Cs
  }else if(model=="BWB"|model==2){
    m=(g1*RH/100)*1/Cs
  }else{print('Please, use BWB or USO or USO_simpl conductance model')}
  a=1/1.6*Rd*m+Rd*k*m-1/1.6*m*x-k*m*x-1/1.6*q-k*q-l 
  b=1/1.6*Gstar*m*x+Gstar*k*m*x-1/1.6*Rd*Cs*m+Rd*Cs*k*m+1/1.6*Rd*m*y+ Rd*k*m*y+1/1.6*Cs*m*x-Cs*k*m*x+1/1.6*Cs*q-Cs*k*q-1/1.6*q*y-k*q*y+Cs*l- l*y+ Rd-x
  c=-1/1.6*Gstar*Cs*m*x+Gstar*Cs*k*m*x-1/1.6*Rd*Cs*m*y+Rd*Cs*k*m*y+1/1.6*Cs*q*y-Cs*k*q*y+Cs*l*y+x*Gstar+Rd*y  
  
  Ci2=(-b-(b^2-4*a*c)^0.5)/(2*a)
  Ci1=(-b+(b^2-4*a*c)^0.5)/(2*a)
  return(pmax(Ci1,Ci2)) 
}

#' @title Leaf to air water vapour pressure deficit calculation
#' @description This function calculates the leaf water pressure deficit (VPDleaf or Ds) using the temperature of the leaf, the temperature of the air and its relative humidity. The saturation vapor pressure of water is approximated usung the Tetens 1930 equation. 
#' @param Tleaf Temperature of the leaf in Kelvin
#' @param Tair Temperature of the air in Kelvin
#' @param RH Humidity of the air (0 to 100)
#' @references Tetens O. 1930. Uber einige meteorologische Begriffe. Z. geophys 6: 297-309.
#' @return Ds in Pascal
#' @export
#' @examples f.ds(Tleaf=273.16 + 30, Tair=273.16+28, RH=70)
f.ds<-function(Tleaf,Tair,RH){
  ds=(0.6108*exp(17.27*(Tleaf-273.16)/(Tleaf-273.16+237.3))-0.6108*exp(17.27*(Tair-273.16)/(Tair-273.16+237.3))*RH/100)*1000
  return(ds)
}


#' @title Conductance model for stomatal conductance to water vapour
#' @description Semi-empirical model of the leaf conductance to water vapour
#' @param A Net assimilation in micromol.m-2.s-1, i-e, the assimilation in presence of respiration
#' @param cs CO2 at the surface of the leaf in ppm
#' @param ds Leaf surface to air vapour pressure deficit in Pa
#' @param Rd Dark respiration for the non linear model. Be carefull that the net assimilation is implemented in all the models except the non linear where the gross assimilation is used (A + Rd)
#' @param RH Humidity at the surface of the leaf (0 - 100). ds or RH as to be specified
#' @param g0 Constant of the USO model, representing the conductance when A is 0, in mol.m-2.s-1
#' @param g1 Slope parameter, between 1.14 and 3.58 KPa^0.5 (Wu et al., 2019)
#' @param model Stomatal model ("USO", "USO_simpl" or "BWB")
#' @export
#' @details USO : #gs=g0+1.6*(1+g1/(ds/1000)^0.5)*(A)/cs
#' USO_simpl : gs=g0+1.6*(g1/(ds/1000)^0.5)*(A)/cs
#' BWB : gs=g0+g1*(A*RH/100)/cs
#' @return This function returns the optimal stomatal conductance to water vapour in mol.m-2.s-1
#' @references Medlyn, B.E., Duursma, R.A., Eamus, D., Ellsworth, D.S., Colin Prentice, I., Barton, C.V.M., Crous, K.Y., de Angelis, P., Freeman, M. and Wingate, L. (2012), Reconciling the optimal and empirical approaches to modelling stomatal conductance. Glob Change Biol, 18: 3476-3476. doi:10.1111/j.1365-2486.2012.02790.x
#'  Wu, J, Serbin, SP, Ely, KS, et al. The response of stomatal conductance to seasonal drought in tropical forests. Glob Change Biol. 2020; 26: 823– 839. https://doi.org/10.1111/gcb.14820
#'  Leuning, R., Kelliher, F. M., De Pury, D. G. G., & Schulze, E. D. (1995). Leaf nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to canopies. Plant, Cell & Environment, 18(10), 1183-1200.
#'  Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In Progress in photosynthesis research (pp. 221-224). Springer, Dordrecht.
#' @examples gs=f.gs(A=30,cs=400,ds=1500,g0=0.01,g1=2)
f.gs<-function(A,Cs,ds=NULL,RH=NULL,Rd=NULL,Gstar=NULL,g0,g1,model="USO"){
  if(model=="USO"|model==0){
    gs=g0+1.6*(1+g1/(ds/1000)^0.5)*(A)/Cs
  } else if(model=="USO_simpl"|model==1){
    gs=g0+1.6*(g1/(ds/1000)^0.5)*(A)/Cs
  } else if(model=="BWB"|model==2){
    gs=g0+g1*(A*RH/100)/Cs
  }
  else{print(paste("Parameter model =",model,"is not in the list of implemented models"))}
  return(gs)
}
