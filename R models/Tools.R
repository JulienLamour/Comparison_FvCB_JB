###########################################################################
### This code includes miscellaneous functions for gas exchange studies ###
###########################################################################


#########################################################################################
###  Standard variable names for gas exchange intruments, based on Ely et al. 2021    ###
###  https://doi.org/10.1016/j.ecoinf.2021.101232.                                    ###
#########################################################################################
LICOR6400_column=c('','Obs','Photo','Ci','CO2S','CO2R','Cond','Press','PARi','RH_S','Tleaf')
ESS_column=c('SampleID','Record','A','Ci','CO2s','CO2r','gsw','Patm','Qin','RHs','Tleaf')

LICOR6400toESS=as.data.frame(t(ESS_column))
colnames(LICOR6400toESS)=LICOR6400_column


###################################################################################
###   Function to plot AQ curves                                                ###
###   This function is used to do the quality analysis of all the AQ curves     ###
###   It displays the SampleID_num of the curve in the title of the plot        ###
###   In the left y axis we display the photosynthesis rate                     ###
###   while in the right y axis we display the intercellular CO2                ###
###   concentration. Bad values flagged in the input dataframe                  ###
###   are displayed in red (photosynthesis points) or red (Ci points)           ###
###################################################################################
#' @param PDFfilename Filename of the pdf for saving all the individual curves
#' @param curated_data Data frame with as columns A (photosynthesis rate), Qin (incident irradiance), Ci (intercellular CO2 concentraiton), QC (quality check variable with values "ok" or "bad"), and SampleID_num (the unique identifier of the curve)
#'
#' @return The function does not return anything but create a pdf file.
#' @export
#'
#' @examples
f.QaQc_AQ_curves<-function(PDFfilename = '1_QA_QC_AQin.pdf', curated_data = NA){
  pdf(file = PDFfilename)
  par(mar = c(5, 4, 4, 4) + 0.3)
  for(SampleID_num in sort(unique(curated_data$SampleID_num))){
    plot(x=curated_data[curated_data$SampleID_num==SampleID_num,'Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num,'A'],main=unique(curated_data[curated_data$SampleID_num==SampleID_num,'SampleID_num']),cex=2,xlab='Qin',ylab='A')
    text(x=curated_data[curated_data$SampleID_num==SampleID_num,'Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num,'A'], labels=curated_data[curated_data$SampleID_num==SampleID_num,'Record'],cex=0.7)
    if(any(curated_data[curated_data$SampleID_num==SampleID_num,'QC']=='bad')){
      points(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','A'],cex=2,col='red')
      text(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','A'], labels=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Record'],cex=0.7,col='red')
    }
    par(new = TRUE)
    plot(x=curated_data[curated_data$SampleID_num==SampleID_num,'Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num,'Ci'],xlab='',ylab='',axes=FALSE,col="blue",cex=1.5)
    axis(side=4, at = pretty(range(curated_data[curated_data$SampleID_num==SampleID_num,'Ci'])),col="blue",col.axis="blue")
    mtext("Ci", side=4, line=3,col="blue")
    text(x=curated_data[curated_data$SampleID_num==SampleID_num,'Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num,'Ci'], labels=curated_data[curated_data$SampleID_num==SampleID_num,'Record'],cex=0.6,col="blue")
    if(any(curated_data[curated_data$SampleID_num==SampleID_num,'QC']=='bad')){
      points(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Ci'],cex=1.5,col= "pink")
      text(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Qin'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Ci'], labels=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Record'],cex=0.6,col='pink')
    }
  }
  dev.off()
}

###################################################################################
###   This function imports the gas exchange data from LICOR 6800 excel files   ###
###################################################################################
#' @title Import Licor 6800 file
#' @description This functions allows to import the excel files produced by LICOR as a data.frame.
#' IMPORTANT: The excel files must be opened and saved before using this function (the Excel calculations are not done until the file is open, so the calculated colums will show 0s if not saved before being imported)
#' @param nskip_header Number of lines to skip in the Excel files to find the column names
#' @param nskip_data Number of lines to skip in the Excel files to find the data
#' @param do.print Print the 5 top lines of the file?
#' @param file File path
#' @param column_display Column you want to display after the import to verify if it worked correctly
#'
#' @examples
f.import_licor6800<-function(nskip_header=16,nskip_data=18,do.print=TRUE,file,column_display=c('A','gsw','Qin','Ci','Species','Canopy','Pheno_Age','Barcode','file')){
  if(do.print){print(file)}
  header=make.names(as.data.frame(readxl::read_excel(path = file,skip = nskip_header,n_max = 1,.name_repair = 'minimal',col_names = FALSE)))##'minimal' to speed up the import
  data_6800=as.data.frame(readxl::read_excel(path = file,skip = nskip_data,col_names = header,.name_repair = 'minimal'))
  data_6800[,'date']=data_6800[1,'date']
  data_6800=cbind(data_6800,file=rep(file,nrow(data_6800)))
  if(do.print){print(head(data_6800[,column_display]))}
  if(all(data_6800$A==0)){print(paste('Open and save the file first',file))}
  return(data_6800)
}


