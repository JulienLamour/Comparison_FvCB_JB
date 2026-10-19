## This code is used to do the step 1 of the data curation process.
## It is used to analyse and check the data quality. Bad points are removed as well as bad curves.

# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Set the working directory to the 'Schmiege_et_al_2021' folder where the data is located
setwd(file.path(path,'/Datasets/A-Ci/Schmiege_et_al_2021'))

# Load the ACi data that were processed in step 0
load('0_curated_data.Rdata',verbose=TRUE)

# Adding a Quality Check column with two values: "ok" or "bad", where "bad" correspond to bad data 
# that is not used to fit the ACi curves in the next step.
curated_data$QC='ok'

# Displaying Ci values to spot bad data
hist(curated_data$Ci) ## No below zero value
hist(curated_data$gsw) ## No below zero value
hist(curated_data$Tleaf)
hist(curated_data$Qin)

# I create a "Quality check" table where I manually inform bad points. 
# When I start the quality check the QC table is empty: QC_table=cbind.data.frame(SampleID_num=c(),Record=c()) 
# I plot the data in the next step (line starting with pdf below). 
# I then manually check the plots and write down the SampleID_num and the Record corresponding to bad points
# This works well if there are few points to remove.
# You can use another method if you need.

QC_table=cbind.data.frame(SampleID_num=c(3,4,6,7,10,12,16,24,27,29,53,53,57,57,59),
                          Record=c(2,2,2,5,10,3,6,15,18,4,10,15,2,16,3)) 
# Here I flag the bad curves by writing down the SampleID_num of the bad curves.
ls_bad_curve=c(1,40:43,58,60,61)
ls_bad_curve=c(ls_bad_curve,8,21,32,50,59) ## I removed these curves with weird behavior during fitting. I suppose the light was to low to be at saturation? It appears that Aj was not limiting in these curves
# Here I flag all the bad curves and bad points
curated_data[paste(curated_data$SampleID_num,curated_data$Record)%in%paste(QC_table$SampleID_num,QC_table$Record),'QC']='bad'
curated_data[curated_data$SampleID_num%in%ls_bad_curve,'QC']='bad'

# Here I remove the transient CO2 values that are used to wait for Rubisco activation
##Function to find and Remove duplicated values
Remove_ACi_duplicate<-function(ACi_curve){
  ##Rounding the Ci values
  ACi_curve=round(ACi_curve,-1)
  Remove=rep('NO',length(ACi_curve))
  ## We Remove the last point of the curve if it is a duplicate
  if(length(ACi_curve)%in%which(duplicated(ACi_curve))){Remove[length(ACi_curve)]='YES'}
  ## We Remove all the duplicates except the last one, so theoretically, if there are several 400 CO2 levels, only the last one is kept
  Remove[which(duplicated(ACi_curve[1:(length(ACi_curve))],fromLast = TRUE))]='YES'
  Remove[1:5]='NO'
  Remove[which(duplicated(ACi_curve[1:5],fromLast = TRUE))]='YES'
  return(Remove)
}

Duplicated_CO2=ave(curated_data$CO2s,curated_data$SampleID_num,FUN=Remove_ACi_duplicate)
curated_data[Duplicated_CO2=="YES","QC"]="bad"

# Creating a PDF file where all the ACi curves are plot. The points in red are bad points or curves and the points in black are ok.
pdf(file='1_QA_QC_ACi.pdf',)
for(SampleID_num in sort(unique(curated_data$SampleID_num))){
plot(x=curated_data[curated_data$SampleID_num==SampleID_num,'Ci'],y=curated_data[curated_data$SampleID_num==SampleID_num,'A'],main=unique(curated_data[curated_data$SampleID_num==SampleID_num,'SampleID_num']),cex=2,xlab='Ci',ylab='A')
  text(x=curated_data[curated_data$SampleID_num==SampleID_num,'Ci'],y=curated_data[curated_data$SampleID_num==SampleID_num,'A'], labels=curated_data[curated_data$SampleID_num==SampleID_num,'Record'],cex=0.7)
if(any(curated_data[curated_data$SampleID_num==SampleID_num,'QC']=='bad')){
  points(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Ci'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','A'],cex=2,col='red')
  text(x=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Ci'],y=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','A'], labels=curated_data[curated_data$SampleID_num==SampleID_num&curated_data$QC=='bad','Record'],cex=0.7,col='red')
}
}
dev.off()

# Keeping only the good points and curves
curated_data=curated_data[curated_data$QC=='ok',]

# Saving the Quality checked data to be used in the next step.
save(curated_data,file='1_QC_ACi_data.Rdata')