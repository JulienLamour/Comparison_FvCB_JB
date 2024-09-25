## This code is used to do the step 1 of the data curation process.
## It is used to analyse and check the data quality. Bad points are removed as well as bad curves.

# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(spectratrait)

# Find the path of the top relative directory
path=here()

# Set the working directory to the 'Rogers_et_al_2019' folder where the data is located
setwd(file.path(path,'/Datasets/A-Ci/Rogers_et_al_2019'))

# Load the ACi data that were processed in step 0
load('0_curated_data.Rdata',verbose=TRUE)


# I only keep the curves with a spectra associated and that constitutes the main dataset

ecosis_id <- "bf41fff2-8571-4f34-bd7d-a3240a8f7dc8"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)

curated_data=curated_data[curated_data$SampleID%in%dat_raw$Sample_ID,]

# Adding a Quality Check column with two values: "ok" or "bad", where "bad" correspond to bad data 
# that is not used to fit the ACi curves in the next step.
curated_data$QC='ok'

# Displaying Ci values to spot bad data
curated_data$QC='ok'

# The authors already flagged bad data, so I re-use their classification here and complete it.
curated_data[!is.na(curated_data$QCauthors)&curated_data$QCauthors==1,'QC']='bad'
hist(curated_data$Ci) ##
curated_data[curated_data$Ci<0,'QC']='bad'
hist(curated_data$gsw) ## No below zero value
curated_data[curated_data$gsw<0,'QC']='bad'

# I create a "Quality check" table where I manually inform bad points. 
# When I start the quality check the QC table is empty: QC_table=cbind.data.frame(SampleID_num=c(),Record=c()) 
# I plot the data in the next step (line starting with pdf below). 
# I then manually check the plots and write down the SampleID_num and the Record corresponding to bad points
# This works well if there are few points to remove.
# You can use another method if you need.


QC_table=cbind.data.frame(SampleID_num=c(46,46,46,53,70,71,72,73,110,110,111,157),
                          Record=c(16,6,7,1,1,16,1,16,16,17,1,18)) ## 
ls_bad_curve=c(43,45,59,74,136,139,147,148,149,150,157,158,159,160,161,162,176,177,178)

# Here I flag all the bad curves and bad points
curated_data[paste(curated_data$SampleID_num,curated_data$Record)%in%paste(QC_table$SampleID_num,QC_table$Record),'QC']='bad'
curated_data[curated_data$SampleID_num%in%ls_bad_curve,'QC']='bad'

############ Here I manually consider some points 'good' so the estimation of Jmax is correct
curated_data[curated_data$SampleID_num==52&curated_data$Record>12,'QC']='ok'
curated_data[curated_data$SampleID_num==72&curated_data$Record>11,'QC']='ok'
curated_data[curated_data$SampleID_num==73&curated_data$Record>12,'QC']='ok'
curated_data[curated_data$SampleID_num==78&curated_data$Record>12,'QC']='ok'
curated_data[curated_data$SampleID_num==79&curated_data$Record>13,'QC']='ok'
curated_data[curated_data$SampleID_num==103&curated_data$Record>13,'QC']='ok'
curated_data[curated_data$SampleID_num==104&curated_data$Record>13,'QC']='ok'
curated_data[curated_data$SampleID_num==115&curated_data$Record>11,'QC']='ok'
curated_data[curated_data$SampleID_num==125&curated_data$Record>12,'QC']='ok'
curated_data[curated_data$SampleID_num==127&curated_data$Record>12,'QC']='ok'


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
curated_data=curated_data[curated_data$QC=="ok",]

# There are some duplicate curves in the same leaves that I remove:
dupp=tapply(curated_data$Replicate,INDEX = curated_data$SampleID,FUN = function(x)length(unique(x)))
dupp=dupp[dupp>1]
curated_data=curated_data[-which(curated_data$SampleID%in%names(dupp)&curated_data$Replicate!=3),]

# Saving the Quality checked data to be used in the next step.
save(curated_data,file='1_QC_ACi_data.Rdata')
