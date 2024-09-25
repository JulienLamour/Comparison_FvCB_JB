## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the ACi data from Rogers et al 2017.
## The dataset is available publicly with DOI 10.5440/1336809

# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names (variable ESS_column) from the "Tools.R" file
source(file.path(path,'/R models/Tools.R')) 

# Set the working directory to the 'Rogers_et_al_2019' folder where the data is located
setwd(file.path(path,'/Datasets/A-Ci/Rogers_et_al_2019'))

# Import the authors original raw data
original_data=read.csv('NGEE-Arctic_A-Ci_curves_2012-2015.csv',skip = 8)

# Select specific columns from the original data to create the curated data set
curated_data=original_data[,c("Sample_Barcode", "Sample_Barcode","Photo","Ci","CO2S","CO2R","Cond","Press","PARi","RH_S","Tleaf","QC","Replicate")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.

colnames(curated_data)=c(ESS_column,'QCauthors','Replicate')

curated_data[curated_data$Replicate==-9999,'Replicate']=1
curated_data$Dataset="Rogers_et_al_2019"

# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(paste(curated_data$SampleID,curated_data$Replicate)))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')

