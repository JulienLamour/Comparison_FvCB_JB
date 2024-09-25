## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Sun et al. 2012
## The dataset was sent by Ülo Niinemets


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Sun_et_al_2012' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Sun_et_al_2012'))

# Import the authors original raw data
original_data=read.csv('Populus-light-response-curves-GCB-paper.csv')
#SampleInfo=read.csv(file = "")
original_data$Dataset="Sun_et_al_2012"
original_data$Species="Populus tremula x P. tremuloides"
original_data$Species_type="Agricultural"
original_data$Biome="Temperate"
original_data$SunShade="NA"
original_data$PFT="Crop"
original_data$SampleID=paste(original_data$Measurement.date,original_data$Measurement.series,original_data$Clone,original_data$growth.CO2,original_data$leaf.code,original_data$Ca)
# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record and gsw columns were missing (column 2 and 7, respectively).
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","Ci","Ca","Ca","SampleID","SampleID","PFD","SampleID","Tleaf","SampleID","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Remove","Species","Species_type","Biome","SunShade","Dataset","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$gsw = NA
curated_data$Record = NA
curated_data$CO2s = curated_data$RHs = curated_data$Patm = curated_data$Remove = NA



curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))


# Removing dupplicates values

# Here I remove the duplicated light values
##Function to find and Remove duplicated values
Remove_AQ_duplicate<-function(AQ_curve){
  ##Rounding the ci values
  AQ_curve=round(AQ_curve,-1)
  Remove=rep('NO',length(AQ_curve))
  ## We Remove all the duplicates except the last one
  Remove[which(duplicated(AQ_curve[1:(length(AQ_curve))],fromLast = TRUE))]='YES'
  return(Remove)
}
curated_data$QC="ok"
Duplicated_Q=ave(curated_data$Qin,curated_data$SampleID_num,FUN=Remove_AQ_duplicate)
curated_data[Duplicated_Q=="YES","QC"]="bad"

curated_data=curated_data[curated_data$QC=="ok",]
# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
