## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from  https://doi.org/10.1111/pce.14448
## The dataset was sent by S. Schmiege


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Schmiege_et_al_2022' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Schmiege_et_al_2022'))

# Import the authors original raw data
original_data=read.csv('2024_05_01_rawAQcurves_adult_trees.csv')
original_data$Dataset="Schmiege_et_al_2022"
original_data$Species_type="Wild"
original_data[original_data$location=="AK","Biome"]="Arctic"
original_data[original_data$location=="BRF","Biome"]="Temperate"
original_data$SunShade="Shade"
original_data[original_data$canopy=="high","SunShade"]="Sun"
original_data$Species="Picea glauca"
original_data$PFT="Gymnosperm"

original_data$SampleID=original_data$concat
for(leaf in unique(original_data$SampleID)){
  original_data[original_data$SampleID==leaf,"Record"]=1:nrow( original_data[original_data$SampleID==leaf,])
}

# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "Record","A","Ci","CO2_s","CO2_r","gsw","gsw","Qin","RHcham","TleafEB","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")
curated_data$Patm=NA


# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
