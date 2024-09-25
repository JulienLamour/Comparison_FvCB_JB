## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Zhou et al 2023
## The dataset was sent by Zhou

# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Zhou_et_al_2023' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Zhou_et_al_2023'))

# Import the authors original raw data
original_data=read.csv('Zhou_et_al_2023_AQ_data.csv')
original_data$Dataset="Zhou_et_al_2023"
original_data$Species_type="Agricultural"
original_data$Species="Oryza sativa L."
original_data$Biome="Tropical"
original_data$SunShade="Sun"
original_data$PFT="Crop"

## I only keep the control genotypes that were not modified
original_data=original_data[original_data$Genotype=="C",]

## And I only keep the data recorded at ambient O2 (21%)
original_data=original_data[original_data$O2==21,]

##Creation of a SampleID column
original_data$SampleID=paste(original_data$Background,original_data$Stage,original_data$Replicate)

## Creation of a record numer
original_data$record=NA
for (SampleID in unique(original_data$SampleID)){
  original_data[original_data$SampleID==SampleID,"record"] = 1 : nrow (original_data[original_data$SampleID==SampleID,])
}

## Only keeping the light curves part of the data
original_data=original_data[original_data$record<11,]
original_data=original_data[-which(original_data$record==10&original_data$Iinc>100),]
# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "record","A","Ci","Ca","Ca","gs","gs","Iinc","Iinc","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")

# Putting NA in the columns without information
curated_data$RHs=NA
curated_data$CO2r=NA
curated_data$Patm=NA

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
