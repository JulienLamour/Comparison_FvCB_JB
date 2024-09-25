## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from .....
## The dataset is available publicly :


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Lamour_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Lamour_et_al_2024'))

# Import the authors original raw data
original_data=read.csv('Brazil_2023_Aq_data.csv')
SampleInfo=read.csv('SampleInfo.csv')
TreeInfo=read.csv('TreeInfo.csv')
SampleInfo=merge(x=SampleInfo,y=TreeInfo,by="TreeID",all.x=TRUE)
original_data=merge(x=original_data,y=SampleInfo,by="SampleID",all.x=TRUE)
original_data$Dataset="Lamour_et_al_2024"
original_data$Species_type="Wild"
original_data$Biome="Tropical"
original_data$SunShade="Sun"
original_data$PFT="Tropical"

## I only keep the "good" data flagged by the authors

original_data=original_data[original_data$Remove=="NO",]
# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "record","A","Ci","CO2s","CO2r","gsw","Patm","Qin","RHs","Tleaf","Remove","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Remove","Species","Species_type","Biome","SunShade","Dataset","PFT")


# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))
curated_data$gsw=curated_data$gsw/1000
# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
