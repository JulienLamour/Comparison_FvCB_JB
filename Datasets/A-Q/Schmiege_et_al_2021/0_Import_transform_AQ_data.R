## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from https://academic.oup.com/treephys/article/41/2/223/5911517
## The dataset is available publicly : https://doi.org/10.5061/dryad.1g1jwstss 


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Schmiege_et_al_2021' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Schmiege_et_al_2021'))

# Import the authors original raw data
original_data=read.csv('SSchmiege_Table3_LightResponse_curves_raw.csv')
original_data$Dataset="Schmiege_et_al_2021"
original_data$Species_type="Wild"
original_data$Biome="Tropical"
original_data$SunShade="Sun"
original_data$PFT="Gymnosperm"

Species=cbind.data.frame(spp=c("pike","pida","pikr","dael","daim","nawa","pone"),
                         Species=c("Pinus kesiya","Pinus dalatensis","Pinus krempfii","Dacrydium elatum","Dacrycarpus imbricatus","Nageia wallichiana","Podocarpus neriifolius")) 
original_data=merge(original_data,Species,by="spp",all.x=TRUE)
original_data$SampleID=paste(original_data$spp,original_data$Date.actual,original_data$ID,original_data$Rep)
for(leaf in unique(original_data$SampleID)){
  original_data[original_data$SampleID==leaf,"record"]=1:nrow( original_data[original_data$SampleID==leaf,])
}

# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "record","A","Ci","CO2_s","CO2_r","gsw","gsw","Qin","H2O_s","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")
curated_data$RHs=NA
curated_data$Patm=NA


# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
