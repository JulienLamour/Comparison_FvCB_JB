## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Niinemets et al. 2015: Niinemets, Ül., Sun, Z., and Talts, E. (2015) Controls of the quantum yield and saturation light of isoprene emission in different-aged aspen leaves. Plant Cell Environ, 38: 2707–2720. doi: 10.1111/pce.12582.
## The dataset was sent by Ülo Niinemets


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Niinemets_et_al_2015' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Niinemets_et_al_2015'))

# Import the authors original raw data
original_data=read.csv('Populus-light-response-curves-PCE-paper.csv')
original_data$Dataset="Niinemets_et_al_2015"
original_data$Species="Populus tremula x P. tremuloides"
original_data$Species_type="Agricultural"
original_data$Biome="Temperate"
original_data$SunShade="NA"
original_data$PFT="Crop"
original_data$SampleID=paste(original_data$Measurement.date,original_data$Measurement.series,original_data$Clone,original_data$growth.CO2,original_data$leaf.code)
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

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
