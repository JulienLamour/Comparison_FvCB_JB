## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Mao et al. 2026.
## 10.48130/forres-0026-0023 (sent by the authors)


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Mao_et_al_2026' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Mao_et_al_2026'))

# Import the authors original raw data
original_data=read.csv('AQ_data.csv')
original_data$SampleID = paste(original_data$Species,original_data$Replicate, original_data$type)
original_data$Dataset="Mao_et_al_2026"
original_data$Species_type="Ornamental"
original_data$Biome="Subtropical"
original_data$SunShade="Sun"
original_data$PFT="Subtropical"

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record and Patm columns were missing
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","ci","ca","CO2abs","GH2O","SampleID","PARtop","rh","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")

curated_data$gsw = curated_data$gsw/1000

# The columns gsw and Record were still filled with SampleID info. I replace by NA values

curated_data$Record = NA
curated_data$Patm = NA

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')


for(SampleID in unique(curated_data$SampleID)){
  plot(x = curated_data[curated_data$SampleID == SampleID,"Qin"], y = curated_data[curated_data$SampleID == SampleID,"A"], main = SampleID)
}
