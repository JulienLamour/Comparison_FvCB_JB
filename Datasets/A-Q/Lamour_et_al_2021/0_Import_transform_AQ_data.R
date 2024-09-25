## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from .....
## The dataset is available publicly :


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Lamour_et_al_2021' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Lamour_et_al_2021'))

# Import the authors original raw data
original_data=read.csv('PA-SLZ_2020_Aq_data.csv')
SampleInfo=read.csv(file = "PA-SLZ_2020_SampleDetails.csv")
SampleInfo$Dataset="Lamour_et_al_2021"
SampleInfo$Species=SampleInfo$Species_Name
SampleInfo$Species_type="Wild"
SampleInfo$Biome="Tropical"
original_data=merge(original_data,SampleInfo,by="SampleID")
original_data[original_data$Vertical_Elevation>5,"SunShade"]="Shade"
original_data[original_data$Vertical_Elevation<3,"SunShade"]="Sun"
original_data$PFT="Tropical"

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record and gsw columns were missing (column 2 and 7, respectively).
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","Ci","CO2s","CO2r","SampleID","Patm","Qin","RHs","Tleaf","Remove","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Remove","Species","Species_type","Biome","SunShade","Dataset","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$gsw = NA
curated_data$Record = NA

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))
curated_data=curated_data[curated_data$Remove=="NO",]

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
