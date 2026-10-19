## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the ACi data from Davidson et al 2024.
## The dataset is available publicly : Davidson K ; Ely K ; Anderson J ; Rogers A ; Serbin S (2024): Photosynthetic and stomatal response, and leaf traits, of two pine species, Talladega National Forest, Alabama, 2023. Crossing the boundary: characterization of vegetation properties to improve the coupled model representation of the land-atmosphere interface, ESS-DIVE repository. Dataset. doi:10.15485/2332922 


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(hms)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names (variable ESS_column) from the "Tools.R" file
source(file.path(path,'/R models/Tools.R')) 

# Set the working directory to the 'Davidson_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Ci/Davidson_et_al_2024'))

# Import the authors original raw data
original_data=read.csv('AL_2023_ACi_AllData.csv')

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record and CO2r columns were missing.
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","Ci","CO2s","SampleID","gsw","Patm","Qin","RHs","Tleaf","POSTIX_Date")]
curated_data$Dataset="Davidson_et_al_2024"
# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Date_time","Dataset")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$Record = NA
curated_data$CO2r = NA
curated_data$time=as_hms(as.POSIXct(x = curated_data$Date_time))
curated_data=curated_data[order(curated_data$SampleID,curated_data$time),]

# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
