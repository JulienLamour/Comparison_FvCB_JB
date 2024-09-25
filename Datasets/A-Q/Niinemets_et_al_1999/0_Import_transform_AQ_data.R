## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Niinemets et al. 1999 (data sent by the author)


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Niinemets_et_al_1999' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Niinemets_et_al_1999'))

# Import the authors original raw data
original_data=read.csv('GCB-1999-Light_curves.csv')
original_data$Species="Quercus suber"
original_data$Dataset="Niinemets_et_al_1999"
original_data$Species_type="Wild"
original_data$Biome="Greenhouse"
original_data$SunShade="Sun"
original_data$PFT="Mediterranean"


# Select specific columns from the original data to create the curated data set
print(ESS_column)
curated_data=original_data[,c("SampleID", "SampleID","A","Ci","Ca","Ca","gsw","gsw","Qin","Qin","Tleaf","ALVPD","Species","Species_type","Biome","SunShade","PFT")]

# Rename the columns of the curated dataset with the ESS standard
colnames(curated_data)=c(ESS_column,"VPDleaf","Species","Species_type","Biome","SunShade","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$Record = NA
curated_data$Patm = NA
curated_data$RHs = NA
curated_data$VPDleaf = curated_data$VPDleaf*0.1 
curated_data$gsw = curated_data$gsw/1000
curated_data$Dataset="Niinemets_et_al_1999"

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
