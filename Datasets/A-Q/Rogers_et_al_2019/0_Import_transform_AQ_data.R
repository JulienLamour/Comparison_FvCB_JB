## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Rogers et al 2017 
## The dataset is available publicly :10.5440/1482338


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Rogers_et_al_2019' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Rogers_et_al_2019'))

# Import the authors original raw data
original_data=read.csv('NGA175_AQcurves_Barrow_2016.csv')

original_data$Dataset="Rogers_et_al_2019"
original_data$Species_type="Wild"
original_data$Biome="Arctic"
original_data$SunShade="Sun"
original_data$PFT="Tundra"

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record, gsw and CO2r columns were missing (column 2 and 7, respectively).
# I then rename this columns and populate them.

curated_data=original_data[,c("Sample_ID", "Sample_ID","A","Ci","CO2s","Sample_ID","Sample_ID","Patm","Qin","RHs","Tleaf","QC","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "QC" which is not part of the standard and that will be usefull to analyse the
# data quality. 

colnames(curated_data)=c(ESS_column,"QC","Species","Species_type","Biome","SunShade","Dataset","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$gsw = NA
curated_data$Record = NA
curated_data$CO2r = NA

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# gsw is not given in this dataset, I calculate it based on fick's law
curated_data$gsw=1.6*curated_data$A/(curated_data$CO2s-curated_data$Ci)

# Removing the data that the authors flagged as bad data
curated_data=curated_data[is.na(curated_data$QC),]

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
