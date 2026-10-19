## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data collected by Burnett et al. 2019 


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Burnett_et_al_2019' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Burnett_et_al_2019'))

# Import the authors original raw data
original_data=read.csv('GE_sample_AQ.csv')  # https://doi.org/10.1111/pce.13574
original_data$SunShade="Sun"
original_data$PFT="Crop"
original_data2=read_xlsx(sheet=1,"Zucchini_AQ.xlsx") 
original_data2$Species_type="Agricultural"
original_data2$Biome="Greenhouse"
original_data2$SunShade="Sun"
original_data2$PFT="Crop"

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record column was missing (column 2).
# I then rename this columns and populate them.

curated_data=original_data[,c("Species", "obs","photo","ci","co2s","co2r","cond","press","pari","rhs","tleaf","Species","Species_type","Biome","SunShade","PFT")]
curated_data2=original_data2[,c("Leaf", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","Species","Species_type","Biome","SunShade","PFT")]

# Rename the columns of the curated dataset with the ESS standard
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","PFT")
colnames(curated_data2)=c(ESS_column,"Species","Species_type","Biome","SunShade","PFT")
curated_data=rbind.data.frame(curated_data,curated_data2)
curated_data$Dataset="Burnett_et_al_2019"

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Removing the last duplicated point of the curves that was measured at high irradiance just after measurement in the dark
QC_table=cbind.data.frame(SampleID_num=c(1,2,3,4),
                          Record=c(19,21,23,3)) 
curated_data$QC="ok"
curated_data[paste(curated_data$SampleID_num,curated_data$Record)%in%paste(QC_table$SampleID_num,QC_table$Record),'QC']='bad'
curated_data=curated_data[curated_data$QC=='ok',]

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
