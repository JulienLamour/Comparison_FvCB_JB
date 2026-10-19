## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data collected in BRF by Davidson et al
## The dataset is available publicly :https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.19137


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R')) 

# Set the working directory to the 'Davidson_et_al_2023' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Davidson_et_al_2023'))

# Import the authors original raw data
original_data=read.csv('BRF_Light_Response_Curves.csv')
TreeID=read.csv('BRF_Leaf_Trait_Data.csv')
SampleInfo=read.csv(('BRF_Tree_Data.csv'))
SampleInfo$Dataset="Davidson_et_al_2023"
SampleInfo$Species=paste(SampleInfo$Genus,SampleInfo$Species)
SampleInfo$Species_type="Wild"
SampleInfo$Biome="Temperate"
SampleInfo$SunShade="Sun"
SampleInfo$PFT="Temperate"
SampleInfo=merge(TreeID,SampleInfo,by="TreeID")
original_data=merge(x=original_data,y=SampleInfo,by="SampleID")

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record column was missing (column 2).
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","Ci","CO2s","CO2s","gsw","Patm","Qin","RHs","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$CO2r = NA
curated_data$Record = NA

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin,decreasing = TRUE),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Converting the conductance into mol m-2 s-1
curated_data$gsw=curated_data$gsw/1000

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
