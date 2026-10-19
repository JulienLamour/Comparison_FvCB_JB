## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Fang_et_al_2026.
## The dataset was sent directly by Liang Fang


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Fang_et_al_2026' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Fang_et_al_2026'))

# Import the authors original raw data
original_data=read.csv('AQ curves of 32 tree species.csv')

original_data$Dataset="Fang_et_al_2026"
original_data$Species=original_data$Scientific.name
original_data$Species_type="Wild"
original_data$Biome="Subtropical"
original_data$SunShade="Sun"
original_data$PFT="Subtropical"
original_data$SampleID = paste(original_data$Dataset,original_data$Species, original_data$Replicate)

# Select specific columns from the original data to create the curated data set
# Here I selected several times the column SampleID. This is because the Record and gsw columns were missing (column 2 and 7, respectively).
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A","Ci","Ca","SampleID","gsw","SampleID","Qin","SampleID","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")

# The columns gsw and Record were still filled with SampleID info. I replace by NA values
curated_data$CO2r = NA
curated_data$Patm = NA
curated_data$RHs = NA
curated_data$Record = NA

hist(curated_data$Ci)
curated_data = curated_data[curated_data$Ci>0,]
for(SampleID_num in unique(curated_data$SampleID_num)){
  AQleaf = curated_data[curated_data$SampleID_num == SampleID_num,]
  plot(x = AQleaf$Qin, y = AQleaf$A, main = SampleID_num)
}

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

curated_data = curated_data[-which(curated_data$SampleID_num%in%c(5,11,23,44,45,58,67)),]

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]
# Populate the 'Record' column with a unique number for each observation of each leaf 'SampleID'
for(SampleID in curated_data$SampleID){
  curated_data[curated_data$SampleID==SampleID,'Record']=1:nrow(curated_data[curated_data$SampleID==SampleID,])
}



# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')