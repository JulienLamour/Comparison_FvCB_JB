## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Perez et al. 2024
## The dataset is available publicly :  https://doi.org/10.18167/DVN1/PQAJTV


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Perez_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Perez_et_al_2024'))

# Import the authors original raw data
ls_files=dir(recursive = TRUE) ## This lists all the files that are in the same folder as this R code.

#####################
### Import files  ###
#####################
li6800_colnames=c("id", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","file")

ls_files=dir(recursive = TRUE) ## This lists all the files that are in the same folder as this R code.

## This remove a lot of files that are not LICOR6400 files. You can change the filters as you want.
ls_files_AQ=ls_files[which(grepl(x=ls_files,pattern=".csv")&grepl(x=ls_files,pattern="Licor"))]

### Import all the AQ curves and create a dataframe with all the curves
AQ_data=data.frame()
for(file in ls_files_AQ){
  colname_file=read.csv2(file,skip = 12,nrows = 1)
  file_data=read.csv2(file,skip = 15,header = FALSE,col.names = colname_file,dec = ".")
  if(is.character(file_data$A)){ ## Some files are encoded with "." as decimal and others with ","
    file_data=read.csv2(file,skip = 15,header = FALSE,col.names = colname_file,dec = ",")}
  file_data$file=file
  file_data=file_data[,li6800_colnames]
  AQ_data=rbind.data.frame(AQ_data,file_data)
}


original_data=AQ_data

original_data$Species="Oryza sativa L."
original_data$Species_type="Agricultural"
original_data$Biome="Tropical"
original_data$SunShade="Sun"
original_data$Dataset="Perez_et_al_2024"
original_data$PFT="Crop"

colnames(original_data)=c(ESS_column,"file","Species","Species_type","Biome","SunShade","Dataset","PFT")

curated_data=original_data
# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
