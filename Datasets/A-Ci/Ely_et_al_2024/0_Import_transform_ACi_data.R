## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the ACi data from Ely et al 2024.
## The dataset is available publicly :Ely K ; Yang D ; Anderson J ; Serbin S ; Rogers A (2024): Plant physiology, shrub size, thaw depth and soil water content, Seward Peninsula, Alaska, 2023. Next-Generation Ecosystem Experiments (NGEE) Arctic, ESS-DIVE repository. Dataset. doi:10.15485/2341585


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names (variable ESS_column) from the "Tools.R" file
source(file.path(path,'/R models/Tools.R')) 

# Set the working directory to the 'Ely_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Ci/Ely_et_al_2024'))

# Import the authors original raw data
colnames_keep=c("hhmmss","SampleID","Species", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","file")
ls_files=dir(file.path(path,'/Datasets/A-Ci/Ely_et_al_2024/LI-COR_files'))
ls_files_Licor=ls_files[which(grepl(x=ls_files,pattern=".xlsx")&grepl(x=ls_files,pattern="ACi"))]
original_data=data.frame()
for(file in ls_files_Licor){
  data_file=f.import_licor6800(file = file.path(path,'/Datasets/A-Ci/Ely_et_al_2024/LI-COR_files',file),column_display=c('A','gsw','Qin','Ci'))
  if(is.null(data_file$SampleID)){data_file$SampleID=data_file$Barcode}
  data_file$file=file
  original_data=rbind.data.frame(original_data,data_file[,colnames_keep])
}


# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","Species","file")]

# Rename the columns of the curated dataset with the ESS standard
colnames(curated_data)=c(ESS_column,c("Species","file"))


curated_data=curated_data[order(curated_data$SampleID,curated_data$Record),]

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))
curated_data$Record=as.numeric(curated_data$Record)
curated_data$Dataset="Ely_et_al_2024"

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
