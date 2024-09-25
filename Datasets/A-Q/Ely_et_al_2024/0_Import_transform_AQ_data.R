## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Rogers et al. 2023



# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)



# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Ely_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Ely_et_al_2024'))

# Import the authors original raw data
colnames_keep=c("hhmmss","SampleID","Species", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","file")
ls_files=dir(file.path(path,'/Datasets/A-Q/Ely_et_al_2024/LI-COR_files'))
ls_files_Licor=ls_files[which(grepl(x=ls_files,pattern=".xlsx")&grepl(x=ls_files,pattern="AQ"))]
original_data=data.frame()
for(file in ls_files_Licor){
  data_file=f.import_licor6800(file = file.path(path,'/Datasets/A-Q/Ely_et_al_2024/LI-COR_files',file),column_display=c('A','gsw','Qin','Ci'))
  if(is.null(data_file$SampleID)){data_file$SampleID=data_file$Barcode}
  data_file$file=file
  data_file=data_file[1:(nrow(data_file)-1),]
  original_data=rbind.data.frame(original_data,data_file[,colnames_keep])
}

# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "obs","A","Ci","CO2_s","CO2_r","gsw","Pa","Qin","RHcham","Tleaf","Species","file")]

# Rename the columns of the curated dataset with the ESS standard
colnames(curated_data)=c(ESS_column,c("Species_Code","file"))
curated_data$Dataset="Ely_et_al_2024"
curated_data$Biome="arctic"
curated_data$Species_type="Wild"
curated_data$SunShade="Sun"
curated_data$PFT="Tundra"
SpeciesInfo=read.csv("SewPenSpecies2023.csv")
SpeciesInfo$Species=paste(SpeciesInfo$Genus,SpeciesInfo$Species)
curated_data=merge(curated_data,SpeciesInfo[,c("USDA.Plants.symbol","Species")],by.x="Species_Code",by.y="USDA.Plants.symbol")

curated_data=curated_data[order(curated_data$SampleID,curated_data$Qin),]

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))
curated_data$Record=as.numeric(curated_data$Record)

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
