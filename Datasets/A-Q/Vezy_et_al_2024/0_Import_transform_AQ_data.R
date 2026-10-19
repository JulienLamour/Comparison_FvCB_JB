## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from Vezy et al. 2024
## The dataset is available publicly : https://doi.org/10.5281/zenodo.12704284


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Vezy_et_al_2024' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Vezy_et_al_2024'))

# Import the authors original raw data
ls_files=dir(recursive = TRUE) ## This lists all the files that are in the same folder as this R code.

## TRemoving one file with issues
ls_files_Walz=ls_files[which(grepl(x=ls_files,pattern=".csv")&!grepl(x=ls_files,pattern="P5F70408.csv"))]

col_walz=c("Date","Time","Code","Object","Area","Status","Comment","CO2abs","dCO2ZP","dCO2MP","H2Oabs","dH2OZP","dH2OMP","Flow","Pamb","Aux1","Aux2","Tcuv","Tleaf","Tamb","Tmin","PARtop","PARbot","PARamb","Imp","rh","E","VPD","GH2O","A","ci","ca","wa","Fo","Fm","Fv/Fm","F","Fm'","Fo'","Fo'calc","Yield","ETR","qP","qL","qN","NPQ","Y(NPQ)","ETR-Fac")
col_walz=make.names(col_walz)

## Import and bind the data from all the files
data_Walz=data.frame()
for(file in ls_files_Walz){
  print(file)
  data_file=read.csv(file,skip=2,col.names = col_walz,sep=";")
  data_file$file=file
  status="CO2 Curve"
  for(line in 1:nrow(data_file)){
    if((data_file[line,"Comment"])!=""){status=data_file[line,"Comment"]}
    data_file[line,"Comment"]=status
    data_file$record=line
  }
  data_Walz=rbind.data.frame(data_Walz,data_file)
}
data_Walz$GH2O=data_Walz$GH2O/1000

## Here, we only keep the light curves
data_Walz=data_Walz[data_Walz$Comment=="ligth Curve",]
## and we only keep the last point of each light level, when fluorescence was measured
data_Walz=data_Walz[!is.na(data_Walz$Yield),]
## Giving an observation name for each file
for (file in unique(data_Walz$file)){
  data_Walz[data_Walz$file==file,"record"]=1:nrow(data_Walz[data_Walz$file==file,])
}


original_data=data_Walz
original_data$Dataset="Vezy_et_al_2024"
original_data$Species="Elaeis guineensis"
original_data$Species_type="Agricultural"
original_data$Biome="Tropical"
original_data$SunShade="Sun"
original_data$PFT="Crop"

# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("file", "record","A","ci","ca","ca","GH2O","Pamb","PARtop","rh","Tleaf","Species","Species_type","Biome","SunShade","Dataset","PFT")]

# Rename the columns of the curated dataset with the ESS standard
# I also kept the column "Remove" which is not part of the standard and that will be usefull to analyse the
# data quality. It corresponds to the author flagging of bad data. THis column doesnt have to be included
# in other datasets.
colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT")
curated_data$CO2r=NA

# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
