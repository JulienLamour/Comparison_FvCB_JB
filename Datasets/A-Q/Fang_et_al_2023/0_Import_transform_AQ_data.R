## This codes is used to do the step 0 of the data curation process.
## It imports and transforms the AQ data from https://onlinelibrary.wiley.com/doi/epdf/10.1111/gcb.16488
## The dataset is available publicly at http://doi.org/10.4121/21304917.


# Load the 'here' package to easily reference files and folders path, which is robust and independent on your platform (Linux, Windows..)
library(here)
library(readxl)

# Find the path of the top relative directory
path=here()

# Loading the ESS standard column names that will be used for all datasets
source(file.path(path,'/R models/Tools.R'))

# Set the working directory to the 'Fang_et_al_2023' folder where the data is located
setwd(file.path(path,'/Datasets/A-Q/Fang_et_al_2023'))

# Import the authors original raw data
original_data=read.csv('GasExchange_ChlorophyllFluorescence.csv')
original_data$Dataset="Fang_et_al_2023"
original_data$Species="Triticum aestivum L."
original_data$Species_type="Agricultural"
original_data$Biome="Temperate"
original_data$SunShade="Sun"
original_data$PFT="Crop"
original_data$SampleID=paste(original_data$Experiment,original_data$Genotype,original_data$Treatment,original_data$Replicate)
#I only keep the data at 21%O2
orginal_data=original_data[original_data$'O2....'=="21",]

# Select specific columns from the original data to create the curated data set
# I then rename this columns and populate them.

curated_data=original_data[,c("SampleID", "SampleID","A..mmol.m2.s.","Ci..ppm.","Ca..ppm.","Ca..ppm.","gs..mol.m2.s.","gs..mol.m2.s.","Iinc..mmol.m2.s.","Tleaf....","Tleaf....","Species","Species_type","Biome","SunShade","Dataset","PFT","Treatment","Genotype")]

# Rename the columns of the curated dataset with the ESS standard

colnames(curated_data)=c(ESS_column,"Species","Species_type","Biome","SunShade","Dataset","PFT","Treatment","Genotype")
curated_data$CO2r=NA
curated_data$Record=NA
curated_data$Patm=NA
curated_data$RHs=NA

# I look for the light curves. They start with a light irradiance at 2000
# I also name them with an ID

posi_2000=which(round(curated_data$Qin)==2000)
diff(posi_2000)
AQ_posi=as.vector(sapply(X = posi_2000,FUN = function(x){seq(x,x+9,1)}))
ID_curve=as.vector(sapply(X = posi_2000,FUN = function(x){rep((x-1)/26+1,10)}))
curated_data=curated_data[AQ_posi,]
curated_data$SampleID=paste(curated_data$SampleID,ID_curve)


# Convert the 'SampleID' column to a numerical factor which is simpler to use and refer to for analyzing the data
curated_data$SampleID_num=as.numeric(as.factor(curated_data$SampleID))

# Adding a SampleID_num column
for(leaf in unique(curated_data$SampleID)){
  curated_data[curated_data$SampleID==leaf,"Record"]=1:nrow( curated_data[curated_data$SampleID==leaf,])
}


# Saving the curated dataset
save(curated_data,file='0_curated_data.Rdata')
