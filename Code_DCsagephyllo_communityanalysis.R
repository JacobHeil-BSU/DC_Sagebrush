#Dry Creek sagebrush phyllosphere project
#Analysis of the microbial communities on sagebrush leaves
#Authors: Jacob Heil, Leonora Bittleston

setwd("~/Sagebrush_Results/Ch1/Code/DC_Sagebrush")

#Install and load packages
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}

#Load in Data
#Community Data
data_ASV <- read.csv("~/Sagebrush_Results/Ch1/Data_Community/Data_DCsp_ASV.csv", header = T, row.names = 1, check.names = F)
#Metadata
metadata_plant <- read.csv("~/Sagebrush_Results/Ch1/Data_Meta/Data_DCsp_meta_plant.csv", header = T, row.names = 1, check.names = F)
#############INSERT TAXONOMY AND ANY OTHER DATA HERE
#Clean data
data_ASV <- subset(data_ASV, row.names(data_ASV) %in% row.names(metadata_plant))
#The data sets should be 180 observations

#Assign new names to ASVs
data_ASV <- data_ASV[,reorder(colnames(data_ASV), order(colnames(data_ASV)))] #put ASVs in numeric order
ITS_names <- paste("ITS",sprintf('%0.4d', 1:5473), sep = "") #generate new ASV names
namingconv_ITS <- data.frame("Original" = colnames(data_ASV), "New" = ITS_names) 
colnames(data_ASV) <- ITS_names
