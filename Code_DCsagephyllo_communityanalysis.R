#Dry Creek sagebrush phyllosphere project
#Analysis of the microbial communities on sagebrush leaves
#Authors: Jacob Heil, Leonora Bittleston

##Working Directory
setwd("~/Sagebrush_Results/Ch1/Code/DC_Sagebrush")

##Install and load packages
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}
if (!require("effects")) {install.packages("effects"); require("effects")}
if (!require("BiocManager")) {install.packages("BiocManager"); require("BiocManager")}
#09-13-2022, phyloseq not working on R 4.2.1
#if (!require("phyloseq")) {BiocManager::install("phyloseq"); require("phyloseq")}

##Data
#Community Data
data_ASV <- read.csv("~/Sagebrush_Results/Ch1/Data_Community/Data_DCsp_ASV.csv", header = T, row.names = 1, check.names = F)
#Metadata
metadata_plant <- read.csv("~/Sagebrush_Results/Ch1/Data_Meta/Data_DCsp_meta_plant.csv", header = T, row.names = 1, check.names = F)
#############INSERT TAXONOMY AND ANY OTHER DATA HERE
#Clean data
data_ASV <- subset(data_ASV, row.names(data_ASV) %in% row.names(metadata_plant))
#The data sets should be 180 observations

##Naming conventions
#Assign new names to ASVs
data_ASV <- data_ASV[,reorder(colnames(data_ASV), order(colnames(data_ASV)))] #put ASVs in numeric order
ITS_names <- paste("ITS",sprintf('%0.4d', 1:5473), sep = "") #generate new ASV names
namingconv_ITS <- data.frame("Original" = colnames(data_ASV), "New" = ITS_names) 
colnames(data_ASV) <- ITS_names

##Rarefaction
#Rarefy to lowest ASV count for a single sample
set.seed(578) #always set seed to ensure same rarefaction
#find the abundance of the least abundant sample
print(min(rowSums(data_ASV))) #lowest ASV count in ITS data, 4457
#rarefaction curve
rarecurve(data_ASV, step = 10, label = FALSE, main = "ITS rarefaction")
abline(v = 4457, col = "red", lwd = 2)
#Rarefy
data_ASV.r <- data.frame(rrarefy(data_ASV, sample = 4457))
#Remove any ASVs below desired threshold (default set at 0)
data_ASV.r <- data_ASV.r[,colSums(data_ASV.r) > 0]

##Alpha Diversity
#Calculate alpha diversity, add to metadata
metadata_plant$richness.ITS <- rowSums(data_ASV.r > 0) #richness
metadata_plant$abundance.ITS <- rowSums(data_ASV.r) #abundance
metadata_plant$shannon.ITS <- diversity(data_ASV.r) #shannon diversity
metadata_plant$effective.ITS <- round(exp(metadata_plant$shannon.ITS)) #effective number of species

#Add in offset measurements
metadata_plant <- cbind(metadata_plant, "R_offset" = as.numeric(c(rep("NA",6), metadata_plant$Richness[1:150])), "A_offset" = as.numeric(c(rep("NA",6),metadata_plant$Abundance[1:150])), "S_offset" = as.numeric(c(rep("NA",6),metadata_plant$Shannon[1:150])))

#GLM Richness
boxplot(richness.ITS ~ Date, data = metadata_plant)
richdate_glm <- glm(richness.ITS ~ Date, data = metadata_plant)
summary(richdate_glm)

boxplot(richness.ITS ~ Date, data = metadata_plant)
Richness_glm <- glm(richness.ITS ~ R_offset + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C, data = metadata_plant, family = "poisson")
summary(Richness_glm)
plot(allEffects(Richness_glm))

#GLM diversity
boxplot(shannon.ITS ~ Date, data = metadata_plant)
divdate_glm <- glm(shannon.ITS ~ Date, data = metadata_plant)
summary(divdate_glm)

Shannon_glm <- glm(shannon.ITS ~ R_offset + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C, data = metadata_plant)
summary(Shannon_glm)
plot(allEffects(Shannon_glm))

##### Beta Diversity Analyses
#Calculate bray-curtis
dist.ITS <- vegdist(data_ASV.r, method = "bray")

