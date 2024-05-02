#Leaf age and air temperature drive change in the sagebrush phyllosphere over the course of a full year

#Set up environment ####
##Install and load packages
pkgs <- c("vegan", "ggplot2", "effects", "AER", "MASS", "BiocManager", "brms", "phyloseq", "ggcorrplot", "factoextra", "ANCOMBC", "this.path", "dbstats", "ggeffects")

# Load and install required packages
for (i in pkgs) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

### Working Directory 
setwd(this.path::here())

#Load and clean data ####
#Community Data from sequencing
data_ASV <- read.csv("../../Data_Community/Data_DCsp_ASV.csv", header = TRUE, row.names = 1, check.names = FALSE)
#Metadata
metadata_plant <- read.csv("../../Data_Meta/Data_DCsp_meta_plant.csv", header = T, row.names = 1, check.names = F)
#Taxonomy
data_tax <- read.csv("newNCBIfulltax.csv", row.names = 1, header = TRUE)
data_tax <- data_tax[reorder(row.names(data_tax), order(row.names(data_tax))),] 
#Remove outliers
#metadata_plant <- metadata_plant[! row.names(metadata_plant) %in% c("DC31PA", "DC16EC", "DC15ED", "DC18PB"),]
#Clean data
data_ASV <- subset(data_ASV, row.names(data_ASV) %in% row.names(metadata_plant))
data_ASV <- data_ASV[,colSums(data_ASV) > 0]
data_tax <- subset(data_tax, row.names(data_tax) %in% colnames(data_ASV))
#reorder date factor to reflect chronology of sampling dates
metadata_plant$Date <- factor(metadata_plant$Date, levels= unique(metadata_plant$Date))
#reorder season factor to reflect chronology of sampling dates
metadata_plant$Season <- factor(metadata_plant$Season, levels = unique(metadata_plant$Season))
#reorder type factor to reflect chronology of sampling dates
metadata_plant$Type <- factor(metadata_plant$Type, levels = c("Persistent", "Ephemeral"))
#Assign new names to ASVs
data_ASV <- data_ASV[,reorder(colnames(data_ASV), order(colnames(data_ASV)))] 
#put ASVs in numeric order
ITS_names <- paste("ITS",sprintf('%0.4d', 1:4998), sep = "") #generate new ASV names
namingconv_ITS <- data.frame("Original" = colnames(data_ASV), "New" = ITS_names) 
colnames(data_ASV) <- ITS_names
row.names(data_tax) == namingconv_ITS$Original
row.names(data_tax) <- namingconv_ITS$New

#Absolute Quantification ####
#Function to compute DNA copy number in a sample
DNAcopynumber <- function(concentration, bp){
  (concentration * 6.022e+23) / (bp * 1e+9 * 650)
}

#Calculate the copy number for each sample and load into meta
metadata_plant <- cbind(metadata_plant, "Copy_Number" = sapply(metadata_plant$`Concentration_mean(ng/ul)`, DNAcopynumber, bp = 290))

#Calculate copy number for each ASV in each sample
data_ASV_abs <- c()
for (a in 1:nrow(data_ASV)){
  abs <- (data_ASV[a,] / sum(data_ASV[a,]) * metadata_plant$Copy_Number[a])
  data_ASV_abs <- rbind(data_ASV_abs, abs)
  remove(abs, a)
}

#Remove the first sample and Remove any ASVs below desired threshold (set at 0) and
data_ASV_abs <- data_ASV_abs[2:180,]
data_ASV_abs <- round(data_ASV_abs)
data_ASV_abs <- data_ASV_abs[,colSums(data_ASV_abs) > 0]
#match meta to data
metadata_plant <- metadata_plant[2:180,]
data_tax <- subset(data_tax, row.names(data_tax) %in% colnames(data_ASV_abs))
#adjust [Chalara] clidemiae
data_tax[data_tax$species %in% "[Chalara] clidemiae",]$genus <- "Chalara"

#Compute Alpha Diversity ####
#Calculate alpha diversity, add to metadata
metadata_plant$richness.ITS <- rowSums(data_ASV_abs > 0) #richness
metadata_plant$abundance.ITS <- rowSums(data_ASV_abs) #abundance
metadata_plant$shannon.ITS <- diversity(data_ASV_abs) #shannon diversity
metadata_plant$effective.ITS <- round(exp(metadata_plant$shannon.ITS)) #effective number of species
#Create Phyloseq object (for ANCOMBC)####
#add asv names to end of species for ancombc at asv level
data_tax2 <- data_tax
data_tax2$species <- paste(data_tax2$species, row.names(data_tax2))

#Combine all tables into a phyloseq object
data_tax_physeq <- tax_table(data_tax2)
row.names(data_tax_physeq) <- row.names(data_tax)
colnames(data_tax_physeq) <- colnames(data_tax)
data_phyloseq_all <- phyloseq(otu_table = otu_table(t(data_ASV), taxa_are_rows = TRUE), tax_table = data_tax_physeq, sample_data(metadata_plant))

#Basic Descriptive Stats####
#distribution of concentration and copy number
hist(metadata_plant$`Concentration_mean(ng/ul)`, main = "Concentration by sample", xlab = "Concentration")
summary(metadata_plant$`Concentration_mean(ng/ul)`)
hist(metadata_plant$Copy_Number, main = "Copy number by sample", xlab = "Copy number")
summary(metadata_plant$Copy_Number)
hist(rowSums(data_ASV_abs), xlab = "Copy Number", main = "Absolute Abundance, copy # per sample")
summary(rowSums(data_ASV_abs))
hist(colSums(data_ASV_abs), xlab = "Copy Number", main = "Absolute Abundance, copy # per ASV in whole data set")
summary(colSums(data_ASV_abs))

#Basic ASV identity stats including taxa barplots ####
data_ASV_PA_sample <- sort((colSums(data_ASV_abs > 0) / 179), decreasing = TRUE)
hist(data_ASV_PA_sample, breaks = 6, labels = TRUE, ylim = c(0,5100), xlab = "presence frequency", main = "Presence/Absence of all ASVs by sample")
text(x = .3, y = 4000, labels = "ASVs w/ 90+% presence:", pos = 4)
text(x = .3, y = 3700, labels = "ITS0185, Aureobasidium leucospermi, 100% presence", pos = 4)
text(x = .3, y = 3400, labels = "ITS0041, [Chalara] clidemiae, 99.4% presence", pos = 4)
text(x = .3, y = 3100, labels = "ITS2632, Endoconidioma populi, 99.4% presence", pos = 4)
text(x = .3, y = 2800, labels = "ITS0204, Taphrina sp., 97.2% presence", pos = 4)
text(x = .3, y = 2500, labels = "ITS1136, Kabatina thujae, 90.5% presence", pos = 4)
text(x = .3, y = 2200, labels = "ITS4184, Phaeococcomyces eucalypti, 90.5% presence", pos = 4)
polygon(x = c(.3,.3,.91,.91), y = c(4200, 2000, 2000, 4200))

#Presence/absence in all dates
data_ASV_PA_date <- c()
for (i in unique(metadata_plant$Date)){
  a <- subset(data_ASV_abs, metadata_plant$Date %in% i)
  a <- colSums(a) > 0
  data_ASV_PA_date <- rbind(data_ASV_PA_date, a)
  remove(a, i)
}
row.names(data_ASV_PA_date) <- unique(metadata_plant$Date)
data_ASV_PA_date <- sort((colSums(data_ASV_PA_date > 0) / 26), decreasing = TRUE)
hist(data_ASV_PA_date, breaks = 6, labels = TRUE, ylim = c(0,5000), xlab = "presence frequency", main = "Presence/Absence of all ASVs by date")
text(x = .3, y = 5000, cex = 1, labels = "ASVs w/ 90+% presence:", pos = 4)
text(x = .3, y = 4800, cex = 1, labels = "ITS1219, Alternaria alstroemeriae, 100% presence", pos = 4)
text(x = .3, y = 4600, cex = 1, labels = "ITS0185, Aureobasidium leucospermi, 100% presence", pos = 4)
text(x = .3, y = 4400, cex = 1, labels = "ITS0041, [Chalara] clidemiae, 100% presence", pos = 4)
text(x = .3, y = 4200, cex = 1, labels = "ITS0089, Comoclathris lini, 100% presence", pos = 4)
text(x = .3, y = 4000, cex = 1, labels = "ITS4019, Didymocyrtis banksiae, 100% presence", pos = 4)
text(x = .3, y = 3800, cex = 1, labels = "ITS2632, Endoconidioma populi, 100% presence", pos = 4)
text(x = .3, y = 3600, cex = 1, labels = "ITS1136, Kabatina thujae, 100% presence", pos = 4)
text(x = .3, y = 3400, cex = 1, labels = "ITS0997, Kabatina thujae, 100% presence", pos = 4)
text(x = .3, y = 3200, cex = 1, labels = "ITS0567, Parafenestella mackenziei, 100% presence", pos = 4)
text(x = .3, y = 3000, cex = 1, labels = "ITS4184, Phaeococcomyces eucalypti, 100% presence", pos = 4)
text(x = .3, y = 2800, cex = 1, labels = "ITS4385, Phaeococcomyces eucalypti, 100% presence", pos = 4)
text(x = .3, y = 2600, cex = 1, labels = "ITS0204, Taphrina sp., 100% presence", pos = 4)
text(x = .3, y = 2400, cex = 1, labels = "ITS1525, Taphrina antarctica, 100% presence", pos = 4)
text(x = .3, y = 2200, cex = 1, labels = "ITS2227, Aureobasidium leucospermi, 96% presence", pos = 4)
text(x = .3, y = 2000, cex = 1, labels = "ITS1176, Chaetothyrium agathis, 96% presence", pos = 4)
text(x = .3, y = 1800, cex = 1, labels = "ITS1815, Comoclathris lini, 96% presence", pos = 4)
text(x = .3, y = 1600, cex = 1, labels = "ITS3758, Neocelosporium eucalypti, 96% presence", pos = 4)
text(x = .3, y = 1400, cex = 1, labels = "ITS1223, Phaeococcomyces eucalypti, 96% presence", pos = 4)
text(x = .3, y = 1200, cex = 1, labels = "ITS2811, Staurosphaeria rhamnicola, 96% presence", pos = 4)
text(x = .3, y = 1000, cex = 1, labels = "ITS2872, Chaetothyrium agathis, 96% presence", pos = 4)
polygon(x = c(.3,.3,.95,.95), y = c(5100, 900, 900, 5100))

#Presence absence of genera by date
data_ASV_PA_genus <- data.frame("tax" = data_tax$genus, t(data_ASV_abs))
data_ASV_PA_genus <- data.frame(t(aggregate(. ~ data_ASV_PA_genus[,1], data_ASV_PA_genus[,2:ncol(data_ASV_PA_genus)], sum)))
colnames(data_ASV_PA_genus) <- data_ASV_PA_genus[1,]
data_ASV_PA_genus <- data_ASV_PA_genus[2:180,]
data_ASV_PA_genus[,1:364] <- lapply(data_ASV_PA_genus[,1:364], as.numeric)
data_ASV_PA_genus_date <- c()
for (i in unique(metadata_plant$Date)){
  a <- subset(data_ASV_PA_genus, metadata_plant$Date %in% i)
  a[,1:364] <- lapply(a[,1:364], as.numeric)
  a <- colSums(a) > 0
  data_ASV_PA_genus_date <- rbind(data_ASV_PA_genus_date, a)
  remove(a, i)
}
row.names(data_ASV_PA_genus_date) <- unique(metadata_plant$Date)
data_ASV_PA_genus_date <- sort((colSums(data_ASV_PA_genus_date > 0) / 26), decreasing = TRUE)
hist(data_ASV_PA_genus_date, breaks = 6, labels = TRUE, ylim = c(0,300), xlab = "presence frequency", main = "Presence/Absence of all genera by date")

#Generate taxa barplots
for (a in 3:ncol(data_tax)){
  data <- data.frame("tax" = data_tax[,a], t(data_ASV_abs))
  data <- aggregate(. ~ data[,1], data[,2:ncol(data)], sum) #aggregate all ASVs by taxonomic level
  row.names(data) <- data[,1]; data <- data[,2:ncol(data)]
  data <- data[order(rowSums(data[,2:ncol(data)]),decreasing = T),]
  data_prop <- c()
  for (x in 1:ncol(data)){
    data_prop <- cbind(data_prop, data[,x] / colSums(data)[x])}
  row.names(data_prop) <- row.names(data)
  colnames(data_prop) <- colnames(data)
  data_prop_top <- data_prop[1:10,] # isolate top 10 most abundant
  other <- colSums(data_prop[11:nrow(data_prop),], na.rm = TRUE) #bin all others in one group
  data_prop_final <- rbind(data_prop_top, "other" = other)
  data_prop_final <- as.matrix(data_prop_final)
  colnames(data_prop_final) <- colnames(data_prop)
  customcol <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
                 "cyan","darkred","purple","mediumblue","palegoldenrod",
                 "lightgoldenrod","indianred","yellow","purple4","darkgreen",
                 "lightsalmon","yellow3","purple2","lightblue","firebrick",
                 "navy","red4","red","darkmagenta","mediumvioletred",
                 "violetred2","skyblue","dodgerblue4")
  barplot(data_prop_final, col=customcol, main = paste("Taxa barplot at", colnames(data_tax)[a], "level", sep = " "),
          legend.text = T, axes = F, cex.names = .8, las = 2, border=NA, space=-0.1,
          args.legend = list(legend = rev(rownames(data_prop_final)), x = "topleft", bty = "n", text.col = "white"))
  remove(data, a, x, customcol)
}

#Generate genus barplot at date
  genus_bar <- data.frame("tax" = data_tax$genus, t(data_ASV_abs))
  genus_bar <- aggregate(. ~ genus_bar[,1], genus_bar[,2:ncol(genus_bar)], sum) #aggregate all ASVs by taxonomic level
  row.names(genus_bar) <- genus_bar[,1]; genus_bar <- genus_bar[,2:ncol(genus_bar)]
  genus_bar <- genus_bar[order(rowSums(genus_bar[,2:ncol(genus_bar)]),decreasing = T),]
  genus_bar_prop <- c()
  for (x in unique(metadata_plant$Date)){
  genus_bar_x <- subset(t(genus_bar), row.names(t(genus_bar)) %in% row.names(metadata_plant[metadata_plant$Date == x,]))
  genus_bar_prop <- rbind(genus_bar_prop, colSums(genus_bar_x) / sum(genus_bar_x))
  }
  row.names(genus_bar_prop) <- unique(metadata_plant$Date)
  taxa_p <- c("Alternaria", "Aureobasidium", "Chaetothyrium", "Chalara", "Cladosporium", "Comoclathris", "Didymocyrtis", "Endoconidioma", "Filobasidium", "Kabatina", "Meristemomyces", "Neocelosporium", "Parafenestella", "Penidiella", "Phaeococcomyces", "Preussia", "Robertozyma", "Staurosphaeria", "Taphrina", "Thyrostroma")
  genus_bar_prop_top <- genus_bar_prop[,colnames(genus_bar_prop) %in% taxa_p] # isolate top 10 most abundant
  other <- rowSums(genus_bar_prop[,-which(colnames(genus_bar_prop) %in% taxa_p)], na.rm = TRUE) #bin all others in one group
  genus_bar_prop_final <- cbind(genus_bar_prop_top, "other" = data.frame(other))
  genus_bar_prop_final <- t(genus_bar_prop_final)
  #colnames(data_prop_final) <- colnames(data_prop)
  customcol <- customcol <- c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
                              "cyan","darkred","purple","mediumblue","palegoldenrod",
                              "indianred","yellow","purple4","darkgreen","lightgoldenrod",
                              "lightsalmon","yellow3","purple2","lightblue","firebrick",
                              "navy","red4","red","darkmagenta","mediumvioletred",
                              "violetred2","skyblue","dodgerblue4")
levels_g <- row.names(genus_bar_prop_final)
barplot(genus_bar_prop_final, col=customcol, main = "Genus barplot by date", legend.text = T, axes = F, cex.names = .8, las = 2, border=NA, space=-0.1, args.legend = list(legend = rev(levels_g), x = "topleft", bty = "n", text.col = "black"))

#look at genera not present at all dates
metadata_plant <- cbind(metadata_plant, "others" = data_prop_final["other",])
boxplot(others ~ Date, metadata_plant)

#PCA ####
data_prop1 <- as.data.frame(data_prop[names(sort(rowSums(data_prop), decreasing = TRUE)),])
##PCA of explanatory variable
#subset variables of interest
metadata_PCA <- metadata_plant[,c(8,17,18,32,33,38)]
metadata_PCA <- na.omit(metadata_PCA)
metadata_PCA <- scale(metadata_PCA)
corr_PCA <- cor(metadata_PCA)
ggcorrplot(corr_PCA)
actual_PCA <- princomp(corr_PCA)
summary(actual_PCA)
fviz_eig(actual_PCA, addlabels = TRUE)
fviz_pca_var(actual_PCA, col.var = "black")
fviz_cos2(actual_PCA, choice = "var", axes = 1:2)
#how much each variable is represented in a given component
#high cos2 green
#mid cos3 orange
#low cos3 black
#~8% variation exlained by y
fviz_pca_var(actual_PCA, col.var = "y",
             gradient.cols = c("black", "yellow", "blue"),
             repel = TRUE)

#~78% variation explained by x
fviz_pca_var(actual_PCA, col.var = "cos2",
             gradient.cols = c("black", "yellow", "blue"),
             repel = TRUE)

#Compare culture and sequence data and response variables to each other ####
cor.test(metadata_plant$culture_richness, metadata_plant$richness.ITS, method = "pearson")
plot(culture_richness ~ richness.ITS, data = metadata_plant, main = "Correlation between culture richness and sequencing richness", ylab = "culture richness", xlab = "sequencing richness")
abline(lm(culture_richness ~ richness.ITS, data = metadata_plant), col = "red")

cor.test(metadata_plant$culture_abundance, metadata_plant$abundance.ITS, method = "pearson")
plot(culture_abundance ~ abundance.ITS, data = metadata_plant, main = "Correlation between culture abundance and sequencing abundance", ylab = "culture abundance", xlab = "sequencing abundance")
abline(lm(culture_abundance ~ abundance.ITS, data = metadata_plant), col = "red")

#compare richness and abundance
cor.test(metadata_plant$abundance.ITS, metadata_plant$richness.ITS, method = "pearson")
plot(abundance.ITS ~ richness.ITS, data = metadata_plant, main = "Correlation between Abundance and Richness (sequencing data)", ylab = "Abundance", xlab = "Richness")
abline(lm(abundance.ITS ~ richness.ITS, data = metadata_plant), col = "red")

#compare effective and abundance
cor.test(metadata_plant$abundance.ITS, metadata_plant$effective.ITS, method = "pearson")
plot(abundance.ITS ~ effective.ITS, data = metadata_plant, main = "Correlation between Abundance and Richness (sequencing data)", ylab = "Abundance", xlab = "Effective Number of Speciess")
abline(lm(abundance.ITS ~ effective.ITS, data = metadata_plant), col = "red")

#compare richness and effective
cor.test(metadata_plant$effective.ITS, metadata_plant$richness.ITS, method = "pearson")
plot(effective.ITS ~ richness.ITS, data = metadata_plant, main = "Correlation between Abundance and Richness (sequencing data)", ylab = "Abundance", xlab = "Richness")
abline(lm(effective.ITS ~ richness.ITS, data = metadata_plant), col = "red")

#compare abundance and concentration
cor.test(metadata_plant$abundance.ITS, metadata_plant$`Concentration_mean(ng/ul)`, method = "pearson")
plot(abundance.ITS ~ `Concentration_mean(ng/ul)`, data = metadata_plant, main = "Correlation between Abundance and Concentration (sequencing data)", ylab = "Abundance", xlab = "Concentration")
abline(lm(abundance.ITS ~ `Concentration_mean(ng/ul)`, data = metadata_plant), col = "red")


#response variables on ephemeral vs. persistent leaves ####
#richness
boxplot(richness.ITS ~ Type, data = metadata_plant, ylab = "Richness")
stripchart(richness.ITS ~ Type, data = metadata_plant, method = "jitter", pch = 19, col = 4, vertical = TRUE, add = TRUE)
t.test(richness.ITS ~ Type, data = metadata_plant)
#abundance
boxplot(log(abundance.ITS + 1) ~ Type, data = metadata_plant, ylab = "Abundance")
stripchart(log(abundance.ITS + 1) ~ Type, data = metadata_plant, method = "jitter", pch = 19, col = 4, vertical = TRUE, add = TRUE)
t.test(abundance.ITS ~ Type, data = metadata_plant)
#Effective number of species
boxplot(effective.ITS ~ Type, data = metadata_plant, ylab = "Effective Number of Species")
stripchart(effective.ITS ~ Type, data = metadata_plant, method = "jitter", pch = 19, col = 4, vertical = TRUE, add = TRUE)
t.test(effective.ITS ~ Type, data = metadata_plant)

#Response variable Change over time ####
#Concentration over time 
ggplot(data = metadata_plant, mapping = aes(x = Date, y = log(`Concentration_mean(ng/ul)`), fill = Type)) + 
  geom_boxplot() + scale_fill_manual(values = c("grey", "white")) + geom_jitter(pch = 1) +
  labs(y = "Concentration", x = "", title = "Concentration over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#richness vs time as continuous variable 
cor.test(metadata_plant$richness.ITS, metadata_plant$Date_cont)

#richness over time boxplot
ggplot(data = metadata_plant, mapping = aes(x = Date, y = richness.ITS, fill = Type)) + 
  geom_boxplot() + scale_fill_manual(values = c("#B4EEB4","#698B69")) +
  labs(y = "ASV Richness", x = "", title = "ASV Richness over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#abundance vs time as continuous variable 
#raw abund
cor.test(metadata_plant$abundance.ITS, metadata_plant$Date_cont)
cor.test(metadata_plant[2:nrow(metadata_plant),]$abundance.ITS, metadata_plant[2:nrow(metadata_plant),]$Date_cont)
#log abund
cor.test(log(metadata_plant$abundance.ITS), metadata_plant$Date_cont)

#abundance over time
ggplot(data = metadata_plant, mapping = aes(x = Date, y = log(abundance.ITS), fill = Type)) +
  geom_boxplot(data = NULL) + scale_fill_manual(values = c("#B4EEB4","#698B69")) +
  labs(y = "log Total Abundance", x = "", title = "Total Abundance over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#richness vs time as continuous variable 
cor.test(metadata_plant$effective.ITS, metadata_plant$Date_cont)

ggplot(data = metadata_plant, mapping = aes(x = Date, y = effective.ITS, fill = Type)) + 
  geom_boxplot() + scale_fill_manual(values = c("orange","skyblue")) +
  labs(y = "log ASV Abundance", x = "", title = "ASV Abundance over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#GLMs with time as explanatory variable ####
#no separation by leaf type
richdate_glm <- glm(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant, family = "poisson")
dispersiontest(richdate_glm) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
richdate_glm <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(richdate_glm)
View(anova(richdate_glm))
plot(richdate_glm)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant$Date <- factor(metadata_plant$Date, levels = c("3/12/2021", "3/26/2021", "4/9/2021", "4/23/2021", "5/7/2021", "5/22/2021", "6/3/2021", "6/18/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "8/27/2021", "9/10/2021", "9/24/2021", "10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
richdate_glm <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(richdate_glm)
richdate_summ <- data.frame(summary(richdate_glm_ephemeral)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(richdate_summ, richdate_summ$Pr...z.. > 5e-2))

#Richness over time ephemeral leaves only
metadata_plant_ephemeral <- subset(metadata_plant, metadata_plant$Type == "Ephemeral")
richdate_glm_ephemeral <- glm(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral, family = "poisson")
dispersiontest(richdate_glm_ephemeral) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
richdate_glm_ephemeral <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(richdate_glm_ephemeral)
plot(richdate_glm_ephemeral)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_ephemeral$Date <- factor(metadata_plant_ephemeral$Date, levels = c("7/2/2021", "4/23/2021", "5/22/2021", "5/7/2021", "6/3/2021", "6/18/2021", "7/16/2021", "7/30/2021"))
richdate_glm_ephemeral <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(richdate_glm_ephemeral)
richdate_summ <- data.frame(summary(richdate_glm_ephemeral)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(richdate_summ, richdate_summ$Pr...z.. > .05))

#Richness over time persistent leaves only
metadata_plant_persistent <- subset(metadata_plant, metadata_plant$Type == "Persistent")
richdate_glm_persistent <- glm(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent, family = "poisson")
dispersiontest(richdate_glm_persistent) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
richdate_glm_persistent <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(richdate_glm_persistent)
plot(richdate_glm_persistent)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_persistent$Date <- factor(metadata_plant_persistent$Date, levels = c("9/24/2021", "3/12/2021", "3/26/2021", "4/23/2021", "4/9/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "8/27/2021", "9/10/2021", "10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
richdate_glm_persistent <- glm.nb(richness.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(richdate_glm_persistent)
richdate_summ <- data.frame(summary(richdate_glm_persistent)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(richdate_summ, richdate_summ$Pr...z.. > .05))

#abundance over time boxplot
ggplot(data = metadata_plant, mapping = aes(x = Date, y = log(abundance.ITS), fill = Type)) + 
  geom_boxplot() + scale_fill_manual(values = c("orange","skyblue")) +
  labs(y = "log ASV abundance", x = "", title = "ASV abundance over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#GLM abundance
#no separation by leaf type
abundate_glm <- glm(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant, family = "poisson")
dispersiontest(abundate_glm) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
abundate_glm <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(abundate_glm)
View(anova(abundate_glm))
plot(abundate_glm)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant$Date <- factor(metadata_plant$Date, levels = c("3/12/2021", "3/26/2021", "4/9/2021", "4/23/2021", "5/7/2021", "5/22/2021", "6/3/2021", "6/18/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "8/27/2021", "9/10/2021", "9/24/2021", "10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
abundate_glm <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(abundate_glm)
abundate_summ <- data.frame(summary(abundate_glm_ephemeral)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(abundate_summ, abundate_summ$Pr...z.. > 5e-2))

#abundance over time ephemeral leaves only
metadata_plant_ephemeral <- subset(metadata_plant, metadata_plant$Type == "Ephemeral")
abundate_glm_ephemeral <- glm(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral, family = "poisson")
dispersiontest(abundate_glm_ephemeral) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
abundate_glm_ephemeral <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(abundate_glm_ephemeral)
plot(abundate_glm_ephemeral)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_ephemeral$Date <- factor(metadata_plant_ephemeral$Date, levels = c("4/23/2021", "5/22/2021", "5/7/2021", "6/3/2021", "6/18/2021", "7/2/2021", "7/16/2021", "7/30/2021"))
abundate_glm_ephemeral <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(abundate_glm_ephemeral)
abundate_summ <- data.frame(summary(abundate_glm_ephemeral)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(abundate_summ, abundate_summ$Pr...z.. > .05))

#abundance over time persistent leaves only
metadata_plant_persistent <- subset(metadata_plant, metadata_plant$Type == "Persistent")
abundate_glm_persistent <- glm(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent, family = "poisson")
dispersiontest(abundate_glm_persistent) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
abundate_glm_persistent <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(abundate_glm_persistent)
plot(abundate_glm_persistent)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_persistent$Date <- factor(metadata_plant_persistent$Date, levels = c("8/27/2021", "3/12/2021", "3/26/2021", "4/9/2021", "4/23/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "9/10/2021", "9/24/2021","10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
abundate_glm_persistent <- glm.nb(abundance.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(abundate_glm_persistent)
abundate_summ <- data.frame(summary(abundate_glm_persistent)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(abundate_summ, abundate_summ$Pr...z.. > .05))

#Effective number of species over time boxplot
ggplot(data = metadata_plant, mapping = aes(x = Date, y = effective.ITS, fill = Type)) + 
  geom_boxplot() + scale_fill_manual(values = c("#B4EEB4","#698B69")) +
  labs(y = "Effective number of species", x = "", title = "Effective number of species over time") +
  theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

#Air temperature over time
ggplot(data = metadata_plant, mapping = aes(x = Date, y = AirTemperature_C, group = 1)) + 
         geom_line(color = "red", lwd = 2, alpha = .5) + theme_classic()

#GLM abundance
#no separation by leaf type
effectivedate_glm <- glm(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant, family = "poisson")
dispersiontest(effectivedate_glm) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
effectivedate_glm <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(effectivedate_glm)
View(anova(effectivedate_glm))
plot(effectivedate_glm)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant$Date <- factor(metadata_plant$Date, levels = c("3/12/2021", "3/26/2021", "4/9/2021", "4/23/2021", "5/7/2021", "5/22/2021", "6/3/2021", "6/18/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "8/27/2021", "9/10/2021", "9/24/2021", "10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
effectivedate_glm <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant)
summary(effectivedate_glm)
effectivedate_summ <- data.frame(summary(effectivedate_glm)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(effectivedate_summ, effectivedate_summ$Pr...z.. > 5e-2))

#abundance over time ephemeral leaves only
metadata_plant_ephemeral <- subset(metadata_plant, metadata_plant$Type == "Ephemeral")
effectivedate_glm_ephemeral <- glm(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral, family = "poisson")
dispersiontest(effectivedate_glm_ephemeral) #w/ poisson, data is over dispersed (>1)
#since the model is over dispersed with the poisson approach, use a negative binomial distribution
effectivedate_glm_ephemeral <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(effectivedate_glm_ephemeral)
plot(effectivedate_glm_ephemeral)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_ephemeral$Date <- factor(metadata_plant_ephemeral$Date, levels = c("4/23/2021", "5/7/2021", "5/22/2021",  "6/3/2021", "6/18/2021", "7/2/2021", "7/16/2021", "7/30/2021"))
effectivedate_glm_ephemeral <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_ephemeral)
summary(effectivedate_glm_ephemeral)
effectivedate_summ <- data.frame(summary(effectivedate_glm_ephemeral)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(effectivedate_summ, effectivedate_summ$Pr...z.. > .05))

#abundance over time persistent leaves only
metadata_plant_persistent <- subset(metadata_plant, metadata_plant$Type == "Persistent")
effectivedate_glm_persistent <- glm(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent, family = "poisson")
dispersiontest(effectivedate_glm_persistent) #w/ poisson, data is NOT over dispersed (>1)
#since the model is not over dispersed, use poisson
#effectivedate_glm_persistent <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(effectivedate_glm_persistent)
plot(effectivedate_glm_persistent)

#explore which categories (dates) are similar to each other by changing the first date
metadata_plant_persistent$Date <- factor(metadata_plant_persistent$Date, levels = c("3/12/2021", "3/26/2021", "4/9/2021", "4/23/2021", "7/2/2021", "7/16/2021", "7/30/2021", "8/13/2021", "8/27/2021", "9/10/2021", "9/24/2021","10/8/2021", "10/22/2021", "11/5/2021", "11/19/2021", "12/3/2021", "12/17/2021", "12/31/2021", "1/14/2022",  "1/28/2022", "2/10/2022", "2/24/2022"))
effectivedate_glm_persistent <- glm.nb(effective.ITS ~ Date + (1|Plant_num), data = metadata_plant_persistent)
summary(effectivedate_glm_persistent)
effectivedate_summ <- data.frame(summary(effectivedate_glm_persistent)[["coefficients"]])
#Which dates are non-sginificant?
row.names(subset(effectivedate_summ, effectivedate_summ$Pr...z.. > .05))

#Combined models ####
#Richness model
#Load from file
richness_bayes <- readRDS("richness_bayes.rds")
#Richness_bayes <- brm(richness.ITS ~  ~ scale(LeafAge) + scale(I(LeafAge^2)) + Type + scale(AirTemperature_C) + scale(I(AirTemperature_C^2)) + scale(Precipitation_mm) + scale(WindSpeed_m.s) + scale(Delta15N) + scale(Delta13C) + (1 | Plant), data = metadata_plant)
summary(richness_bayes)
pp_check(richness_bayes)
pp_check(richness_bayes, type="stat", stat="mean")
mcmc_plot(richness_bayes, variable = c("^b_", "^sd_"), regex=TRUE)
parameters = c("b_scaleLeafAge", "b_scaleILeafAgeE2", "b_TypePersistent", "b_scaleAirTemperature_C", "b_scaleIAirTemperature_CE2", "b_scalePrecipitation_mm", "b_scaleWindSpeed_m.s", "b_scaleDelta15N", "b_scaleDelta13C", "sd_Plant__Intercept")
mcmc_areas(richness_bayes, pars = parameters,
           prob = 0.9, # 90% intervals 
           prob_outer = 0.99, # 99%
           point_est = "mean") +
  theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")+ theme(text = element_text(size = 30))
#marginal effects plots
plot(marginal_effects(richness_bayes))

#Abundance model
#read from file 
abundance_bayes <- readRDS("abundance_bayes.rds")
#abundance_bayes <- brm(abundance.ITS ~  ~ scale(LeafAge) + scale(I(LeafAge^2)) + Type + scale(AirTemperature_C) + scale(I(AirTemperature_C^2)) + scale(Precipitation_mm) + scale(WindSpeed_m.s) + scale(Delta15N) + scale(Delta13C) + (1 | Plant), data = metadata_plant)
summary(abundance_bayes)
plot(marginal_effects(abundance_bayes))
marg_abund <- marginal_effects(abundance_bayes)
#probability of direction
posterior_abund <- as.data.frame(abundance_bayes) # extract posterior samples from model fit
head(posterior_abund) # This will show all the posteriors for all variables
hist(posterior_abund$b_TypePersistent)
length(which(posterior_abund$b_TypePersistent < 0)) / nrow(posterior_abund) # swap out > for < where needed.

#Effective Number of Species model
#read from file 
effective_bayes <- readRDS("effective_bayes.rds")
#effective_bayes <- brm(effective.ITS ~ scale(LeafAge) + scale(I(LeafAge^2)) + Type + scale(AirTemperature_C) + scale(I(AirTemperature_C^2)) + scale(Precipitation_mm) + scale(WindSpeed_m.s) + scale(Delta15N) + scale(Delta13C) + (1 | Plant), data = metadata_plant)
summary(effective_bayes)
plot(marginal_effects(effective_bayes))


#plot posteriors for all models
par(mfrow = c(length(parameters),1), mai = c(0,1,0,1), mar=c(0,13,0,1))
#layout(matrix(rep(x = 1:10, times = 3), ncol = 2), mai = c(0,1,0,1), mar=c(0,1,0,1))
for (a in 1:length(parameters)){
  message(sprintf("Processing No. %s of %s",a,length(parameters)))
  if (a < length(parameters))
  {
    dens1 <- density(as_draws_array(richness_bayes, variable = parameters[a]))
    print(mean(dens1$x))
    dens2 <- density(as_draws_array(abundance_bayes, variable = parameters[a]))
    dens3 <- density(as_draws_array(effective_bayes, variable = parameters[a]))
    plot(dens1, col = "#CC6677",
         xlim = c(-2, 3), xaxt = "n", xlab="", yaxt = "n", ylab = "", main="", frame.plot = FALSE)
    polygon(dens1, col = alpha("#CC6677", .25))
    xval <- abs(dens1$x-mean(dens1$x))
    segments(x0 = mean(dens1$x), y0 = 0, x1 = mean(dens1$x), y1 = dens1$y[which(xval == min(xval))], col = "#CC6677", lwd = 3)
    mtext(strsplit(parameters[a], "_")[[1]][2], side = 2, las = 1)
    lines(dens2, col = "#88CCEE")
    polygon(dens2, col = alpha("#88CCEE", .25))
    xval <- abs(dens2$x-mean(dens2$x))
    segments(x0 = mean(dens2$x), y0 = 0, x1 = mean(dens2$x), y1 = dens2$y[which(xval == min(xval))], col = "#88CCEE", lwd = 3)
    lines(dens3, col = "#DDCC77")
    polygon(dens3, col = alpha("#DDCC77", .25))
    xval <- abs(dens3$x-mean(dens3$x))
    segments(x0 = mean(dens3$x), y0 = 0, x1 = mean(dens3$x), y1 = dens3$y[which(xval == min(xval))], col = "#DDCC77", lwd = 3)
    abline(v = 0, lty = 2, lwd = 1)
    abline(v = -3, lty = 1, lwd = 1)
  }else if (a == length(parameters)){
    dens1 <- density(as_draws_array(richness_bayes, variable = parameters[a]))
    print(mean(dens1$x))
    dens2 <- density(as_draws_array(abundance_bayes, variable = parameters[a]))
    dens3 <- density(as_draws_array(effective_bayes, variable = parameters[a]))
    par(mar=c(3,13,0,1))
    plot(dens1, col = "#CC6677",
         xlim = c(-2, 3), xlab="", yaxt = "n", ylab = "", main="", frame.plot = FALSE)
    polygon(dens1, col = alpha("#CC6677", .25))
    xval <- abs(dens1$x-mean(dens1$x))
    segments(x0 = mean(dens1$x), y0 = 0, x1 = mean(dens1$x), y1 = dens1$y[which(xval == min(xval))], col = "#CC6677", lwd = 3)
    mtext(strsplit(parameters[a], "_")[[1]][2], side = 2, las = 1)
    lines(dens2, col = "#88CCEE")
    polygon(dens2, col = alpha("#88CCEE", .25))
    xval <- abs(dens2$x-mean(dens2$x))
    segments(x0 = mean(dens2$x), y0 = 0, x1 = mean(dens2$x), y1 = dens2$y[which(xval == min(xval))], col = "#88CCEE", lwd = 3)
    lines(dens3, col = "#DDCC77")
    polygon(dens3, col = alpha("#DDCC77", .25))
    xval <- abs(dens3$x-mean(dens3$x))
    segments(x0 = mean(dens3$x), y0 = 0, x1 = mean(dens3$x), y1 = dens3$y[which(xval == min(xval))], col = "#DDCC77", lwd = 3)
    abline(v = 0, lty = 2, lwd = 1)
    abline(v = -3, lty = 1, lwd = 1)
  }
}

#plot intercepts for all models
dens1 <- density(as_draws_array(richness_bayes, variable = "b_Intercept"))
dens2 <- density(as_draws_array(abundance_bayes, variable = "b_Intercept"))
dens3 <- density(as_draws_array(effective_bayes, variable = "b_Intercept"))
plot(dens1, col = "#CC6677",
     xlim = c(0,20), xlab="", yaxt = "n", ylab = "", main="", frame.plot = FALSE)
polygon(dens1, col = alpha("#CC6677", .25))
xval <- abs(dens1$x-mean(dens1$x))
segments(x0 = mean(dens1$x), y0 = 0, x1 = mean(dens1$x), y1 = dens1$y[which(xval == min(xval))], col = "#CC6677", lwd = 3)
mtext("Intercepts", side = 2, las = 1)
lines(dens2, col = "#88CCEE")
polygon(dens2, col = alpha("#88CCEE", .25))
xval <- abs(dens2$x-mean(dens2$x))
segments(x0 = mean(dens2$x), y0 = 0, x1 = mean(dens2$x), y1 = dens2$y[which(xval == min(xval))], col = "#88CCEE", lwd = 3)
lines(dens3, col = "#DDCC77")
polygon(dens3, col = alpha("#DDCC77", .25))
xval <- abs(dens3$x-mean(dens3$x))
segments(x0 = mean(dens3$x), y0 = 0, x1 = mean(dens3$x), y1 = dens3$y[which(xval == min(xval))], col = "#DDCC77", lwd = 3)
abline(v = 0, lty = 1, lwd = 1)

#plot marginal effects
plot(conditional_effects(abundance_bayes))[[1]]
rich_eff <- conditional_effects(richness_bayes)
abund_eff <-conditional_effects(abundance_bayes)
effective_eff <- conditional_effects(effective_bayes)

plot(estimate__ ~ LeafAge, data = rich_eff$LeafAge, type = "l", col = "#CC6677", lwd = 2, ylim = c(0,100), ylab = "Richness")
polygon(x = c(rich_eff$LeafAge$LeafAge, rev(rich_eff$LeafAge$LeafAge)),
        y = c(rich_eff$LeafAge$upper__, rich_eff$LeafAge$lower__),
        col =  adjustcolor("#CC6677", alpha.f = 0.10), border = NA)

plot(estimate__ ~ LeafAge, data = effective_eff$LeafAge, type = "l", col = "#DDCC77", lwd = 2, yaxt = "n")
polygon(x = c(effective_eff$LeafAge$LeafAge, rev(effective_eff$LeafAge$LeafAge)),
        y = c(effective_eff$LeafAge$upper__, effective_eff$LeafAge$lower__),
        col =  adjustcolor("#DDCC77", alpha.f = 0.10), border = NA)


#response prediction
print(ggpredict(model = richness_bayes, terms = "LeafAge [0:21 by=1]", type="fixed", ci_lvl = .95, re.formula=NA), n = Inf)

posterior_predict(object = richness_bayes, re_formula = "LeafAge")

#generate data table where all variables are held at mean except leaf age which varies from 0 to 21
LeafAge_varied <- new_data(model = richness_bayes, terms = "LeafAge [0:21] by=1")

LeafAge_varied <- richness_bayes %>%
  data_grid(LeafAge = seq_range(LeafAge,n=210),
            Type = mean(Type),
            AirTemperature_C = mean(AirTemperature_C),
            Precipitation_mm = mean(Precipitation_mm),
            WindSpeed_m.s = mean(WindSpeed_m.s),
            Delta15N = mean(Delta15N),
            Delta13C = mean(Delta13C))

colnames(richness_bayes$data)


rich_bayes_effects <- conditional_effects(richness_bayes) 
#maximum - rich
max(rich_bayes_effects$AirTemperature_C$estimate__) #check max estimate value (77.8359)
rich_bayes_effects$AirTemperature_C[which(rich_bayes_effects$AirTemperature_C$estimate__ == max(rich_bayes_effects$AirTemperature_C$estimate__)),]$AirTemperature_C #get value for leaf age at max estimate (9.969697)
#minimum - rich
min(rich_bayes_effects$AirTemperature_C$estimate__) #check min estimate value
rich_bayes_effects$AirTemperature_C[which(rich_bayes_effects$AirTemperature_C$estimate__ == min(rich_bayes_effects$AirTemperature_C$estimate__)),]$AirTemperature_C #get value for leaf age at min estimate

effective_bayes_effects <- conditional_effects(effective_bayes)
#maximum - ens
max(effective_bayes_effects$AirTemperature_C$estimate__) #check max estimate value 
effective_bayes_effects$AirTemperature[which(effective_bayes_effects$AirTemperature_C$estimate__ == max(effective_bayes_effects$AirTemperature_C$estimate__)),]$AirTemperature_C #get value for leaf age at max estimate 
#minimum - ens
min(effective_bayes_effects$AirTemperature_C$estimate__) #check min estimate value
effective_bayes_effects$AirTemperature_C[which(effective_bayes_effects$AirTemperature_C$estimate__ == min(effective_bayes_effects$AirTemperature_C$estimate__)),]$AirTemperature_C #get value for leaf age at min estimate

abundance_bayes_effects <- conditional_effects(abundance_bayes)
#maximum - ens
max(abundance_bayes_effects$LeafAge$estimate__) #check max estimate value 
abundance_bayes_effects$LeafAge[which(abundance_bayes_effects$LeafAge$estimate__ == max(abundance_bayes_effects$LeafAge$estimate__)),]$LeafAge #get value for leaf age at max estimate 
#minimum - ens
min(abundance_bayes_effects$LeafAge$estimate__) #check min estimate value
abundance_bayes_effects$LeafAge[which(abundance_bayes_effects$LeafAge$estimate__ == min(abundance_bayes_effects$LeafAge$estimate__)),]$LeafAge #get value for leaf age at min estimate



plot(test$LeafAge$estimate__ ~ test$LeafAge$LeafAge, type = "l", lwd = 2, col = "red")
polygon(x = c(rich_bayes_effects$LeafAge$LeafAge, rev(rich_bayes_effects$LeafAge$LeafAge)),
        y = c(test$LeafAge$estimate__ - test$LeafAge$se__, rev(test$LeafAge$estimate__ + test$LeafAge$se__)),
        col = adjustcolor("red", alpha.f=0.1), border = NA)


#draw from model based on new data set
LeafAge_effect <- add_epred_draws(richness_bayes, 
                                  newdata = LeafAge_varied,
                                  ndraws = 1000,
                                  re_formula = NA)

#calculate mean of predicted response values from draws
tapply(LeafAge_effect$.epred, LeafAge_effect$LeafAge, mean)


brms::marginal_effects(richness_bayes)


#Beta Diversity Analysis ####
#Calculate bray-curtis
dist.ITS <- vegdist(data_ASV_abs, method = "bray", na.rm = TRUE)
#Run dbRDA beta diversity model
dbrda_all <- dbrda(dist.ITS ~ scale(LeafAge) + scale(I(LeafAge^2)) + Type + scale(AirTemperature_C) + scale(I(AirTemperature_C^2)) + scale(Precipitation_mm) + scale(WindSpeed_m.s) + scale(Delta15N) + scale(Delta13C), data = metadata_plant, Condition = metadata_plant$Plant, na.action = na.exclude, dist = "bray")
#extract important dbRDA things for figures
sum_all <- summary(dbrda_all)
biplot <- sum_all$biplot
sites_all <- as.data.frame(sum_all$sites)
#Create vector specifying pch by leaf type (for figure)
pch_type <- as.numeric(lapply(substr(row.names(sites_all),5,5), function(x) {
  y <- gsub("P", "16", x)
  gsub("E", "17", y)
}))
#Create gradient of leaf age colors
leafcol <- colorRampPalette(c("darkseagreen1", "darkgreen"))(length(unique(metadata_plant$LeafAge)))
#Function to plot color gradient for legend
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
#create figure of color gradient (add as legend in illustrator)
color.bar(leafcol, min = 0, max = 22, nticks = 3)
#Plot dbRDA
plot(dbRDA2 ~ dbRDA1, data = sites_all, pch = pch_type, type = "p", cex = 1.5, col = leafcol[as.factor(metadata_plant$LeafAge)])
arrows(x0 = rep(0,5), y0 = rep(0,5), x1 = dbrda_all$CCA$biplot[c(1:2,4:5,9),1]*3, y1 = dbrda_all$CCA$biplot[c(1:2,4:5,9),2]*3, lwd = 5, col = adjustcolor("black", alpha = .5))
text(dbrda_all$CCA$biplot[1,1]*3,dbrda_all$CCA$biplot[1,2]*3, "Leaf Age", pos = 3, offset = 1.3)
text(dbrda_all$CCA$biplot[2,1]*3 - .2,dbrda_all$CCA$biplot[2,2]*3 + .05, "Leaf Age^2", pos = 3)
text(dbrda_all$CCA$biplot[4,1]*3,dbrda_all$CCA$biplot[4,2]*3, "Air \nTemp. (C)", pos = 4)
text(dbrda_all$CCA$biplot[5,1]*3,dbrda_all$CCA$biplot[5,2]*3, "Air \nTemp.^2 (C)", pos = 4)
text(dbrda_all$CCA$biplot[9,1]*3,dbrda_all$CCA$biplot[9,2]*3, "δ13C", pos = 1)
legend(x = "topleft", legend= c("ephemeral","persistent"), pch = c(17,16), pt.cex = 1.5)

#Run permutational anova on dbRDA model
permanova_all <- anova.cca(dbrda_all, by = "margin", permutations = 999)
permanova_all

#ANCOM-BC ####
#ANCOM-BC at genus level
results_phyloseq_genus2 <- ancombc2(data = data_phyloseq_all, tax_level = "genus", fix_formula = "LeafAge + Type + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C", p_adj_method = "holm", prv_cut = 0.10, group = "Plant", struc_zero = TRUE, neg_lb = FALSE, alpha = 0.05)

View(results_phyloseq_genus2$res)

results_phyloseq_genus_res <- subset(results_phyloseq_genus2$res,  'diff_(Intercept)' == T | diff_LeafAge == T | diff_AirTemperature_C == T | diff_TypeEphemeral == T)
nrow(results_phyloseq_genus_res)
View(results_phyloseq_genus_res)

for (a in 1:length(results_phyloseq_genus_res$taxon)){
  if (results_phyloseq_genus_res$diff_LeafAge[a] == TRUE){
counts <- subset(data_tax, data_tax$genus == strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5])
counts <- data_ASV_abs[colnames(data_ASV_abs) %in% row.names(counts)]
counts <- cbind("LeafAge" = metadata_plant$LeafAge, log(data.frame(rowSums(counts))))
counts[counts == -Inf] <- 0
plot(rowSums.counts. ~ LeafAge, data = counts, 
        main = paste(strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5], "sp.", sep = " "), ylab = "log abundance")
abline(lm(rowSums.counts. ~ LeafAge, data = counts), col = "red")
  }
  else if (results_phyloseq_genus_res$diff_TypeEphemeral[a] == TRUE){
    counts <- subset(data_tax, data_tax$genus == strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5])
    counts <- data_ASV_abs[colnames(data_ASV_abs) %in% row.names(counts)]
    counts <- cbind("Type" = metadata_plant$Type, log(data.frame(rowSums(counts))))
    counts[counts == -Inf] <- 0
    boxplot(rowSums.counts. ~ Type, data = counts, 
            main = paste(strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5], "sp.", sep = " "), ylab = "log abundance")
  }
  else if (results_phyloseq_genus_res$diff_AirTemperature_C[a] == TRUE){
    counts <- subset(data_tax, data_tax$genus == strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5])
    counts <- data_ASV_abs[colnames(data_ASV_abs) %in% row.names(counts)]
    counts <- cbind("Temperature_C" = metadata_plant$AirTemperature_C, log(data.frame(rowSums(counts))))
    counts[counts == -Inf] <- 0
    plot(rowSums.counts. ~ Temperature_C, data = counts, 
            main = paste(strsplit(results_phyloseq_genus_res$taxon[a],"_")[[1]][5], "sp.", sep = " "), ylab = "log abundance")
    abline(lm(rowSums.counts. ~ Temperature_C, data = counts), col = "red")
  }
}

#Plot log fold changes
ANCOMBC_plot_df <- data.frame(rbind(
  cbind(strsplit(results_phyloseq_genus_res[1,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[1,]$lfc_AirTemperature_C, results_phyloseq_genus_res[1,]$se_AirTemperature_C, "Air_Temperature"),
  cbind(strsplit(results_phyloseq_genus_res[1,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[1,]$lfc_Delta13C, results_phyloseq_genus_res[1,]$se_Delta13C, "Delta_13C"),
  cbind(strsplit(results_phyloseq_genus_res[2,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[2,]$lfc_AirTemperature_C, results_phyloseq_genus_res[2,]$se_AirTemperature_C, "Air_Temperature"),
  cbind(strsplit(results_phyloseq_genus_res[3,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[3,]$lfc_AirTemperature_C, results_phyloseq_genus_res[3,]$se_AirTemperature_C, "Air_Temperature"),
  cbind(strsplit(results_phyloseq_genus_res[4,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[4,]$lfc_AirTemperature_C, results_phyloseq_genus_res[4,]$se_AirTemperature_C, "Air_Temperature"),
  cbind(strsplit(results_phyloseq_genus_res[4,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[4,]$lfc_LeafAge, results_phyloseq_genus_res[4,]$se_LeafAge, "Leaf_Age"),
  cbind(strsplit(results_phyloseq_genus_res[5,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[5,]$lfc_TypeEphemeral, results_phyloseq_genus_res[5,]$se_TypeEphemeral, "Leaf_Type"),
  cbind(strsplit(results_phyloseq_genus_res[6,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[6,]$lfc_AirTemperature_C, results_phyloseq_genus_res[6,]$se_AirTemperature_C, "Air_Temperature"),
  cbind(strsplit(results_phyloseq_genus_res[7,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[7,]$lfc_TypeEphemeral, results_phyloseq_genus_res[7,]$se_TypeEphemeral, "Leaf_Type"),
  cbind(strsplit(results_phyloseq_genus_res[8,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[8,]$lfc_TypeEphemeral, results_phyloseq_genus_res[8,]$se_TypeEphemeral, "Leaf_Type"),
  cbind(strsplit(results_phyloseq_genus_res[9,]$taxon, "_")[[1]][5], results_phyloseq_genus_res[9,]$lfc_TypeEphemeral, results_phyloseq_genus_res[9,]$se_TypeEphemeral, "Leaf_Type")))
colnames(ANCOMBC_plot_df) <- c("species", "lfc", "se", "variable")
ANCOMBC_plot_df <- ANCOMBC_plot_df[order(ANCOMBC_plot_df$variable),]

c("cadetblue4","royalblue3","darkblue","tomato1","dodgerblue2",
  "cyan","darkred","purple","mediumblue","palegoldenrod",
  "indianred","yellow","purple4","darkgreen","lightgoldenrod",
  "lightsalmon","yellow3","purple2","lightblue","firebrick",
  "navy","red4","red","darkmagenta","mediumvioletred",
  "violetred2","skyblue","dodgerblue4")

##LFC chart
par(mar=c(2,1,1,1))
layout(matrix(c(1,1,1,1,1,2,3,4,4,4,4), nrow = 11, ncol = 1, byrow = TRUE))
#air temp
dotchart(as.numeric(ANCOMBC_plot_df$lfc)[1:5], col = c("cadetblue4","black","dodgerblue2","purple2","darkred") , lcolor = "white", pch = 16, pt.cex = 1.25, xlim = c(-.1, .15))
text(ANCOMBC_plot_df$species[1:5], x = as.numeric(ANCOMBC_plot_df$lfc)[1:5], y = c(1.25, 2.25, 3.25, 4.25, 5.25), col = c("cadetblue4","black","dodgerblue2","purple2","darkred"))
legend("topleft", legend = "Air Temperature", lty = 0, bty = "n")
abline(v = 0, lty = 2, lwd = 2)
arrows(x0 = as.numeric(ANCOMBC_plot_df$lfc)[1:5] + as.numeric(ANCOMBC_plot_df$se)[1:5], y0 = 1:5, 
       x1 = as.numeric(ANCOMBC_plot_df$lfc)[1:5] - as.numeric(ANCOMBC_plot_df$se)[1:5], y1 = 1:5,
       length=0.05, angle=90, code=3, lwd = 2, col = c("cadetblue4","black","dodgerblue2","purple2","darkred"))
#delta 13c
dotchart(as.numeric(ANCOMBC_plot_df$lfc)[6], lcolor = "white", col = "cadetblue4", pch = 16, pt.cex = 1.25, xlim = c(-0.8,0), ylim = 1.5)
text(ANCOMBC_plot_df$species[6], x = -0.39, y = 1, pos = 4, col = "cadetblue4",)
text("Delta 13C", x = -.01, y = 1, lty = 0, pos = 2, bty = "n")
abline(v = 0, lty = 2, lwd = 2)
arrows(x0 = as.numeric(ANCOMBC_plot_df$lfc)[6] + as.numeric(ANCOMBC_plot_df$se)[6], y0 = 1, 
       x1 = as.numeric(ANCOMBC_plot_df$lfc)[6] - as.numeric(ANCOMBC_plot_df$se)[6], y1 = 1,
       length=0.05, angle=90, code=3, lwd = 2, col = "cadetblue4")
#Leaf Age
dotchart(as.numeric(ANCOMBC_plot_df$lfc)[7], lcolor = "white", col = "purple2", pch = 16, pt.cex = 1.25, xlim = c(-0.2,0))
text(ANCOMBC_plot_df$species[7], x = -0.11, y = 1, pos = 4, col = "purple2")
text("Leaf Age", x = -.01, y = 1, lty = 0, pos = 2, bty = "n")
abline(v = 0, lty = 2, lwd = 2)
arrows(x0 = as.numeric(ANCOMBC_plot_df$lfc)[7] + as.numeric(ANCOMBC_plot_df$se)[7], y0 = 1, 
       x1 = as.numeric(ANCOMBC_plot_df$lfc)[7] - as.numeric(ANCOMBC_plot_df$se)[7], y1 = 1,
       length=0.05, angle=90, code=3, lwd = 2, col = "purple2")
#Leaf Type
dotchart(as.numeric(ANCOMBC_plot_df$lfc)[8:11], lcolor = "white", pch = 16, pt.cex = 1.25, xlim = c(-2,3))
text(ANCOMBC_plot_df$species[8:11], x = as.numeric(ANCOMBC_plot_df$lfc)[8:11], y = c(1.25, 2.25, 3.25, 4.25))
legend("topleft", legend = "Leaf Type", lty = 0, bty = "n")
abline(v = 0, lty = 2, lwd = 2)
arrows(x0 = as.numeric(ANCOMBC_plot_df$lfc)[8:11] + as.numeric(ANCOMBC_plot_df$se)[8:11], y0 = 1:4, 
       x1 = as.numeric(ANCOMBC_plot_df$lfc)[8:11] - as.numeric(ANCOMBC_plot_df$se)[8:11], y1 = 1:4,
       length=0.05, angle=90, code=3, lwd = 2)
