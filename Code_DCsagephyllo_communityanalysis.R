#Dry Creek sagebrush phyllosphere project
#Analysis of the microbial communities on sagebrush leaves
#Authors: Jacob Heil, Leonora Bittleston

##Working Directory 
setwd("~/Sagebrush_Results/Ch1/Code/DC_Sagebrush")

##Install and load packages
if (!require("vegan")) {install.packages("vegan"); require("vegan")}
if (!require("effects")) {install.packages("effects"); require("effects")}
if (!require("ggplot2")) {install.packages("ggplot2"); require("ggplot2")}
if (!require("effects")) {install.packages("effects"); require("effects")}
if (!require("BiocManager")) {install.packages("BiocManager"); require("BiocManager")}
#09-13-2022, phyloseq not working on R 4.2.1 on Windows
if (!require("phyloseq")) {BiocManager::install("phyloseq"); require("phyloseq")}

##Data
#Community Data 
data_ASV <- read.csv("~/Sagebrush_Results/Ch1/Data_Community/Data_DCsp_ASV.csv", header = T, row.names = 1, check.names = F)
#Metadata
metadata_plant <- read.csv("~/Sagebrush_Results/Ch1/Data_Meta/Data_DCsp_meta_plant.csv", header = T, row.names = 1, check.names = F)
#############INSERT TAXONOMY AND ANY OTHER DATA HERE
#Remove outliers
metadata_plant <- metadata_plant[! row.names(metadata_plant) %in% c("DC20PF", "DC35PD"),]
#Clean data
data_ASV <- subset(data_ASV, row.names(data_ASV) %in% row.names(metadata_plant))
#The data sets should be 178 observations

##Naming conventions
#Assign new names to ASVs
data_ASV <- data_ASV[,reorder(colnames(data_ASV), order(colnames(data_ASV)))] 
#put ASVs in numeric order
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
metadata_plant$R_offset <- as.numeric(c(rep("NA",6), metadata_plant$richness.ITS[1:172]))
metadata_plant$A_offset <- as.numeric(c(rep("NA",6),metadata_plant$abundance.ITS[1:172]))
metadata_plant$S_offset <- as.numeric(c(rep("NA",6),metadata_plant$shannon.ITS[1:172]))

#compare culture and sequence data
cor.test(metadata_plant$culture_richness, metadata_plant$richness.ITS, method = "pearson")
plot(culture_richness ~ richness.ITS, data = metadata_plant)
abline(lm(culture_richness ~ richness.ITS, data = metadata_plant), col = "red")

#GLM Richness
boxplot(richness.ITS ~ Date, data = metadata_plant)
richdate_glm <- glm(richness.ITS ~ Date, data = metadata_plant)
summary(richdate_glm)

boxplot(richness.ITS ~ Season, data = metadata_plant)
richseason_glm <- glm(richness.ITS ~ Season, data = metadata_plant)
summary(richseason_glm)

boxplot(richness.ITS ~ Plant, data = metadata_plant)
richplant_glm <- glm(richness.ITS ~ Plant, data = metadata_plant)
summary(richplant_glm) #plant not significant

Richness_glm <- glm(richness.ITS ~ R_offset + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C, data = metadata_plant, family = "poisson")
summary(Richness_glm)
plot(allEffects(Richness_glm))

#GLM diversity
boxplot(shannon.ITS ~ Date, data = metadata_plant)
divdate_glm <- glm(shannon.ITS ~ Date, data = metadata_plant)
summary(divdate_glm)

boxplot(shannon.ITS ~ Season, data = metadata_plant)
divseason_glm <- glm(shannon.ITS ~ Season, data = metadata_plant)
summary(divseason_glm)

boxplot(shannon.ITS ~ Plant, data = metadata_plant)
divplant_glm <- glm(shannon.ITS ~ Plant, data = metadata_plant)
summary(divplant_glm) #not significant

Shannon_glm <- glm(shannon.ITS ~ R_offset + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C, data = metadata_plant)
summary(Shannon_glm)
plot(allEffects(Shannon_glm))

##### Beta Diversity Analyses
#Calculate bray-curtis
dist.ITS <- vegdist(data_ASV.r, method = "bray")

#Beta dispersion
betadisp_season <- betadisper(dist.ITS, metadata_plant$Season, type = "centroid")
boxplot(betadisp_season, xlab = "")
## Permutation test for pairwise comparisons
permutest(betadisp_season, pairwise = TRUE, permutations = 99)
anova(betadisp_season) #significant

betadisp_date <- betadisper(dist.ITS, metadata_plant$Date, type = "centroid")
boxplot(betadisp_date, xlab = "")
permutest(betadisp_date, pairwise = TRUE, permutations = 99)
anova(betadisp_date) #significant

betadisp_plant <- betadisper(dist.ITS, metadata_plant$Plant, type = "centroid")
boxplot(betadisp_plant)
permutest(betadisp_plant, pairwise = TRUE, permutations = 99)
anova(betadisp_plant) #Plant has the highest variance

#PERMANOVA
perm_date <- adonis2(dist.ITS ~ Date + Plant, data = metadata_plant, by = "margin") 
perm_date

perm_season <- adonis2(dist.ITS ~ Season + Plant, data = metadata_plant, by = "margin") 
perm_season

#Mantel
mant <- mantel(dist.ITS, vegdist(metadata_plant[,"AirTemperature_C"], method="euclidean"), permutations=999)
mant

factorlist <- colnames(metadata_plant)[c(21:22,27)]
mantdf.ITS <- c()
for (i in factorlist){
  #Mantel
  mant <- mantel(dist.ITS, vegdist(metadata_plant[[i]], method="euclidean"), permutations=999) 
  #Generate mantel stat tables for all species entered
  mantdf.ITS[[i]] <- c("Statistic" = mant[["statistic"]],"Significance"=mant[["signif"]])
}
mantdf.ITS <- data.frame(mantdf.ITS)
mantdf.ITS
p.adjust(mantdf.ITS[2,], method = "holm")

#NMDS
set.seed(600)
nmds1 <- vegan::metaMDS(dist.ITS, k = 2, trymax=1000)

#Season
plot(nmds1$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col = c("green","orange", "brown","blue")[as.factor(metadata_plant$Season)],
     pch=19, main = "Season")
ordihull(nmds1$points, groups = metadata_plant$Season, col = c("green","orange", "brown","blue"))

#Plant
plot(nmds1$points, xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col = rainbow(6)[as.factor(metadata_plant$Plant)],
     pch=19, main = "Plant")
ordihull(nmds1$points, groups = metadata_plant$Plant, col = rainbow(6))


###Beta-diversity subsetted by season

#Spring
data_ASV_spring <- subset(data_ASV.r, metadata_plant$Season %in% "Spring")
meta_spring <- subset(metadata_plant, metadata_plant$Season %in% "Spring")

dist_spring <- vegdist(data_ASV_spring, method = "bray")

bdisp_spring <- betadisper(dist_spring, meta_spring$Plant, type = "centroid")
anova(bdisp_spring)

permutest(bdisp_spring, pairwise = TRUE, permutations = 99)

perm_spring <- adonis2(dist_spring ~ Plant, data = meta_spring, by = "margin") 
perm_spring

factorlist <- colnames(metadata_plant)[c(21:22,27)]
mantdf_spring <- c()
for (i in factorlist){
  #Mantel
  mant <- mantel(dist_spring, vegdist(meta_spring[[i]], method="euclidean"), permutations=999) 
  #Generate mantel stat tables for all species entered
  mantdf_spring[[i]] <- c("Statistic" = mant[["statistic"]],"Significance"=mant[["signif"]])
}
mantdf_spring <- data.frame(mantdf_spring)
mantdf_spring
p.adjust(mantdf_spring[2,], method = "holm")

#NMDS
set.seed(600)
nmds_spring <- vegan::metaMDS(dist_spring, k = 2, trymax=500)
ordiplot(nmds_spring, type = "text")

#Summer
data_ASV_Summer <- subset(data_ASV.r, metadata_plant$Season %in% "Summer")
meta_Summer <- subset(metadata_plant, metadata_plant$Season %in% "Summer")

dist_Summer <- vegdist(data_ASV_Summer, method = "bray")

bdisp_Summer <- betadisper(dist_Summer, meta_Summer$Plant, type = "centroid")
anova(bdisp_Summer)
boxplot(bdisp_Summer)

permutest(bdisp_Summer, pairwise = TRUE, permutations = 99)

perm_Summer <- adonis2(dist_Summer ~ Plant, data = meta_Summer, by = "margin") 
perm_Summer

factorlist <- colnames(metadata_plant)[c(21:22,27)]
mantdf_Summer <- c()
for (i in factorlist){
  #Mantel
  mant <- mantel(dist_Summer, vegdist(meta_Summer[[i]], method="euclidean"), permutations=999) 
  #Generate mantel stat tables for all species entered
  mantdf_Summer[[i]] <- c("Statistic" = mant[["statistic"]],"Significance"=mant[["signif"]])
}
mantdf_Summer <- data.frame(mantdf_Summer)
mantdf_Summer
p.adjust(mantdf_Summer[2,], method = "holm")

#NMDS
set.seed(600)
nmds_Summer <- vegan::metaMDS(dist_Summer, k = 3, trymax=1000)
plot(nmds_Summer$points[,1:2], xlab="NMDS Axis 1", ylab="NMDS Axis 2", 
     col = rainbow(6)[as.factor(meta_Summer$Plant)],
     pch=19, main = "Plant")
ordihull(nmds_Summer$points[,1:2], groups=meta_Summer$Plant, col = rainbow(6))

plot(nmds_Summer$points[,2:3], xlab="NMDS Axis 2", ylab="NMDS Axis 3", 
     col = rainbow(6)[as.factor(meta_Summer$Plant)],
     pch=19, main = "Plant")
ordihull(nmds_Summer$points[,2:3], groups=meta_Summer$Plan, col = rainbow(6))

plot(nmds_Summer$points[,c(1,3)], xlab="NMDS Axis 1", ylab="NMDS Axis 3", 
     col = rainbow(6)[as.factor(meta_Summer$Plant)],
     pch=19, main = "Plant")
ordihull(nmds_Summer$points[,c(1,3)], groups=meta_Summer$Plan, col = rainbow(6))

#Fall
data_ASV_Fall <- subset(data_ASV.r, metadata_plant$Season %in% "Fall")
meta_Fall <- subset(metadata_plant, metadata_plant$Season %in% "Fall")

dist_Fall <- vegdist(data_ASV_Fall, method = "bray")

bdisp_Fall <- betadisper(dist_Fall, meta_Fall$Plant, type = "centroid")

permutest(bdisp_Fall, pairwise = TRUE, permutations = 99)

anova(bdisp_Fall)

perm_Fall <- adonis2(dist_Fall ~ Plant, data = meta_Fall, by = "margin") 
perm_Fall

factorlist <- colnames(metadata_plant)[c(21:22,27)]
mantdf_Fall <- c()
for (i in factorlist){
  #Mantel
  mant <- mantel(dist_Fall, vegdist(meta_Fall[[i]], method="euclidean"), permutations=999) 
  #Generate mantel stat tables for all species entered
  mantdf_Fall[[i]] <- c("Statistic" = mant[["statistic"]],"Significance"=mant[["signif"]])
}
mantdf_Fall <- data.frame(mantdf_Fall)
mantdf_Fall
p.adjust(mantdf_Fall[2,], method = "holm")

#NMDS
set.seed(600)
nmds_Fall <- vegan::metaMDS(dist_Fall, k = 2, trymax=500)
ordiplot(nmds_Fall, type = "text")

#Winter
data_ASV_Winter <- subset(data_ASV.r, metadata_plant$Season %in% "Winter")
meta_Winter <- subset(metadata_plant, metadata_plant$Season %in% "Winter")

dist_Winter <- vegdist(data_ASV_Winter, method = "bray")

bdisp_Winter <- betadisper(dist_Winter, meta_Winter$Plant, type = "centroid")

permutest(bdisp_Winter, pairwise = TRUE, permutations = 99)

anova(bdisp_Winter)

perm_Winter <- adonis2(dist_Winter ~ Plant, data = meta_Winter, by = "margin") 
perm_Winter

factorlist <- colnames(metadata_plant)[c(21:22,27)]
mantdf_Winter <- c()
for (i in factorlist){
  #Mantel
  mant <- mantel(dist_Winter, vegdist(meta_Winter[[i]], method="euclidean"), permutations=999) 
  #Generate mantel stat tables for all species entered
  mantdf_Winter[[i]] <- c("Statistic" = mant[["statistic"]],"Significance"=mant[["signif"]])
}
mantdf_Winter <- data.frame(mantdf_Winter)
mantdf_Winter
p.adjust(mantdf_Winter[2,], method = "holm")

#NMDS
set.seed(600)
nmds_Winter <- vegan::metaMDS(dist_Winter, k = 2, trymax=500)
ordiplot(nmds_Winter, type = "text")

plot(nmds_Winter$points, xlab="NMDS Axis 1", ylab="NMDS Axis 3", 
     col = rainbow(6)[as.factor(meta_Winter$Plant)],
     pch=19, main = "Plant")
ordispider(nmds_Winter$points, groups=meta_Winter$Plant, col = rainbow(6))

## ANCOM ####
#Detect ASVs playing significant role in driving trends for specified factor

#Load Function 'ANCOM.main' from S. Mandal (2017)####
ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  asv.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(asv.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

#Load Function 'ANCOM.w' from S. Mandal (2017) ####
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",row.names=NULL), check.names = F)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by=otu_data$Sample.ID),row.names=NULL, check.names=F)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      datalr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = datalr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = datalr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = datalr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = datalr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  
  ### Code to extract surrogate p-value ####
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}
##### End of functions ####

#ANCOM - 16s & 18s, host species as factor#####
#Create "Sample.ID" column all data tables
#ANCOM requires that data be formatted so that first *column* is named "Sample.ID"
#Sample IDs as row names does not count!
data_ancom <- data.frame("Sample.ID" = row.names(data_ASV.r), data_ASV.r, check.names = F)
meta_ancom <- data.frame("Sample.ID" = row.names(metadata_plant), metadata_plant)

ANCOM.ITS <- ANCOM.main(data_ancom, meta_ancom, F, F, "Date", NULL, NULL, F, NULL, 2, .05, .9)

ANCOM.ITS

#Find original ASV names for significant ASVs
namingconv_ITS[which(namingconv_ITS$New == ANCOM.ITS$W.taxa[1,1]),]  
namingconv_ITS[which(namingconv_ITS$New == ANCOM.ITS$W.taxa[2,1]),]  
namingconv_ITS[which(namingconv_ITS$New == ANCOM.ITS$W.taxa[3,1]),]  

ANCOM_plots <- data.frame("Plant" = metadata_plant$Date, "ITS2838" = data_ASV[,ANCOM.ITS$W.taxa[1,1]], "ITS0060" = data_ASV[,ANCOM.ITS$W.taxa[2,1]], "ITS2967" = data_ASV[,ANCOM.ITS$W.taxa[3,1]])

boxplot(ITS2838 ~ Plant, data = ANCOM_plots)
boxplot(ITS0060~ Plant, data = ANCOM_plots)
boxplot(ITS2967 ~ Plant, data = ANCOM_plots)
