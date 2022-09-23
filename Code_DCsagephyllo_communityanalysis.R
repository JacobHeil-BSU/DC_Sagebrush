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
metadata_plant$R_offset <- as.numeric(c(rep("NA",6), metadata_plant$richness.ITS[1:172]))
metadata_plant$A_offset <- as.numeric(c(rep("NA",6),metadata_plant$abundance.ITS[1:172]))
metadata_plant$S_offset <- as.numeric(c(rep("NA",6),metadata_plant$shannon.ITS[1:172]))

#GLM Richness
boxplot(richness.ITS ~ Date, data = metadata_plant)
richdate_glm <- glm(richness.ITS ~ Date, data = metadata_plant)
summary(richdate_glm)

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

boxplot(shannon.ITS ~ Plant, data = metadata_plant)
divplant_glm <- glm(shannon.ITS ~ Plant, data = metadata_plant)
summary(divplant_glm) #not significant

Shannon_glm <- glm(shannon.ITS ~ R_offset + AirTemperature_C + Precipitation_mm + WindSpeed_m.s + Delta15N + Delta13C, data = metadata_plant)
summary(Shannon_glm)
plot(allEffects(Shannon_glm))

##### Beta Diversity Analyses
#Calculate bray-curtis
dist.ITS <- vegdist(data_ASV.r, method = "bray")

#PERMANOVA
perm <- adonis2(dist.ITS ~ Date + Plant, data = metadata_plant, by = "margin") 
perm

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
ordiplot(nmds1, type = "text")

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
