
############# Differential abundance testing

folder = "G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/"

library(phyloseq)

ps2 <- readRDS(paste(folder, "phyloseq/Phyloseq_Agglomerate_Species.rds", sep="/"))
ps2h = subset_samples(ps2, sample_data(ps2)$Group3=="Healthy")


#### Wilcoxon text
# https://microbiome.github.io/tutorials/all.html
#BiocManager::install("microbiome")
library(microbiome)
library(patchwork)

library(qqman)

#Parameters to analyse
Param <- c("Run_16ST", "I18_num_ventre", 
           colnames(sample_data(ps2h))[grep("Num", colnames(sample_data(ps2h)))], 
           "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")

#Taxa levels
cluster_taxa <- c("phylum", "class", "order", "family", "Genus", "Species")

#Loop within a loop to test each parameter in each level
for(cl_taxa in cluster_taxa){
  
  if(!cl_taxa=="Species"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
  }else{
    cluster_data <- ps2h
  }
  
  #Testing only the taxa that appears in a minimum of a 10% of samples
  Core10.taxa <- readLines(paste0(folder,"phyloseq/CoreTaxa_",cl_taxa,"_010.txt"))
  
  if(cl_taxa=="phylum"){
    Core10.taxa <- paste0("Ph_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ph_",taxa_names(cluster_data))
  }else if(cl_taxa=="class"){
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    cluster_data <- subset_taxa(cluster_data, class %in% Core10.taxa)
    Core10.taxa <- paste0("Cl_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Cl_",taxa_names(cluster_data))
  }else if(cl_taxa=="order"){
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    cluster_data <- subset_taxa(cluster_data, order %in% Core10.taxa)
    Core10.taxa <- paste0("Or_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Or_",taxa_names(cluster_data))
  }else if(cl_taxa=="family"){
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    cluster_data <- subset_taxa(cluster_data, family %in% Core10.taxa)
    Core10.taxa <- paste0("Fa_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Fa_",taxa_names(cluster_data))
  }else if(cl_taxa=="Genus"){
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    cluster_data <- subset_taxa(cluster_data, Genus %in% Core10.taxa)
    Core10.taxa <- paste0("Ge_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ge_",taxa_names(cluster_data))
  }else{
    cluster_data <- prune_taxa(Core10.taxa, cluster_data)
  }
  
  ## Healthy
  for(pa in Param){
    
    #for parameters with 2 categories: Wilcox.text
    if(!pa%in%c("Run_16ST", "Location_Num", "Smoking_Num", "I18_num_ventre", "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")){
      # Calculate p-values with the two different methods for each taxonomic unit
      pvalue.NonP <- c()
      for (taxa in Core10.taxa) {
        # Create a new data frame for each taxonomic group
        df <- data.frame(Abundance = abundances(cluster_data)[taxa,],
                         Log10_Abundance = log10(1 + abundances(cluster_data)[taxa,]),  
                         Group = meta(cluster_data)[,pa])
        
        pvalue.NonP <- c(pvalue.NonP, wilcox.test(Log10_Abundance ~ Group, data = df)$p.value)  
      }
      # Arrange the results in a data.frame
      pvalues <- data.frame(taxon = Core10.taxa,
                            pvalue.NonP = pvalue.NonP)
      
      # Note that multiple testing occurs.
      # We must correct the p-values.
      # let us apply the standard Benjamini-Hochberg False Discovery Rate (FDR)
      # correction
      pvalues$pvalue.NonP.adjusted <- p.adjust(pvalues$pvalue.NonP)
      
    #for parameters with more than 2 categories: Kruskal.Wallis test  
    }else if(pa%in%c("Run_16ST", "Location_Num", "Smoking_Num", "I18_num_ventre", "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")){
      pvalue.NonP <- c()
      for (taxa in Core10.taxa) {
        # Create a new data frame for each taxonomic group
        df <- data.frame(Abundance = abundances(cluster_data)[taxa,],
                         Log10_Abundance = log10(1 + abundances(cluster_data)[taxa,]),  
                         Group = meta(cluster_data)[,pa])
        
        pvalue.NonP <- c(pvalue.NonP, kruskal.test(Log10_Abundance ~ Group, data = df)$p.value)  
      }
      # Arrange the results in a data.frame
      pvalues <- data.frame(taxon = Core10.taxa,
                            pvalue.NonP = pvalue.NonP)
      
      # Note that multiple testing occurs.
      # We must correct the p-values.
      # let us apply the standard Benjamini-Hochberg False Discovery Rate (FDR)
      # correction
      pvalues$pvalue.NonP.adjusted <- p.adjust(pvalues$pvalue.NonP)
    }else{}
    
    # Diff expression analysis (estimate dispersion, fitting model and testing)
    pvalues <- pvalues[order(pvalues$pvalue.NonP.adjusted),] # FDR is for BH
    write.table(pvalues, paste0(folder,"modeling/02.Wilcoxon/02.Wilcoxon_Healthy_",pa,"_Pvalues_",cl_taxa,"_Core10.txt"), quote=F, sep="\t", row.names=F, col.names=T)
    
    if(sum(!pvalues$pvalue.NonP.adjusted==1, na.rm=T)==0){
      tiff(paste0(folder,"/modeling/02.Wilcoxon/02.Wilcoxon_Healthy_",pa,"_QQ_",cl_taxa,"_Core10.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(pvalues$pvalue.NonP, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      dev.off()
    }else{
      tiff(paste0(folder,"/modeling/02.Wilcoxon/02.Wilcoxon_Healthy_",pa,"_QQ_",cl_taxa,"_Core10.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(pvalues$pvalue.NonP, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      qq(pvalues$pvalue.NonP.adjusted, main="adjusted P")
      dev.off()
    }
  }
}

###### MetagenomeSeq 
# https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

Param <- c("Run_16ST", "I18_num_ventre",
           colnames(sample_data(ps2h))[grep("Num", colnames(sample_data(ps2h)))], 
           "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")
ParamCore <- c("Sex_Num", "Age_Bi_Num", "Study_Num", "Run_16ST", "Location_Num")
Param <- Param[!Param%in%ParamCore]

cluster_taxa <- c("phylum", "class", "order", "family", "Genus", "Species")

##### Only Healthy
for(cl_taxa in cluster_taxa){
  
  Core10.taxa <- readLines(paste0(folder,"phyloseq/CoreTaxa_",cl_taxa,"_010.txt"))
  
  if(cl_taxa=="phylum"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- paste0("Ph_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ph_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="class"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Cl_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Cl_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="order"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Or_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Or_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="family"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Fa_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Fa_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="Genus"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Ge_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ge_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else{
    cluster_data <- ps2h
    Counts <- data.frame(t(otu_table(cluster_data)))
  }
  
  #Counts[1:5,1:10]
  #colnames(Counts)
  
  Data <- data.frame(sample_data(cluster_data))
  Data$ID_icobiobank <- paste0("X", Data$ID_icobiobank)
  row.names(Data) <- Data$ID_icobiobank
  #names(Data)
  #Data[1:5,1:10]
  phenotypeData = AnnotatedDataFrame(Data)
  
  Taxa <- data.frame(tax_table(cluster_data))
  head(Taxa)
  OTUdata = AnnotatedDataFrame(Taxa)
  
  obj <- newMRexperiment(Counts, phenotypeData, OTUdata)
  
  #phenoData(obj)
  #head(pData(obj), 3)
  
  #featureData(obj)
  #head(fData(obj)[,2:7], 3)
  
  #head(MRcounts(obj[, 1:2]))
  
  #featuresToKeep = which(rowSums(obj) >= 100)
  #samplesToKeep = which(pData(obj)$SmokingStatus == "Smoker")
  #obj_smokers = obj[featuresToKeep, samplesToKeep]
  #obj_smokers
  
  # Alternative normalization scaling factors can be accessed or 
  # replaced with the normFactors method:
  #normFactors(obj) <- rnorm(ncol(obj))
  #head(normFactors(obj))
  #libSize(obj) <- rnorm(ncol(obj))
  #head(libSize(obj))
  
  # data can be filtered to maintain a threshold of minimum depth or OTU presence:
  # filterData(mouseData, present = 10, depth = 1000)
  
  # Two MRexperiment-class objects can be merged with the mergeMRexperiments function,
  #newobj = mergeMRexperiments(mouseData, mouseData)
  
  ### Normalization
  #cumNormMat(obj)
  #MRcounts(obj,norm=TRUE,log=FALSE)
  ## Calculating normalization factors
  #p = cumNormStat(obj)
  #obj = cumNorm(obj, p = 0.9)
  ## Exporting data
  # normalized count matrices
  #mat = MRcounts(obj, norm = TRUE, log = TRUE)[1:5, 1:5]
  #exportMat(mat, file = file.path(folder, "tmp.tsv"))
  # sample statistics
  #exportStats(obj[, 1:5], file = file.path(folder, "tmp.tsv"))
  ##
  
  ### Statistical testing
  ## Zero-infated Log-Normal mixture model for each feature
  #obj <- newMRexperiment(Counts, phenotypeData, OTUdata)
  #sum(is.na(pData(obj)$BMI_binned))
  #sum(is.na(pData(obj)$Study))
  #obj = obj[, -which(is.na(pData(obj)$BMI_binned))]
  #obj = filterData(obj, present = 1, depth = 1)
  #obj <- cumNorm(obj, p = 0.5)
  #pd <- pData(obj)
  #mod <- model.matrix(~1 + BMI_binned, data = pd)
  #objres1 = fitFeatureModel(obj, mod)
  #head(MRcoefs(objres1))
  
  ## Zero-infated Gaussian mixture model
  #controls = grep("Extraction.Control", pData(obj)$SampleType)
  #obj = obj[, -controls]
  rareFeatures = which(rowSums(MRcounts(obj) > 0) < 2)
  obj = obj[-rareFeatures, ]
  objp = cumNormStat(obj, pFlag = TRUE, main = "Trimmed data")
  obj = cumNorm(obj, p = objp)
  
  Age_Bi_Num = pData(obj)$Age_Bi_Num
  Sex_Num = pData(obj)$Sex_Num
  Study_Num = pData(obj)$Study_Num
  Run_16ST = pData(obj)$Run_16ST
  Location_Num = pData(obj)$Location_Num
  
  normFactor = normFactors(obj)
  normFactor = log2(normFactor/median(normFactor) + 1)
  
  #Modeling ParamCore
  mod = model.matrix(~Sex_Num + Age_Bi_Num + Study_Num + Run_16ST + Location_Num + normFactor)
  settings = zigControl(maxit = 10, verbose = TRUE)
  fit = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE,
               control = settings)
  # Exporting fits
  #head(MRcoefs(fit, coef = 2))
  #head(MRcoefs(fit))
  
  ## Log Normal permutation test
  coeffOfInterest = 2
  #res = fitLogNormal(obj = obj, mod = mod, useCSSoffset = FALSE,
  #                   B = 10, coef = coeffOfInterest)
  res = fitLogNormal(obj = obj, mod = mod, useCSSoffset = FALSE)
  
  # extract the absolute fold-change estimates
  Coeffs <- data.frame(res$fit$coefficients)
  #head(Coeffs)
  
  StDev_unScaled <- data.frame(res$fit$stdev.unscaled)
  #head(StDev_unScaled)
  StDev_Sex_Num <- Coeffs$Sex_Num / res$t
  ScalingFactor <- StDev_Sex_Num/StDev_unScaled$Sex_Num
  StDev <- StDev_unScaled
  for(i in 1:nrow(StDev)){
    StDev[i,] <- StDev[i,]*ScalingFactor[i]
  }
  
  Ts <- Coeffs/StDev
  #head(Ts)
  
  Ps <- Ts
  for(i in 1:nrow(Ps)){
    for(j in 1:ncol(Ps)){
      Ps[i,j] <- pt(q=abs(Ts[i,j]), df=308, lower.tail=F)*2 
    }
  }
  #head(Ps)
  # adjust for multiple testing 
  Psadj <- Ps
  for(i in 1:ncol(Psadj)){
    Psadj[,i] <- p.adjust(Psadj[,i], method = "fdr")
  }
  #head(Psadj)
  
  for(pa in ParamCore){
    # Parameters of Diff expression anaysis
    results_var <- data.frame(rownames(Coeffs),
                              Coeffs[,colnames(Coeffs)%in%pa],
                              StDev[,colnames(StDev)%in%pa],
                              Ts[,colnames(Ts)%in%pa],
                              Ps[,colnames(Ps)%in%pa],
                              Psadj[,colnames(Psadj)%in%pa])
    
    colnames(results_var) <- c("Taxon", "logFC", "StDev", "T", "P", "Padj")
    
    results_var <- results_var[order(results_var$P),]
    results_var <- results_var[order(results_var$Padj),]
    results_var <- results_var[results_var$Taxon%in%Core10.taxa,]
    write.table(results_var, paste0(folder,"modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), quote=F, sep="\t", row.names=F, col.names=T)
    
    tiff(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
    par(mfrow=c(1,2))
    qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
    qq(results_var$Padj, main="adjusted P")
    dev.off()
    
    #Volcano Plot
    tiff(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
    par(mfrow=c(1,1))
    with(results_var, plot(logFC, -log10(P), pch=20, main=paste0("Volcano plots ",pa," \nraw estimates")))
    with(subset(results_var, Padj<.01 ), points(logFC, -log10(P), pch=20, col="blue"))
    with(subset(results_var, Padj<.01 & abs(logFC)>2), points(logFC, -log10(P), pch=20, col="red"))
    dev.off()
  }
  
  for(pa in Param){
    #Modeling ParamCore
    obj_Var <- obj[, -which(is.na(pData(obj)[,pa]))]
    names(pData(obj_Var)) <- gsub(pa, "Var_Num", names(pData(obj_Var)))
    
    Age_Bi_Num = pData(obj_Var)$Age_Bi_Num
    Sex_Num = pData(obj_Var)$Sex_Num
    Study_Num = pData(obj_Var)$Study_Num
    Run_16ST = pData(obj_Var)$Run_16ST
    Location_Num = pData(obj_Var)$Location_Num
    Var_Num = pData(obj_Var)$Var_Num
    
    normFactor = normFactors(obj_Var)
    normFactor = log2(normFactor/median(normFactor) + 1)
    
    #Modeling Param
    mod = model.matrix(~Sex_Num + Age_Bi_Num + Study_Num + Run_16ST + Location_Num + Var_Num + normFactor)
    settings = zigControl(maxit = 10, verbose = TRUE)
    
    if(class(try(fitZig(obj = obj_Var, mod = mod, useCSSoffset = TRUE,
                        control = settings)))=="try-error"){
      results_var <- data.frame(matrix(rep(NA, length(Core10.taxa)*6), nrow=length(Core10.taxa), ncol=6))
      colnames(results_var) <- c("Taxon", "logFC", "StDev", "T", "P", "Padj")
      results_var$Taxon <- Core10.taxa
      write.table(results_var, paste0(folder,"modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), quote=F, sep="\t", row.names=F, col.names=T)
      
    }else{
      
      fit = fitZig(obj = obj_Var, mod = mod, useCSSoffset = TRUE,
                   control = settings)
      # Exporting fits
      #head(MRcoefs(fit, coef = 2))
      #head(MRcoefs(fit))
      
      ## Log Normal permutation test
      coeffOfInterest = 2
      #res = fitLogNormal(obj = obj, mod = mod, useCSSoffset = FALSE,
      #                   B = 10, coef = coeffOfInterest)
      res = fitLogNormal(obj = obj_Var, mod = mod, useCSSoffset = FALSE)
      
      # extract the absolute fold-change estimates
      Coeffs <- data.frame(res$fit$coefficients)
      #head(Coeffs)
      
      StDev_unScaled <- data.frame(res$fit$stdev.unscaled)
      #head(StDev_unScaled)
      StDev_Sex_Num <- Coeffs$Sex_Num / res$t
      ScalingFactor <- StDev_Sex_Num/StDev_unScaled$Sex_Num
      StDev <- StDev_unScaled
      for(i in 1:nrow(StDev)){
        StDev[i,] <- StDev[i,]*ScalingFactor[i]
      }
      
      Ts <- Coeffs/StDev
      #head(Ts)
      
      Ps <- Ts
      for(i in 1:nrow(Ps)){
        for(j in 1:ncol(Ps)){
          Ps[i,j] <- pt(q=abs(Ts[i,j]), df=(length(Var_Num)-1), lower.tail=F)*2 
        }
      }
      #head(Ps)
      # adjust for multiple testing 
      Psadj <- Ps
      for(i in 1:ncol(Psadj)){
        Psadj[,i] <- p.adjust(Psadj[,i], method = "fdr")
      }
      #head(Psadj)
      
      
      # Parameters of Diff expression anaysis
      results_var <- data.frame(rownames(Coeffs),
                                Coeffs[,colnames(Coeffs)%in%"Var_Num"],
                                StDev[,colnames(StDev)%in%"Var_Num"],
                                Ts[,colnames(Ts)%in%"Var_Num"],
                                Ps[,colnames(Ps)%in%"Var_Num"],
                                Psadj[,colnames(Psadj)%in%"Var_Num"])
      
      colnames(results_var) <- c("Taxon", "logFC", "StDev", "T", "P", "Padj")
      
      results_var <- results_var[order(results_var$P),]
      results_var <- results_var[order(results_var$Padj),]
      results_var <- results_var[results_var$Taxon%in%Core10.taxa,]
      write.table(results_var, paste0(folder,"modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), quote=F, sep="\t", row.names=F, col.names=T)
      
      tiff(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      qq(results_var$Padj, main="adjusted P")
      dev.off()
      
      #Volcano Plot
      tiff(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
      par(mfrow=c(1,1))
      with(results_var, plot(logFC, -log10(P), pch=20, main=paste0("Volcano plots ",pa," \nraw estimates")))
      with(subset(results_var, Padj<.01 ), points(logFC, -log10(P), pch=20, col="blue"))
      with(subset(results_var, Padj<.01 & abs(logFC)>2), points(logFC, -log10(P), pch=20, col="red"))
      dev.off()
      
    }
  }
}

######## ALDEx2
# http://bioconductor.org/packages/release/bioc/html/ALDEx2.html
#BiocManager::install("ALDEx2")
library(ALDEx2)

# Only Healthy
for(cl_taxa in cluster_taxa){
  
  Core10.taxa <- readLines(paste0(folder,"phyloseq/CoreTaxa_",cl_taxa,"_010.txt"))
  
  if(cl_taxa=="phylum"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- paste0("Ph_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ph_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="class"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Cl_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Cl_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="order"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Or_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Or_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="family"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Fa_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Fa_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else if(cl_taxa=="Genus"){
    cluster_data = aggregate_taxa(ps2h, cl_taxa)
    Core10.taxa <- Core10.taxa[!Core10.taxa%in%c("Unknown")]
    Core10.taxa <- paste0("Ge_", Core10.taxa)
    taxa_names(cluster_data) <- paste0("Ge_",taxa_names(cluster_data))
    
    Counts <- data.frame(otu_table(cluster_data))
  }else{
    cluster_data <- ps2h
    Counts <- data.frame(t(otu_table(cluster_data)))
  }
  
  #Counts[1:5,1:10]
  #colnames(Counts)
  
  Data <- data.frame(sample_data(cluster_data))
  Data$ID_icobiobank <- paste0("X", Data$ID_icobiobank)
  row.names(Data) <- Data$ID_icobiobank
  #names(Data)
  #Data[1:5,1:10]
  
  mm <- model.matrix(~Sex_Num + Age_Bi_Num + Study_Num + Run_16ST + Location_Num, Data)
  x <- aldex.clr(Counts, mm, mc.samples=128, denom="all")
  # checking for condition length disabled!
  # using all features for denominator
  # operating in serial mode
  # computing center with all features
  glm.test <- aldex.glm(x, mm)
  #head(glm.test)
  glm.test$Taxon <- rownames(glm.test)
  
  for(pa in ParamCore){
    
    # Parameters of Diff expression anaysis
    results_var <- glm.test[,grep(pa, colnames(glm.test))]
    
    colnames(results_var) <- c("Est", "StDev", "T", "P", "Padj")
    results_var$Taxon <- rownames(results_var)
    results_var <- results_var[order(results_var$P),]
    
    results_var <- results_var[results_var$Taxon%in%Core10.taxa,]
    write.table(results_var, paste0(folder,"modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), quote=F, sep="\t", row.names=F, col.names=T)
    
    if(sum(!results_var$Padj==1)==0){
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      dev.off()
    }else{
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      qq(results_var$Padj, main="adjusted P")
      dev.off()
    }
 
    #Volcano Plot
    tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
    par(mfrow=c(1,1))
    with(results_var, plot(Est, -log10(P), pch=20, main=paste0("Volcano plots ",pa," \nraw estimates")))
    with(subset(results_var, Padj<.01 ), points(Est, -log10(P), pch=20, col="blue"))
    with(subset(results_var, Padj<.01 & abs(Est)>2), points(Est, -log10(P), pch=20, col="red"))
    dev.off()
  }
  
  for(pa in Param){
    Data_Var <- Data[!is.na(Data[,pa]),]
    colnames(Data_Var) <- gsub(pa, "Var_Num", colnames(Data_Var))
    Counts_Var <- Counts[, colnames(Counts)%in%rownames(Data_Var)]
    
    
    mm <- model.matrix(~Sex_Num + Age_Bi_Num + Study_Num + Run_16ST + Location_Num + Var_Num, Data_Var)
    x <- aldex.clr(Counts_Var, mm, mc.samples=128, denom="all")
    # checking for condition length disabled!
    # using all features for denominator
    # operating in serial mode
    # computing center with all features
    glm.test <- aldex.glm(x, mm)
    #head(glm.test)
    glm.test$Taxon <- rownames(glm.test)
    
    # Parameters of Diff expression anaysis
    results_var <- glm.test[,grep("Var_Num", colnames(glm.test))]
    
    colnames(results_var) <- c("Est", "StDev", "T", "P", "Padj")
    results_var$Taxon <- rownames(results_var)
    results_var <- results_var[order(results_var$P),]
    
    results_var <- results_var[results_var$Taxon%in%Core10.taxa,]
    write.table(results_var, paste0(folder,"modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), quote=F, sep="\t", row.names=F, col.names=T)
    
    if(sum(!results_var$Padj==1)==0){
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      dev.off()
    }else{
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_QQ_",cl_taxa,"_Core10ASV.tiff"), units="in", width=10, height=5, res=320)
      par(mfrow=c(1,2))
      qq(results_var$P, main=paste0("QQ plots for ",pa," \n unadjusted P"))
      qq(results_var$Padj, main="adjusted P")
      dev.off()
    }
    
    #Volcano Plot
    tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
    par(mfrow=c(1,1))
    with(results_var, plot(Est, -log10(P), pch=20, main=paste0("Volcano plots ",pa," \nraw estimates")))
    with(subset(results_var, Padj<.01 ), points(Est, -log10(P), pch=20, col="blue"))
    with(subset(results_var, Padj<.01 & abs(Est)>2), points(Est, -log10(P), pch=20, col="red"))
    dev.off()
  }
}


#### Summary
library(VennDiagram)
library(phyloseq)
library(memisc)
library(microbiome)
library(patchwork)
library(qqman)
library(compositions)

folder = "G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/"

ps2 <- readRDS(paste(folder, "phyloseq/Phyloseq_Agglomerate_Species.rds", sep="/"))
summary(as.factor(sample_data(ps2)$Group3))

Param <- c("Run_16ST", "I18_num_ventre", 
           colnames(sample_data(ps2))[grep("Num", colnames(sample_data(ps2)))], 
           "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")

cluster_taxa <- c("phylum", "class", "order", "family", "Genus", "Species")
cl_taxa <- "phylum"

## Only Healthy
for(cl_taxa in cluster_taxa){
  #Identifying associated taxa
  for(pa in Param){

    Taxa <- c()
    
    Wilc <- read.delim(paste0(folder,"/modeling/02.Wilcoxon/02.Wilcoxon_Healthy_",pa,"_Pvalues_",cl_taxa,"_Core10.txt"), header=T)
    Wilc_tx <- Wilc[Wilc$pvalue.wilcoxon.adjusted<0.05, "taxon"]
    
    MgS <- read.delim(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), header=T)
    MgS_tx <- MgS[MgS$Padj<0.05, "Taxon"]
    
    ALD <- read.delim(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), header=T)
    ALD_tx <- ALD[ALD$Padj<0.05, "Taxon"]
    
    Taxa <- c(Wilc_tx, MgS_tx, ALD_tx)
    Taxa <- Taxa[!duplicated(Taxa)]
    Taxa <- Taxa[!is.na(Taxa)]
    
    #Subsetting associated taxa
    Wilc_sig <- Wilc[Wilc$taxon%in%Taxa, c("taxon", "pvalue.NonP", "pvalue.NonP.adjusted")]
    colnames(Wilc_sig) <- c("Taxon","Wilcox_P", "Wilcox_Padj")
    
    MgS_sig <- MgS[MgS$Taxon%in%Taxa, c("Taxon", "logFC", "P", "Padj")]
    colnames(MgS_sig) <- c("Taxon", paste0("MgS_",c("Est", "P", "Padj")))
    
    Test_sig <- merge(Wilc_sig, MgS_sig, by="Taxon")
    
    ALD_sig <- ALD[ALD$Taxon%in%Taxa, c("Taxon", "Est", "P", "Padj")]
    colnames(ALD_sig) <- c("Taxon", paste0("ALD_",c("Est", "P", "Padj")))
    
    Test_sig <- merge(Test_sig, ALD_sig, by="Taxon")
    
    if(nrow(Test_sig)>0){
      write.table(Test_sig, paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,".txt"), quote=F, row.names=F, col.names=T, sep="\t")
    }
  }
}

#Plotting Counts
Files <- list.files(paste0(folder,"/modeling/Summary/"))
Files <- gsub("Test_Sign_Healthy_", "", Files)
Files <- gsub(".txt", "", Files)

for(cl_taxa in cluster_taxa){
  
  for(pa in Param){
    
    if(paste0(pa,"_",cl_taxa)%in%Files){
      
      Test_sig <- read.delim(paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,".txt"), header=T)
      
      #Test_sig_TF <- data.frame(Test_sig<0.1)
      #Top <- Test_sig$Taxon[rowSums(Test_sig_TF[,c("MgS_Padj", "ALD_Padj")])>1]
      #Top <- Top[!is.na(Top)]
      Top <- Test_sig$Taxon
      
      if(cl_taxa=="phylum"){
        cluster_data = aggregate_taxa(ps2, cl_taxa)
        taxa_names(cluster_data) <- paste0("Ph_",taxa_names(cluster_data))
        
        Abund <- otu_table(cluster_data)
      }else if(cl_taxa=="class"){
        cluster_data = aggregate_taxa(ps2, cl_taxa)
        taxa_names(cluster_data) <- paste0("Cl_",taxa_names(cluster_data))
        
        Abund <- otu_table(cluster_data)
      }else if(cl_taxa=="order"){
        cluster_data = aggregate_taxa(ps2, cl_taxa)
        taxa_names(cluster_data) <- paste0("Or_",taxa_names(cluster_data))
        
        Abund <- otu_table(cluster_data)
      }else if(cl_taxa=="family"){
        cluster_data = aggregate_taxa(ps2, cl_taxa)
        taxa_names(cluster_data) <- paste0("Fa_",taxa_names(cluster_data))
        
        Abund <- otu_table(cluster_data)
      }else if(cl_taxa=="Genus"){
        cluster_data = aggregate_taxa(ps2, cl_taxa)
        taxa_names(cluster_data) <- paste0("Ge_",taxa_names(cluster_data))
        
        Abund <- otu_table(cluster_data)
      }else{
        cluster_data <- ps2
        Abund <- t(otu_table(cluster_data))
      }
      
      for(taxa in Top){
        df <- data.frame(Abundance = abundances(cluster_data)[taxa,],
                         Group = meta(cluster_data)[,pa])
        #df$Abundance[which(df$Abundance==0)] <- NA
        #for(i in 1:ncol(Abund)){
        #Abund[,i][which(Abund[,i]==0)] <- NA
        #}
        clr_data <- data.frame(clr(Abund))
        
        WP<- round(print(kruskal.test(Abundance ~ Group, data = df)$p.value), 5)
        p1 <- ggplot(df, aes(x = Group, y = Abundance)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = paste0(taxa, " Abs Abun"), 
               y = paste0("Abundance\n (read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18))
        
        # Let us add the log2(1+x) version:
        df$Log2_Abundance <- log2(1 + df$Abundance)
        WP<- round(print(kruskal.test(Log2_Abundance ~ Group, data = df)$p.value), 5)
        p2 <- ggplot(df, aes(x = Group, y = Log2_Abundance)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = "Log2 Abun", 
               y = paste0("Abundance\n (log2(1+x) read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18)) 
        
        # Let us add the CLR version:
        df$CLR <- t(clr_data[taxa,])
        WP<- round(print(kruskal.test(CLR ~ Group, data = df)$p.value), 5)
        p3 <- ggplot(df, aes(x = Group, y = CLR)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = "Centered LogRatio", 
               y = paste0("Abundance\n (CLR read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18)) 
        
        tiff(paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,"_",taxa,"_Abund_Box.tiff"), units="in", width=10, height=5, res=320)
        print(p1 + p2 + p3)
        dev.off()
        
        #Density plot
        p <- ggplot(df, aes(fill = Group, x = CLR)) +
          geom_density(aes(group = Group), alpha = 0.5)
        
        tiff(paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,"_",taxa,"_Abund_Dens.tiff"), units="in", width=10, height=5, res=320)
        print(p)
        dev.off()
        
        #NAs 
        df <- data.frame(Abundance = abundances(cluster_data)[taxa,],
                         Group = meta(cluster_data)[,pa])
        df$Abundance[which(df$Abundance==0)] <- NA
        for(i in 1:ncol(Abund)){
          Abund[,i][which(Abund[,i]==0)] <- NA
        }
        clr_datah <- data.frame(clr(Abund))
        
        WP<- round(print(kruskal.test(Abundance ~ Group, data = df)$p.value), 5)
        p1 <- ggplot(df, aes(x = Group, y = Abundance)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = paste0(taxa, " Abs Abun"), 
               y = paste0("Abundance\n (read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18))
        
        # Let us add the log2(1+x) version:
        df$Log2_Abundance <- log2(1 + df$Abundance)
        WP<- round(print(kruskal.test(Log2_Abundance ~ Group, data = df)$p.value), 5)
        p2 <- ggplot(df, aes(x = Group, y = Log2_Abundance)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = "Log2 Abun", 
               y = paste0("Abundance\n (log2(1+x) read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18))    
        
        # Let us add the CLR version:
        df$CLR <- t(clr_data[taxa,])
        WP<- round(print(kruskal.test(CLR ~ Group, data = df)$p.value), 5)
        p3 <- ggplot(df, aes(x = Group, y = CLR)) +
          geom_boxplot(aes(group = Group)) +
          labs(title = "Centered LogRatio", 
               y = paste0("Abundance\n (CLR read count)"),
               x = paste0("Wilcox P: ",WP)) + 
          theme(plot.title = element_text(size=18)) 
        
        tiff(paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,"_",taxa,"_Abund_Box_NA.tiff"), units="in", width=10, height=5, res=320)
        print(p1 + p2 + p3)
        dev.off()
        
        #Density plot
        p <- ggplot(df, aes(fill = Group, x = CLR)) +
          geom_density(aes(group = Group), alpha = 0.5)
        
        tiff(paste0(folder,"/modeling/Summary/Test_Sign_Healthy_",pa,"_",cl_taxa,"_",taxa,"_Abund_Dens_NA.tiff"), units="in", width=10, height=5, res=320)
        print(p)
        dev.off()
      }
    }
  }
}


######### VolcanoPlots
library(phyloseq)

folder = "G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/"

ps2 <- readRDS(paste(folder, "phyloseq/Phyloseq_Agglomerate_Species.rds", sep="/"))
summary(as.factor(sample_data(ps2)$Group3))

Param <- c("Run_16ST", "I18_num_ventre", 
           colnames(sample_data(ps2))[grep("Num", colnames(sample_data(ps2)))], 
           "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")
pa <- "Buckland_rMEDScore"

cluster_taxa <- c("phylum", "class", "order", "family", "Genus", "Species")
cl_taxa <- "Species"

## Only Healthy
for(cl_taxa in cluster_taxa){
  #Identifying associated taxa
  for(pa in Param){
    
    Wilc <- read.delim(paste0(folder,"/modeling/02.Wilcoxon/02.Wilcoxon_Healthy_",pa,"_Pvalues_",cl_taxa,"_Core10.txt"), header=T)
    Wilc <- Wilc[, c("taxon", "pvalue.NonP", "pvalue.NonP.adjusted")]
    colnames(Wilc) <- c("Taxon","Wilcox_P", "Wilcox_Padj")
    
    MgS <- read.delim(paste0(folder,"/modeling/08.MetagenomeSeq/08.MetagenomeSeq_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), header=T)
    MgS <- MgS[, c("Taxon", "logFC", "P", "Padj")]
    colnames(MgS) <- c("Taxon", paste0("MgS_",c("Est", "P", "Padj")))
    
    ALD <- read.delim(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_FoldChangeEstimates_",cl_taxa,"_Core10ASV.txt"), header=T)
    ALD <- ALD[, c("Taxon", "Est", "P", "Padj")]
    colnames(ALD) <- c("Taxon", paste0("ALD_",c("Est", "P", "Padj")))
    
    Test <- merge(Wilc, MgS, by="Taxon")
    Test <- merge(Test, ALD, by="Taxon")
    
    write.table(Test, paste0(folder,"/modeling/Summary/Test_Healthy_",pa,"_",cl_taxa,".txt"), quote=F, row.names=F, col.names=T, sep="\t")
   
    blues <- which(Test$ALD_P<0.01)
    if(length(blues)<1){
      #Volcano Plot
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
      par(mfrow=c(1,1))
      with(Test, plot(ALD_Est, -log10(ALD_P), pch=20, main=paste0("Volcano plots ",pa," \n",cl_taxa," estimates")))
      dev.off()
    }else{
      #Volcano Plot
      tiff(paste0(folder,"/modeling/09.ALDEx2/09.ALDEx2_Healthy_",pa,"_Volcano_",cl_taxa,"_Core10ASV.tiff"), units="in", width=5, height=5, res=320)
      par(mfrow=c(1,1))
      with(Test, plot(ALD_Est, -log10(ALD_P), pch=20, main=paste0("Volcano plots ",pa," \n",cl_taxa," estimates")))
      with(subset(Test, ALD_P<.01 ), points(ALD_Est, -log10(ALD_P), pch=20, col="blue"))
      with(subset(Test, ALD_P<.01 & abs(ALD_Est)>2), points(ALD_Est, -log10(ALD_P), pch=20, col="red"))
      text((Test$ALD_Est[blues]-(Test$ALD_Est[blues]*0.1)), -log10(Test$ALD_P)[blues]*0.95, cex = 0.8, Test$Taxon[blues])
      dev.off()
    }
  }
}
head(Test)
summary(Test)



