library(tidyverse) # for plotting and wrangling data
library(SpiecEasi) # Has sparcc and also does clr transforms
library(otuSummary)
library(reshape2) # has the melt funct  ion, which I use to wrangle data
library(psych) # for calculating regular correlations with p values
library(metan)
library(ImageGP)
pass <- function(x){x}


####################
##IF WE WANT TO COMBINE ASV, GENUS AND SPECIES NAMES
#tax_table_data <- as.data.frame(tax_table(physeq))
#tax_table_data$NewColumn4 <- with(tax_table_data, ifelse(is.na(Genus), NA, 
#ifelse(is.na(Species), paste(rownames(tax_table_data), Genus, "spp."), 
#       paste(rownames(tax_table_data), paste(Genus, Species, sep = " ")))))
###################

outcorrsparcc <- sparcc(otu_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))
view(outcorrsparcc)
rownames(outcorrsparcc$Cor) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
colnames(outcorrsparcc$Cor) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
rownames(outcorrsparcc$Cov) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
colnames(outcorrsparcc$Cov) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus

view(outcorrsparcc)
outcorrsparcc$Cor[1:5, 1:5]

plotableSparcc <- outcorrsparcc$Cor %>% reorder_cormat %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()

Sparcc_plot <- plotableSparcc %>% ggplot(aes(x = Var2, y = Var1, fill = value)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
Sparcc_plot


tp0 <- proc.time()
out2 <- sparccboot(otu_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered), R = 100)
tp1 <- proc.time()
tp1 - tp0



outP <- pval.sparccboot(out2)
data.frame(outP$cors, outP$pvals) %>% head



cors <- outP$cors
pvals <- outP$pvals
sparCCpcors <- diag(0.5, nrow = dim(outcorrsparcc$Cor)[1], ncol = dim(outcorrsparcc$Cor)[1])
sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
sparCCpcors <- sparCCpcors + t(sparCCpcors)

sparCCpval <- diag(0.5, nrow = dim(outcorrsparcc$Cor)[1], ncol = dim(outcorrsparcc$Cor)[1])
sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
sparCCpval <- sparCCpval + t(sparCCpval)


rownames(sparCCpcors) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
colnames(sparCCpcors) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
rownames(sparCCpval) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus
colnames(sparCCpval) <- data.frame(tax_table(ps.betadiversity_woSamplesinOnly1Visit_Visit1_filtered))$Genus

sparCCpcors[1:5, 1:5]
sparCCpval[1:5, 1:5]



sparccCor_processed <- sparCCpcors  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(cor = value)
sparccP_processed <- sparCCpval  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(p = value)

# join the two data frames

SparccP <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
  # # remove self correlations
  # filter(Var1 != Var2) %>% 
  # calculate the false discovery rate to adjust for multiple p values
  mutate(fdr = p.adjust(p, method = "BH"))

fdrThresh <- 0.01 # fdr threshold

sparccOkP <- SparccP%>% filter(fdr < fdrThresh) 

SparccP_plot <- SparccP %>% ggplot(aes(x = Var2, y = Var1, fill = cor)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = "*", size =3)

SparccP_plot
