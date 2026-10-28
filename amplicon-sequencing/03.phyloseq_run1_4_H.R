

############## further analysis of microbiome data

# http://joey711.github.io/phyloseq/index.html
# https://f1000research.com/articles/5-1492/v2

#install.packages("shiny")
#shiny::runGitHub("shiny-phyloseq","joey711")

#BiocManager::install("phyloseq")

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(plyr)
library(memisc)

folder = "G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/phyloseq"

seqtab <- readRDS(paste("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/Dada2/seqtab.rds", sep = "/"))
seqtab[1:2,1:2]
rownames(seqtab) <- gsub("_.*", "", rownames(seqtab))

#taxa <- readRDS(paste("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/Dada2/taxa_dada2_species.rds", sep = "/"))
taxid_species <- readRDS("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/Dada2/taxa_DECIPHER_species.rds")
taxid_species[1:3, 1:3]

samdf <- read.delim("G:/IDIBELL/Data/Bcn_datasets/Epi/Bcn_datasets_Epi_16STissue.txt", header=T)

head(samdf)
table(samdf$Tissue_16S_location, samdf$RNAseq_location)
#          right transverse left
# right       52         40   15
# transvers   36         23   50
# left        28         38   30
# rectum       4          2   13
summary(as.factor(abs(samdf$ID_icobiobank)))
samdf$ID_icobiobank <- abs(samdf$ID_icobiobank)
sum(rownames(seqtab)%in%samdf$ID_icobiobank)
#327

fw_run = read.table("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/highcov_fw.txt", header=T)
fw_run <- fw_run[,c("ID","RunN")]
colnames(fw_run) <- c("ID_icobiobank", "Run_16ST")
head(fw_run)
summary(fw_run)
sum(fw_run$ID_icobiobank%in%samdf$ID_icobiobank)
#417
samdf <- merge(samdf, fw_run, by="ID_icobiobank")
#samdf$Sex <- as.factor(samdf$Sex)
rownames(samdf) <- paste(samdf$ID_icobiobank)
sum(rownames(seqtab)%in%rownames(samdf))
#327
head(samdf)
summary(samdf)
summary(samdf$Origen)
table(samdf$Study, samdf$Origen)
table(samdf$Tissue_16S_location, samdf$RNAseq_location)

fitGTR <- readRDS(paste("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/Dada2/fitGTR.rds", sep = "/"))


#We now construct a phyloseq object directly from the dada2 outputs
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxid_species),
               phy_tree(fitGTR$tree))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 23792 taxa and 327 samples ]
#sample_data() Sample Data:       [ 327 samples by 286 sample variables ]
#tax_table()   Taxonomy Table:    [ 23792 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 23792 tips and 23790 internal nodes ]
#refseq()      DNAStringSet:      [ 23792 reference sequences ]
summary(as.factor(sample_data(ps)$Group3))
# Healthy 
# 327 
summary(as.factor(sample_data(ps)$Group2))
# Healthy 
# 327 
ntaxa(ps)
# 23792
taxa_names(ps)[1:10]
#"ASV1"  "ASV2"  "ASV3"  "ASV4"  "ASV5"  "ASV6"  "ASV7"  "ASV8"  "ASV9"  "ASV10"
nsamples(ps)
# 327
sample_names(ps)[1:5]
# "1009082290" "1064772642" "1090324508" "1119309403" "1122979476"
rank_names(ps)
# "domain"  "phylum"  "class"   "order"   "family"  "Genus"   "Species"
sample_variables(ps)
# "ID_icobiobank" ....Sex" "Age".....
otu_table(ps)[1:5, 1:5]
#OTU Table:          [5 taxa and 5 samples]
#taxa are columns
#           ASV1 ASV2 ASV3 ASV4 ASV5
#1009082290 5184 1344  441 2180  220
#1064772642 6385  534 7938 3705 4928
#1090324508  934  173    3  244 1600
#1119309403 1577  929    6    1  142
#1122979476  827 7410 1077  494 1545
tax_table(ps)[1:5,]
#Taxonomy Table:     [5 taxa by 7 taxonomic ranks]:
#      domain    phylum       class        order             family            Genus              Species        
#ASV1 "Bacteria" "Firmicutes" "Clostridia" "Oscillospirales" "Ruminococcaceae" "Faecalibacterium" "prausnitzii"  
#ASV2 "Bacteria" "Firmicutes" "Clostridia" "Oscillospirales" "Ruminococcaceae" "Faecalibacterium" "prausnitzii"  
#ASV3 "Bacteria" "Firmicutes" "Clostridia" "Oscillospirales" "Ruminococcaceae" "Faecalibacterium" "prausnitzii"
#ASV4 "Bacteria" "Firmicutes" "Clostridia" "Lachnospirales"  "Lachnospiraceae" "Ruminococcus"     NA
#ASV5 "Bacteria" "Firmicutes" "Clostridia" "Lachnospirales"  "Lachnospiraceae" "Blautia"          NA             
summary(tax_table(ps))
#sort(table(tax_table(ps)[, "class"], exclude = NULL), decreasing=T)
sort(table(tax_table(ps)[, "phylum"], exclude = NULL), decreasing=T)
#      Firmicutes       Bacteroidetes                <NA>      Actinobacteria 
#           14608                5083                2867                 461 
#  Proteobacteria        Fusobacteria       Cyanobacteria   Verrucomicrobiota 
#             370                 151                  64                  37 
#    Synergistota     Patescibacteria    Campilobacterota    Desulfobacterota 
#              33                  30                  26                  19 
#     Chloroflexi        Deinococcota       Spirochaetota     Elusimicrobiota
#               8                   8                   7                   4 
# Gemmatimonadota    Thermoplasmatota     Acidobacteriota       Euryarchaeota
#               4                   3                   2                   2
# Planctomycetota    Bdellovibrionota         Myxococcota         Sumerlaeota
#               2                   1                   1                   1
phy_tree(ps)
#Phylogenetic tree with 23792 tips and 23790 internal nodes.
#Tip labels:
#  ASV1, ASV2, ASV3, ASV4, ASV5, ASV6, ...
#Unrooted; includes branch lengths.

rm(dna)
rm(fitGTR)
rm(taxid_species)

######### Preprocessing
#### Taxonomic filtering ##
sort(table(tax_table(ps)[, "domain"], exclude = NULL), decreasing=T)
# Bacteria     <NA>  Archaea 
#    23239      548        5 
#Domain of NA are probably artifacts in a dataset like this, and should be removed.
ps0 <- subset_taxa(ps, !is.na(domain) & !domain %in% c("Archaea"))

sort(table(tax_table(ps0)[, "phylum"], exclude = NULL), decreasing=T)
#Phylum of NA are probably artifacts in a dataset like this, and should be removed.
ps0 <- subset_taxa(ps0, !is.na(phylum) & !phylum %in% c("", "uncharacterized"))
ps0
# otu_table()   OTU Table: [ 20920 taxa and 327 samples ]
# Compute prevalence of each feature, store as data.frame
prevdf0 = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf0 = data.frame(Prevalence = prevdf0,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
head(prevdf0)
#       Prevalence TotalAbundance   domain   phylum     class      order         family          Genus            Species
# ASV1  327(n ind)  822987 (n seq)  Bacteria Firmicutes Clostridia Clostridiales Ruminococcaceae Faecalibacterium prausnitzii
# ASV2  326         658745          Bacteria Firmicutes Clostridia Clostridiales Ruminococcaceae Faecalibacterium prausnitzii

# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf0, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#              phylum Mean_Prev Sum_Prev
#1    Acidobacteriota  1.000000      2 (mean observation when observed; total of observations)
#2   Actinobacteriota 17.247289   7951
#3       Bacteroidota  7.027149  35719
#4   Bdellovibrionota  1.000000      1
#5   Campilobacterota 12.500000    325
#6        Chloroflexi  1.625000     13
#7      Cyanobacteria  6.328125    405
#8       Deinococcota  7.000000     56
#9   Desulfobacterota  7.315789    139
#10   Elusimicrobiota  3.000000     12
#11        Firmicutes 10.270058 150025
#12    Fusobacteriota 10.298013   1555
#13   Gemmatimonadota  1.000000      4
#14       Myxococcota  1.000000      1
#15   Patescibacteria 12.033333    361
#16   Planctomycetota  1.000000      2
#17    Proteobacteria  7.205405   2666
#18     Spirochaetota  6.285714     44
#19       Sumerlaeota  3.000000      3
#20      Synergistota  7.212121    238
#21 Verrucomicrobiota 20.540541    760

#### Rarefaction curve
sample_sums(ps0)
set.seed(42)

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(ps0, c('Observed', 'Shannon', 'Simpson'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary$Sample <- gsub("X","",rarefaction_curve_data_summary$Sample)
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(ps0)), by.x = 'Sample', by.y = 'row.names')

Rare_plot <- ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(x = Depth,y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    group = Sample)) + 
  geom_line() + 
  geom_pointrange() + 
  facet_wrap(facets = ~ Measure, scales = 'free_y')
ggsave(paste(folder,"/00.Rarefaction.tiff",sep=""), Rare_plot, width=10, height=4, units="in", dpi=320)

Rare_plot <- ggplot(
  data = rarefaction_curve_data_summary_verbose[rarefaction_curve_data_summary_verbose$Depth<25000,],
  mapping = aes(x = Depth,y = Alpha_diversity_mean,
                ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                group = Sample)) + 
  geom_line() + 
  geom_pointrange() + 
  facet_wrap(facets = ~ Measure, scales = 'free_y')
ggsave(paste(folder,"/00.Rarefaction25000.tiff",sep=""), Rare_plot, width=10, height=4, units="in", dpi=320)

#### Rarefaction
ps1 <- rarefy_even_depth(ps0, 10000, verbose = FALSE)
# otu_table()   OTU Table:         [ 11619 taxa and 310 samples ]
# sample_data() Sample Data:       [ 310 samples by 286 sample variables ]
# tax_table()   Taxonomy Table:    [ 11619 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 11619 tips and 11617 internal nodes ]
# refseq()      DNAStringSet:      [ 11619 reference sequences ]

summary(as.factor(sample_data(ps1)$Group3))
# Healthy 
# 310 
summary(as.factor(sample_data(ps1)$Group2))
# Healthy 
# 310
rm(ps); rm(ps0); rm(seqtab); rm(rarefaction_curve_data)

#### Agglomerate taxa ##
## per species
ps2 = tax_glom(ps1, "Species", NArm = F)
# otu_table()   OTU Table:         [ 604 taxa and 310 samples ]
# sample_data() Sample Data:       [ 310 samples by 286 sample variables ]
# tax_table()   Taxonomy Table:    [ 604 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 604 tips and 603 internal nodes ]
# refseq()      DNAStringSet:      [ 604 reference sequences ]

# Plotting ASV trees
multiPlotTitleTextSize = 10
p1tree = plot_tree(ps1, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration; 11619 taxa: ASV") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left", 
                   title = "604 taxa: species-level") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
tiff(paste(folder,"/01.Tree_ASV_vs_TaxaSpecies.tiff",sep=""), units="in", res=320, width=6, height=8)
grid.arrange(nrow = 1, p1tree, p2tree)
dev.off()

######## Removing Phenos and Outliers

## Outlier detection
# Bray-Curtis distances
ord.nmds.bray <- ordinate(ps2, method="NMDS", distance="bray")
MDS_coord <- as.data.frame(ord.nmds.bray$points)
head(MDS_coord)
summary(MDS_coord)
plot(MDS_coord)

Outlies_coord1 <- rownames(MDS_coord)[which(abs(MDS_coord$MDS1)>sd(MDS_coord$MDS1)*6)]
Outlies_coord2 <- rownames(MDS_coord)[which(abs(MDS_coord$MDS2)>sd(MDS_coord$MDS2)*6)]
Outliers <- unique(c(Outlies_coord1, Outlies_coord2))

ps2 = subset_samples(ps2, !rownames(sample_data(ps2))%in%Outliers)

ord.nmds.bray <- ordinate(ps2, method="NMDS", distance="bray")
saveRDS(ord.nmds.bray, paste(folder, "Phyloseq_Agglomerate_Species_ord.nmds.bray.rds", sep="/"))

MDS_coord2 <- as.data.frame(ord.nmds.bray$points)
plot(MDS_coord2)

tiff(paste0(folder,"/02.MDS_coord_Outliers.tiff"), units="in", res=320, width=8, height=4)
par(mfrow=c(1,2))
plot(MDS_coord)
plot(MDS_coord2)
dev.off()

# ASV abundance
Abun <- as.data.frame(otu_table(ps2))
plotAbu <-qplot(colSums(Abun>0),binwidth=20) +
  xlab("N samples")+
  ylab("N taxa")
tiff(paste(folder,"/03.Taxa_Presence_ASV.tiff",sep=""), res = 300, width=2000, height=1000)
plotAbu
dev.off()

## Prevalence
library(microbiome)
cluster_taxa <- c("phylum", "class", "order", "family", "Genus", "Species")

for(cl_taxa in cluster_taxa){
  
  if(!cl_taxa=="Species"){
    cluster_data = aggregate_taxa(ps2, cl_taxa)
  }else{
    cluster_data <- ps2
  }
  prevdf0 = apply(X = otu_table(cluster_data),
                  MARGIN = ifelse(taxa_are_rows(cluster_data), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf0 = data.frame(Prevalence = prevdf0,
                       TotalAbundance = taxa_sums(cluster_data),
                       tax_table(cluster_data))
  
  prevdf0$Prevalence100 <- prevdf0$Prevalence/nsamples(ps2)
  prevdf0$Taxa <- row.names(prevdf0)
  write.table(prevdf0, paste0(folder, "/Prevalence_taxa_",cl_taxa,".txt"), sep="\t", quote=F, col.names=T, row.names=F)
  
  Core80.taxa <- core_members(cluster_data, prevalence = 80/100, detection = 1/100)
  writeLines(Core80.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_080.txt"))
  Core50.taxa <- core_members(cluster_data, prevalence = 50/100, detection = 1/100)
  writeLines(Core50.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_050.txt"))
  Core20.taxa <- core_members(cluster_data, prevalence = 20/100, detection = 1/100)
  writeLines(Core20.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_020.txt"))
  Core10.taxa <- core_members(cluster_data, prevalence = 10/100, detection = 1/100)
  writeLines(Core10.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_010.txt"))
  Core05.taxa <- core_members(cluster_data, prevalence = 5/100, detection = 1/100)
  writeLines(Core05.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_005.txt"))
  Core01.taxa <- core_members(cluster_data, prevalence = 1/100, detection = 1/100)
  writeLines(Core01.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_001.txt"))
  All.taxa <- core_members(cluster_data, prevalence = 1/10000000, detection = 1/100)
  writeLines(All.taxa, paste0(folder,"/CoreTaxa_",cl_taxa,"_100.txt"))
}

prevdf0 = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf0 = data.frame(Prevalence = prevdf0,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))
# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf0, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#              phylum        1     2
#1    Acidobacteriota  1.00000     1
#2   Actinobacteriota 27.87805  3429
#3       Bacteroidota 63.91964  7159
#4   Bdellovibrionota  1.00000     1
#5   Campilobacterota 11.86667   178
#6        Chloroflexi  1.75000     7
#7      Cyanobacteria 23.22222   209
#8       Deinococcota 11.50000    23
#9   Desulfobacterota 36.00000    36
#10   Elusimicrobiota  7.00000     7
#11        Firmicutes 81.44344 17999
#12    Fusobacteriota 32.62500   522
#13   Gemmatimonadota  2.00000     2
#14   Patescibacteria 41.66667   125
#15   Planctomycetota  1.00000     2
#16    Proteobacteria 11.45946   848
#17     Spirochaetota  4.25000    17
#18       Sumerlaeota  1.00000     1
#19      Synergistota 18.00000   108
#20 Verrucomicrobiota 47.57143   333

#### Prevalence Filtering ##
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps2)
prevalenceThreshold
## 15.45
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf0)[(prevdf0$Prevalence >= prevalenceThreshold)]
ps05 = prune_taxa(keepTaxa, ps2)
# otu_table()   OTU Table: [ 268 taxa and 309 samples ]

# Compute prevalence of each feature, store as data.frame
prevdf05 = apply(X = otu_table(ps05),
               MARGIN = ifelse(taxa_are_rows(ps05), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf05 = data.frame(Prevalence = prevdf05,
                    TotalAbundance = taxa_sums(ps05),
                    tax_table(ps05))
head(prevdf05)
# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf05, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#              phylum         1     2
#1   Actinobacteriota  86.19444  3103
#2       Bacteroidota 120.94737  6894
#3   Campilobacterota  41.33333   124
#4      Cyanobacteria  98.50000   197
#5       Deinococcota  16.00000    16
#6   Desulfobacterota  36.00000    36
#7         Firmicutes 126.83333 17503
#8     Fusobacteriota  92.20000   461
#9    Patescibacteria 114.00000   114
#10    Proteobacteria  40.82353   694
#11      Synergistota  35.00000    70
#12 Verrucomicrobiota  79.25000   317

#Plotting Abundance vs Prevalence
ASV_Prev <- ggplot(prevdf05, aes(TotalAbundance, Prevalence / nsamples(ps05),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
ggsave(paste(folder,"/04.ASV_Prev_Species05.tiff",sep=""), ASV_Prev, width=7, height=4.3, units="in", dpi=320)

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(ps2)
prevalenceThreshold
## 30.9
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf0)[(prevdf0$Prevalence >= prevalenceThreshold)]
ps10 = prune_taxa(keepTaxa, ps2)
# otu_table()   OTU Table: [ 205 taxa and 309 samples ]

# Compute prevalence of each feature, store as data.frame
prevdf10 = apply(X = otu_table(ps10),
                 MARGIN = ifelse(taxa_are_rows(ps10), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf10 = data.frame(Prevalence = prevdf10,
                      TotalAbundance = taxa_sums(ps10),
                      tax_table(ps10))
head(prevdf10)
# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf10, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#              phylum         1     2
#1   Actinobacteriota 104.32143  2921
#2       Bacteroidota 139.18750  6681
#3   Campilobacterota  41.33333   124
#4      Cyanobacteria  98.50000   197
#5   Desulfobacterota  36.00000    36
#6         Firmicutes 154.77982 16871
#7     Fusobacteriota  92.20000   461
#8    Patescibacteria 114.00000   114
#9     Proteobacteria  56.77778   511
#10      Synergistota  46.00000    46
#11 Verrucomicrobiota 100.00000   300

#Plotting Abundance vs Prevalence
ASV_Prev <- ggplot(prevdf10, aes(TotalAbundance, Prevalence / nsamples(ps10),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
ggsave(paste(folder,"/04.ASV_Prev_Species10.tiff",sep=""), ASV_Prev, width=7, height=4.3, units="in", dpi=320)

# Define prevalence threshold as 20% of total samples
prevalenceThreshold = 0.2 * nsamples(ps2)
prevalenceThreshold
## 61.8
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf0)[(prevdf0$Prevalence >= prevalenceThreshold)]
ps20 = prune_taxa(keepTaxa, ps2)
# otu_table()   OTU Table: [ 151 taxa and 309 samples ]

# Compute prevalence of each feature, store as data.frame
prevdf20 = apply(X = otu_table(ps20),
                 MARGIN = ifelse(taxa_are_rows(ps20), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf20 = data.frame(Prevalence = prevdf20,
                      TotalAbundance = taxa_sums(ps20),
                      tax_table(ps20))
head(prevdf20)
# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf20, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#             phylum        1     2
#1  Actinobacteriota 136.61111  2459
#2      Bacteroidota 174.71429  6115
#3     Cyanobacteria  98.50000   197
#4        Firmicutes 191.02439 15664
#5    Fusobacteriota 100.00000   400
#6   Patescibacteria 114.00000   114
#7    Proteobacteria  84.33333   253
#8 Verrucomicrobiota 132.50000   265

#Plotting Abundance vs Prevalence
ASV_Prev <- ggplot(prevdf20, aes(TotalAbundance, Prevalence / nsamples(ps20),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
ggsave(paste(folder,"/04.ASV_Prev_Species20.tiff",sep=""), ASV_Prev, width=7, height=4.3, units="in", dpi=320)

# Define prevalence threshold as 50% of total samples
prevalenceThreshold = 0.5 * nsamples(ps2)
prevalenceThreshold
## 154.5
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf0)[(prevdf0$Prevalence >= prevalenceThreshold)]
ps50 = prune_taxa(keepTaxa, ps2)
# otu_table()   OTU Table: [ 78 taxa and 309 samples ]

# Compute prevalence of each feature, store as data.frame
prevdf50 = apply(X = otu_table(ps50),
                 MARGIN = ifelse(taxa_are_rows(ps50), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf50 = data.frame(Prevalence = prevdf50,
                      TotalAbundance = taxa_sums(ps50),
                      tax_table(ps50))
head(prevdf50)
# total and average prevalences of the features in each phylum.
plyr::ddply(prevdf50, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#             phylum        1     2
#1  Actinobacteriota 230.5000  1383
#2      Bacteroidota 247.1667  4449
#3        Firmicutes 247.7255 12634
#4 Verrucomicrobiota 162.0000   162

#Plotting Abundance vs Prevalence
ASV_Prev <- ggplot(prevdf50, aes(TotalAbundance, Prevalence / nsamples(ps50),color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
ggsave(paste(folder,"/04.ASV_Prev_Species50.tiff",sep=""), ASV_Prev, width=7, height=4.3, units="in", dpi=320)


# Plotting ASV trees
multiPlotTitleTextSize = 10
p1tree = plot_tree(ps05, method = "treeonly",
                   ladderize = "left",
                   title = "268 species-level taxa at 5%") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p2tree = plot_tree(ps10, method = "treeonly",
                   ladderize = "left", 
                   title = "205 species-level taxa at 10%") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps50, method = "treeonly",
                   ladderize = "left", 
                   title = "78 species-level taxa at 50%") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
tiff(paste(folder,"/04.Tree_TaxaSpecies_Prev.tiff",sep=""), units="in", res=320, width=10, height=6)
grid.arrange(nrow = 1, p1tree, p2tree, p3tree)
dev.off()


# Transform data to proportions or transformed to relative abundance
ps2.ra <- transform_sample_counts(ps2, function(otu) otu/sum(otu))

#Counts vs Relative Abundance
plot_abundance=function(physeq,title="",Facet="phylum",Color="phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  mphyseq= psmelt(physeq)  
  mphyseq<- subset(mphyseq, Abundance > 0)
  ggplot(data= mphyseq,mapping=aes_string(x="Sex",y="Abundance",color= Color,fill= Color)) +
    geom_violin(fill=NA) +
    geom_point(size=1,alpha=0.3,position=position_jitter(width=0.3)) +
    facet_wrap(facets= Facet) + 
    scale_y_log10()+
    theme(legend.position="none")
}

plotBefore= plot_abundance(ps2,"")
tiff(paste(folder,"/05.Abundance_Phylum_Counts.tiff",sep=""), res = 300, width=2000, height=1250)
plotBefore
dev.off()
rm(plotBefore)

plotAfter= plot_abundance(ps2.ra,"")
tiff(paste(folder,"/05.Abundance_Phylum_Prop.tiff",sep=""), res = 300, width=2000, height=1250)
plotAfter
dev.off()
rm(plotAfter)

# Combine each plot into one graphic.
plotBefore= plot_abundance(ps2,"")
plotAfter= plot_abundance(ps2.ra,"")
tiff(paste(folder,"/05.Abundance_Phylum_Counts_vs_Prop.tiff",sep=""), res = 300, width=2000, height=2500)
grid.arrange(nrow = 2, plotBefore, plotAfter)
dev.off()
rm(plotBefore)
rm(plotAfter)


#### Counts per sample distribution
summary(rowSums(otu_table(ps2)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10000   10000   10000   10000   10000   10000
plotBefore <- qplot(rowSums(otu_table(ps2)),binwidth=5000) +
  xlab("Raw counts-per-sample")
plotAfter <- qplot(log10(rowSums(otu_table(ps2))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
tiff(paste(folder,"/03.Counts_per_Sample.tiff",sep=""), res = 300, width=4000, height=2000)
grid.arrange(nrow = 1, plotBefore, plotAfter)
dev.off()


saveRDS(ps1, paste(folder, "Phyloseq_ASV_with_Outliers.rds", sep="/"))
saveRDS(ps2, paste(folder, "Phyloseq_Agglomerate_Species.rds", sep="/"))
saveRDS(ps2.ra, paste(folder, "Phyloseq_Agglomerate_Species_RelAbundance.rds", sep="/"))


#Phylum Relative Abundance Bar plot
tiff(paste(folder,"05.BarPlot_Phylum.tiff",sep="/"), units="in", width=15, height=9, res=320)
plot_bar(ps2.ra, fill = "phylum")
dev.off()


### Visualize alpha-diversity:
tiff(paste(folder,"06.Richness_Age_by_Sex.tiff",sep="/"), units="in", width=5, height=3, res=320)
plot_richness(ps2, x="Age", measures=c("Shannon", "Simpson"), color="Sex")
dev.off()
tiff(paste(folder,"06.Richness_BMI_by_TissueLocation.tiff",sep="/"), units="in", width=5, height=3, res=320)
plot_richness(ps2, x="BMI", measures=c("Shannon", "Simpson"), color="Tissue_16S_location")
dev.off()
tiff(paste(folder,"06.Richness_BMI_by_Run16S.tiff",sep="/"), units="in", width=5, height=3, res=320)
plot_richness(ps2, x="BMI", measures=c("Shannon", "Simpson"), color="Run_16ST")
dev.off()

#all samples
Ps2Rich <- estimate_richness(ps2, measures=c("Shannon", "Simpson"))
Ps2Rich$ID_icobiobank <- rownames(Ps2Rich)
head(Ps2Rich)
write.table(Ps2Rich, paste("G:/IDIBELL/Data/Bcn_datasets/16S_tissue/Run1_4_H/BarcUVa_Otus_richness.txt",sep="/"), quote=F, sep="\t", col.names=T, row.names=F)

Samp <- data.frame(sample_data(ps2))
SampRich <- merge(Samp, Ps2Rich, by="ID_icobiobank")
head(SampRich)
summary(as.factor(SampRich$Sex))


Groups <- c("Sex", "Age_Bi", "BMI_Bi", "Smoking", "Study", "Tissue_16S_location", "Run_16ST", 
            "I2_ardor", "I3_ulcera", "I3_hernia_hiato", "I3_helicobacter", "I3_ninguna",
            "I17_ritme", "I18_num_ventre", "I19_laxant", "Diarrhea", "Constipation",
            "Hemorroides", "Fisura_anal",
            "Diabetes2", "Hipertensio2", "Colesterol", "Infart", "Ictus", "Circulatori",
            "Artritis", "Osteoporosis", "Gota", "Anemia", "EPOC", "Piedras_vesic", "Hepatitis", "Tiroides",
            "Infeccion_urinaria",  "Litiasis_renal",
            "Migranya", "Depressio", "Cancer", 
            "A02BC", "A10BA", "A12A","Calci",  
            "B", "B01AC06", 
            "C",
            "C02", "C03", "C07", "C08", "C09", "C10AA",
            "G", "H", "L",
            "I55_antiinflamatorio", "M01A", "Aines", "M05BA", "Bifosfonats", 
            "N02B", "N05BA", "N06A",
            "R", 
            "t_grams_Bi", "t_energy_Bi", "t_water_Bi",
            "t_protein_ani_Bi", "t_protein_veg_Bi", "t_protein_tot_Bi", "tpc_prot_Bi",
            "t_fats_mfa_Bi","tpc_fats_mfa_Bi", "t_fats_pfa_Bi","tpc_fats_pfa_Bi",
            "t_fats_sfa_Bi","tpc_fats_sfa_Bi", "t_fats_total_Bi", "tpc_fats_total_Bi",
            "t_sugars_dig_Bi",  "t_glucid_tot_Bi", "t_polysach_Bi", "t_fiber_Bi", 
            "t_cholesterol_Bi","tpc_cho_Bi", "t_ethanol_Bi", "t_edible_Bi",
            "t_calcio_Bi",  "t_zinc_Bi","t_phosphorus_Bi", "t_iron_Bi", 
            "t_magnesium_Bi", "t_potassium_Bi", "t_sodium_Bi", "t_niacina_Bi",
            "t_folic_Bi", "t_carotenoids_Bi", "t_retinoids_Bi",
            "t_vit_a_Bi", "t_vit_b1_Bi", "t_vit_b12_Bi", "t_vit_b2_Bi",             
            "t_vit_b6_Bi",  "t_vit_c_Bi", "t_vit_d_Bi", "t_vit_e_Bi",
            "Trichopoulou_MDScore", "Buckland_rMEDScore", "Buckland_arMEDScore")
            
gr <- "Sex"
RichAggr <- c()
for(gr in Groups){
  G <- data.frame(summary(as.factor(SampRich[,gr])))
  G$Group <- paste(gr, row.names(G),sep="_")
  colnames(G) <- c("N", "Group")
  
  Aggr <- aggregate(SampRich[, c("Shannon", "Simpson")], list(SampRich[,gr]), mean)
  Aggr$Group.1 <- paste(gr, Aggr$Group.1,sep="_")
  colnames(Aggr) <- c("Group", "Shannon_mean", "Simpson_mean")
  
  RichAggrG <- merge(G, Aggr, by="Group")
  
  Aggr <- aggregate(SampRich[, c("Shannon", "Simpson")], list(SampRich[,gr]), sd)
  Aggr$Group.1 <- paste(gr, Aggr$Group.1,sep="_")
  colnames(Aggr) <- c("Group", "Shannon_sd", "Simpson_sd")
  
  RichAggrG <- merge(RichAggrG, Aggr, by="Group")
  
  res.aov <- aov(SampRich$Shannon ~ SampRich[,gr])
  dfres <- unlist(summary(res.aov))
  RichAggrG$Shannon_P <- c(rep(0, (nrow(RichAggrG)-1)), dfres["Pr(>F)1"])
  
  res.aov <- aov(SampRich$Simpson ~ SampRich[,gr])
  dfres <- unlist(summary(res.aov))
  RichAggrG$Simpson_P <- c(rep(0, (nrow(RichAggrG)-1)), dfres["Pr(>F)1"])
  
  RichAggrG <- RichAggrG[,c("Group", "N", "Shannon_mean", "Shannon_sd", "Shannon_P", "Simpson_mean", "Simpson_sd", "Simpson_P")]
  
  RichAggr <- rbind(RichAggr, RichAggrG)
}
write.table(RichAggr, paste(folder,"06.Richness_by_Groups.txt",sep="/"), quote=F, sep="\t", col.names=T, row.names=F)




###### Distances
tiff(paste(folder,"07.Ordination_by_BMI.tiff",sep="/"), units="in", width=5, height=3, res=320)
plot_ordination(ps2, ord.nmds.bray, color="BMI", title="Bray NMDS")
dev.off()

##PCoA on unweighted UniFrac distance matrix
gp_pcoa_uf <- ordinate(ps2, method="PCoA", distance="unifrac")
plot_ordination(ps2, gp_pcoa_uf, color="Run_16ST", title="Unifrac")

#
## data(enterotype) 
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
# Remove them from the vector methods that require tree
# "unifrac" "wunifrac"    "dpcoa"
#dist_methods <- dist_methods[-(1:3)]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods%in%c("ANY","kulczynski","morisita","chao","cao", "raup"))]
# Loop through each distance method, save each plot to a list, called plist.
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

## Multiple distances by mutliple phenos
for( i in dist_methods ){
  # Calculate distance matrix
  #iDist <- distance(ps.prop, method=i)
  # Calculate ordination
  iMDS  <- ordinate(ps2, "MDS", distance=i)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ps2, iMDS)# , shape="Group3")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

#Combine results
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=Age))
ggsave(paste0(folder,"/07.Distances_by_Age.tiff"), p, units="in", dpi=320,  width=14, height=10)

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=Sex))
ggsave(paste0(folder,"/07.Distances_by_Sex.tiff"), p, units="in", dpi=320,  width=14, height=10)

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=Tissue_16S_location))
ggsave(paste0(folder,"/07.Distances_by_Tissue_16S_location.tiff"), p, units="in", dpi=320,  width=14, height=10)

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=Run_16ST))
ggsave(paste0(folder,"/07.Distances_by_Run_16ST.tiff"), p, units="in", dpi=320,  width=14, height=10)

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=BMI))
ggsave(paste0(folder,"/07.Distances_by_BMI.tiff"), p, units="in", dpi=320,  width=14, height=10)

p = ggplot(df, aes(Axis.1, Axis.2))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Tissue 16S dataset")
p = p + geom_point(aes(colour=Study))
ggsave(paste0(folder,"/07.Distances0_by_Study.tiff"), p, units="in", dpi=320,  width=14, height=10)

#######Graph-based visualization and testing
library(phyloseq)
library(igraph)
library(ggplot2)
library(dplyr)
library(ggnetwork)
library(phyloseqGraphTest)
library(intergraph)
library(gridExtra)
library(reshape2)
library(nlme)
library(tidyverse)

### Minimum Spanning Tree (MST) or K Neareast Neighbors (KNN) graphs with Jaccard dissimilarity or Bray-Curtis dissimilarity.
## Are the different groups coming from the same distribution, based in measured dissimilarity?
summary(sample_data(ps2))

groups <- c("Age_binned", "Sex", "BMI_binned", "Run_16ST", "Tissue_16S_location", "Study", "Origen", "Smoking",
            "A", "B", "C", "G", "H", "J", "L", "M", "N", "R", "Z")

Pvalue_PermTest <- c()
for(gr in groups){
  ps2b.ra = subset_samples(ps2, !is.na(sample_data(ps2)[,gr]))
  
  #permutation test MST with Jaccard distance
  gt_jc_mst <- graph_perm_test(ps2b.ra, gr,
                        distance = "jaccard", type = "mst")
  #permutation test KNN with Jaccard distance
  gt_jc_knn <- graph_perm_test(ps2b.ra, gr,
                        distance = "jaccard", type = "knn", knn = 2)
  #permutation test KNN with Bray-Curtis distance
  gt_bc_knn <- graph_perm_test(ps2b.ra, gr,
                        distance = "bray", type = "knn", knn = 2)
  
  Pvalue_PermTest_gr <- c(gr, gt_jc_mst$pval, gt_jc_knn$pval,gt_bc_knn$pval)
  Pvalue_PermTest <- rbind(Pvalue_PermTest, Pvalue_PermTest_gr)
  
  # MST Graph and Jaccard dissimilarity Permutation Hist
  plotNet1=plot_test_network(gt_jc_mst) + 
    theme(legend.text = element_text(size = 8),legend.title = element_text(size = 9)) 

  plotPerm1=plot_permutations(gt_jc_mst)
  
  # KNN Graph and Bray-Curtis dissimilarity Permutation Hist
  plotNet2=plot_test_network(gt_bc_knn) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
  plotPerm2=plot_permutations(gt_bc_knn)
  grid.arrange(ncol = 2, plotNet2, plotPerm2)
  
  tiff(paste(folder,"/08.Graph_MST_Jaccard_PermHist_NKK_BrayCurtis_PermHist_",gr,".tiff", sep=""), res=300, width=2000, height=1000)
  grid.arrange(ncol = 2, nrow=2, plotNet1, plotPerm1, plotNet2, plotPerm2)
  dev.off()
}

Pvalue_PermTest <- as.data.frame(Pvalue_PermTest)
colnames(Pvalue_PermTest) <- c("Parameter", "PermP_MST_Jc", "PermP_KNN_Jc", "PermP_KNN_BC")
write.table(Pvalue_PermTest, paste(folder,"/08.Graph_MST_Jaccard_PermHist_NKK_BrayCurtis_PermHist_Pvalues.txt", sep=""), quote=F, row.names=F, col.names=T, sep="\t")


install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")

