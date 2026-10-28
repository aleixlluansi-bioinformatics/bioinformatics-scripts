##20210629 SOP Illumina data analyses update

#First of all (Very important!). Every time we create a file, image, Rhistory, etc. they will be locked. To unlock we need to open a new terminal and do:  sudo chmod -R 0777 /home/digestiu...(whatever folder they are)

##depending on in which computer do you want to follow this Standard Operation Procedure (SOP) for illumina sequencing analyses, some initial steps (i.e., 1) will not be necessary
##0_ activate external hard_drive containing all files and folders of research

$ sudo vgchange -ay --activationmode partial; sudo mount /mnt/Dades/

##1_ getting ready. Uncopress data.

$ cd untrimmed_data #cd to the folder of interest where data is stored, will be manipulated and processed. Remember that data must be untrimmed
$ gunzip *.gz
$ ls (#must show all files as .fastq and all other files contained in the folder)

##2_ load R and work. In Linux run R as sudo, whereas in macOS run R normally

(base) digestiu@digespostdoc02:~$ cd /path_to_folder/
(base) digestiu@digespostdoc02:/folder$ sudo R
[sudo] password for digestiu: Idibgi2019


##3_ load packages and redirect work to working folder. All below mentioned packages are already installed, if not, install them using install.package command or biocmanager.

library("dada2");packageVersion("dada2");library("phyloseq");packageVersion("phyloseq");library("tidyverse");packageVersion("tidyverse");library("readxl");library("Biostrings");packageVersion("Biostrings");library("phangorn");packageVersion("phangorn");library("ape");packageVersion("ape");library("microbiome");packageVersion("microbiome");library("vegan");packageVersion("vegan");library("DESeq2");library("genefilter");library("btools");library("xtable");library("scales");library("devtools");library("selbal")
setwd("working_directory") # CHANGE ME to the directory containing the fastq files after unzipping.
path<-"/home/microbioma/Documentos/Microomics-Yolanda/BAL/Reads/Run2/ALL/" # force redirection of files and results to working directory
list.files(path) # list of files in working directory. A way to verify working in desired directory
theme_set(theme_classic())

##4_ Filter and Trim sequences. This must be redone according to your OWN data files
#_Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs<-sort(list.files(path,pattern="_R1_001.fastq",full.names=TRUE))
fnRs<-sort(list.files(path,pattern="_R2_001.fastq",full.names=TRUE))

##5_ Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse = "_"))
sample.names
 

##6_ Inspect read quality profiles of Forward and Reverse reads
#examine quality profiles. [1:3] allows for 1 line 3 plots of the first 3 samples. Images must be saved as screen captures.

plotQualityProfile (fnFs[1:6])
plotQualityProfile (fnRs[1:6])


##7_ Filter and Trim reads and assign the filenames for the filtered fastq.gz files

filtFs<-file.path(path,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs<-file.path(path,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs)<-sample.names
names(filtRs)<-sample.names
out<-filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen=c(260,220),maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE, trimLeft = c(17,21), trimRight = c(17,21), minLen = c(150,150)) #truncLen can vary (e.g, bread-v3v4 280,160; sv-v4 240,120; bread-18Sv9 200,200). Normally, we choose the trunclen parametres according to the green line of the quality plot (for example when it drops down to less than 30).  You can also run FIGARO to knwon trimming positions if quality plot is not clear or we need to optimize the trimming.In 16S not need cutadapt to verify truncLen, but 18S yes.
#FIGARO_ To run figaro tool open a new terminal (out of R) and cd to the folder where figaro tool is downloaded (cd /home/digestiu/figaro/). Then run: figaro -i /home/digestiu..... (folder where the data is) -o /home/digestiu... (folder where you want the output data) -a 428 (amplicon) -f 17(forward primer) -r 20(reverse primer)


##7.1_Report on reads in and out of each analysed file

head(out)

report.PROJECT.table<-data.frame(out)
write.table(report.PROJECT.table,file="report.PROJECT.table.txt",sep="t",row.names=TRUE,col.names=TRUE)

##8_ learn error rates and visualize them. Learning parametric error model. Computational intensive, so only use by default a subset of the data set (the first 1M reads). Images must be saved as screen captures.

errF<-learnErrors(filtFs,multithread=TRUE)
errR<-learnErrors(filtRs,multithread=TRUE)
plotErrors(errF,nominalQ=TRUE) #4x4 plots with curves for each error possibility (A->T, ...)
plotErrors(errR,nominalQ=TRUE) #4x4 plots with curves for each error possibility (A->T, ...)


##9_ Sample inference. Infer the sequence variants in each sample. Data inspection.

dadaFs<-dada(filtFs,err=errF,multithread=TRUE)
dadaRs<-dada(filtRs,err=errR,multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]


##10_ Merge paired reads. We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs script. In your data, most of your reads should successfully merge. If that is not the case, upstream parameters may need to be revisited: Did you trim away the overlap between your reads?

mergers<-mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose=TRUE)

#_Inspect the merger data.frame for the first sample. Return a list of sequences and its abundance, forward, reverse, nmatch, nmismatch, nindel, prefer, accept.

head(mergers[[1]])

##11_ Construct a sequence table. A higher resolution table of the OTU table, so a matrix (rows corresponding to and named as the samples, and columns corresponding to and named as the sequence variants)

seqtab<-makeSequenceTable(mergers)
dim(seqtab)

##11.1_Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

##11.2_Histogram of the number of sequences according to length of the sequence (it's a table). Image should be saved as screen capture

hist(nchar(getSequences(seqtab))),main="Distribution of all PROJECT sequence length")


##12_Chimera detection and removal.

seqtab.nochim<-removeBimeraDenovo(seqtab,method="consensus",multithread=TRUE,verbose=TRUE)

dim(seqtab.nochim)


##12.1_Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab.nochim)))


##12.2_Histogram of the number of sequences according to length of the sequence (it's a table). Image should be saved as screen capture

hist(nchar(getSequences(seqtab.nochim))),main="Distribution of PROJECT non.chimera sequence length")


##13_Fraction of chimera detected

sum(seqtab.nochim)/sum(seqtab)


##14_Track reads through the pipeline. Final check of the process returning a table (head track) of input/filtered/denoised/merged/tabled/nonchim sequences in each sample.

getN<-function(x) sum(getUniques(x))
track<-cbind(out,sapply(dadaFs,getN),sapply(dadaRs,getN),sapply(mergers,getN),rowSums(seqtab.nochim))
colnames(track)<-c("input","filtered","denoisedF","denoisedR","merged","non.chimera")
rownames(track)<-sample.names
head(track)


track.PROJECT.table<-data.frame(track)
write.table(track.PROJECT.table,file="track.PROJECT.table.txt",sep="t",row.names=TRUE,col.names=TRUE)

##15_Taxonomy. It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. The DADA2 package provides a native implementation of the RDP-s naive Bayesian classifier for this purpose. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.

taxa<-assignTaxonomy(seqtab.nochim,"rdp_train_set_18.fa.gz",multithread=TRUE,verbose = T,tryRC = T)
taxa.plus<-addSpecies(taxa,"rdp_species_assignment_18.fa.gz",verbose=TRUE, allowMultiple = T, tryRC = T)



##15.1_save taxonomy table as a file

write.table(taxa,file="taxa.txt",sep="\t",row.names=TRUE,col.names=TRUE) #taxonomy table
write.table(seqtab.nochim,file="seqtab.nochim.txt",sep="\t",row.names=TRUE,col.names=TRUE) #final non.chimeric sequences in a table


##16_Evaluating dada2 accuracy in mock community (positve control, mixture of 20 species from Zymo)

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n"))

mock.ref<-getSequences(file.path(path,"zymoStDNA.fas"))
match.ref<-sum(sapply(names(unqs.mock),function(x) any(grepl(x,mock.ref))))
cat("of those,",sum(match.ref),"were exact mathces to the expected reference sequences.\n")


##17_Phyloseq. R package to interpret data and obtain figures produced with ggplot.

##17.1_Import OTU table

OTU=otu_table(seqtab.nochim,taxa_are_rows=FALSE)


##17.2_Import TAXONOMY table

TAX=tax_table(taxa.plus)

##17.3_Import metadata
#the good option to import metadata

mapping_file_run1 <- mapping_file[mapping_file$SampleID %in% sample.names, ]
rownames(mapping_file_run1)<-mapping_file_run1$SampleID




##18_Construct Phyloseq Object

ps_RUN2<-phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),tax_table(taxa.plus),sample_data(mapping_file_run2))
ps_RUN2


##18.1_visualize dataset
sample_names(ps_RUN2)
  

sample_variables(ps_RUN2) #will vary according to your dataset

rank_names(ps_RUN2)


###DECONTAM
df_cont <- as.data.frame(sample_data(ps_RUN2)) # Put sample_data into a ggplot-friendly data frame
df_cont$LibrarySize <- sample_sums(ps_RUN2)
df_cont <- df_cont[order(df_cont$LibrarySize),]
df_cont$Index <- seq(nrow(df_cont))
ggplot(data=df_cont[1:209,], aes(x=Index, y=LibrarySize, color=Tipo_muestra)) + geom_point()

sample_data(ps_RUN2)$is.neg <- sample_data(ps_RUN2)$Tipo_muestra == "CN"
contamdf.prev05 <- isContaminant(ps_RUN2, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev05$contaminant<-as.factor(contamdf.prev05$contaminant)
table(contamdf.prev05$contaminant)

# Make a phyloseq object of presence-absence in negative controls and true samples
ps_RUN2_pa <- transform_sample_counts(ps_RUN2, function(abund) 1*(abund>0))
ps_RUN2_pa.neg <- prune_samples(sample_data(ps_RUN2_pa)$Tipo_muestra == "CN", ps_RUN2_pa)
ps_RUN2_pa.pos <- prune_samples(sample_data(ps_RUN2_pa)$Tipo_muestra == "BAL", ps_RUN2_pa)
# Make data frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps_RUN2_pa.pos), pa.neg=taxa_sums(ps_RUN2_pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

OTUS_keep <- subset(contamdf.prev05, contaminant=="FALSE", ) # lista de OTU's que nos quedamos
taxa_names<-rownames(OTUS_keep)
ps_RUN2_decontam=prune_taxa(taxa_names, ps_RUN2)


physeq<-subset_samples(ps_RUN2_decontam, Tipo_muestra=="BAL")
physeq<-prune_taxa(taxa_sums(physeq) > 0, physeq)
physeq


##19.2_Remove undesired taxa: chloroplast or mitochondria from phyloseq object ...Examples:


physeq_clean<-subset_taxa(physeq,Phylum!="NA")
physeq_clean<-subset_taxa(physeq_clean,Phylum!="Cyanobacteria")
physeq_clean


dna<-Biostrings::DNAStringSet(taxa_names(physeq_clean))
names(dna)<-taxa_names(physeq_clean)
physeq_clean<-merge_phyloseq(physeq_clean,dna)
taxa_names(physeq_clean)<-paste0("ASV",seq(ntaxa(physeq_clean)))
physeq_clean







library(tidyverse)

library(phyloseq)

library(microbiome)

library(scales)


phy2 <- physeq_clean %>% aggregate_taxa(level = "Family") %>%
  
  microbiome::transform(transform = "compositional")

phy2

phy2 <- physeq_clean %>% aggregate_top_taxa(level = "Family", top = 20) %>%
  
  microbiome::transform(transform = "compositional")

phy2



phy_genus <- tax_glom(physeq_clean, taxrank = "Genus")
phy_genus_rel <- microbiome::transform(phy_genus, transform = "compositional")
top20 <- names(sort(taxa_sums(phy_genus_rel), decreasing = TRUE))[1:20]
# First, make sure taxa names are characters
taxa_names(phy_genus_rel) <- make.unique(taxa_names(phy_genus_rel))

# Mark all taxa not in top20 as "Other"
tax_table(phy_genus_rel)[!taxa_names(phy_genus_rel) %in% top20, "Genus"] <- "Other"

# Agglomerate again to merge all "Other" taxa into one
phy_top20_plus_other <- tax_glom(phy_genus_rel, taxrank = "Genus")


phy2 <- phy_top20_plus_other %>% aggregate_taxa(level = "Genus") %>%
  
  microbiome::transform(transform = "compositional")

phy2


phy_genus_rel %>%  plot_composition(group_by = "Tipo_muestra2")+ scale_y_continuous(labels = percent)


physeq_clean_rel <- transform_sample_counts(physeq_clean, function(x) x / sum(x))
# Extraer datos
ps_data <- psmelt(physeq_clean_rel)  
# psmelt da un data.frame con columnas: Sample, OTU, Abundance, Taxonomy + variables de sample_data

head(ps_data)

# Agrupar por Tipo_muestra y Genus y sumar abundancia total
# Agrupar por Tipo_muestra2 y Genus y sumar abundancia total
top20_by_tipo <- ps_data %>%
  group_by(Tipo_muestra2, Genus) %>%
  summarise(total_abund = sum(Abundance), .groups = "drop") %>%
  arrange(Tipo_muestra2, desc(total_abund)) %>%
  group_by(Tipo_muestra2) %>%
  mutate(rank = row_number()) %>%
  mutate(Genus_grouped = ifelse(rank > 20, "Other", Genus)) %>%
  select(-rank)

top20_by_tipo <- ps_data %>%
  left_join(top20_by_tipo %>% select(Tipo_muestra2, Genus, Genus_grouped),
            by = c("Tipo_muestra2", "Genus")) %>%
  mutate(Genus_final = Genus_grouped) %>%
  group_by(Sample, Tipo_muestra2, Genus_final) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")





library(scales)

library(RColorBrewer)
library(scales)

# Niveles con "Other" al final
genus_levels <- unique(top20_by_tipo$Genus_final)
genus_levels <- c(setdiff(genus_levels, "Other"), "Other")
top20_by_tipo$Genus_final <- factor(top20_by_tipo$Genus_final, levels = genus_levels)

# Número de géneros menos "Other"
n_taxa <- length(genus_levels) - 1

# Paleta Paired, que es colorida y contrastada
base_colors <- brewer.pal(min(12, n_taxa), "Paired")

# Si necesitas más colores, expandimos interpolando
if (n_taxa > length(base_colors)) {
  base_colors <- colorRampPalette(base_colors)(n_taxa)
}

# Asignar gris claro sólo para "Other"
genus_colors <- c(base_colors, "Other" = "#d3d3d3")
names(genus_colors) <- genus_levels

# Plot
ggplot(top20_by_tipo, aes(x = Sample, y = Abundance, fill = Genus_final)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  facet_wrap(~ Tipo_muestra2, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = genus_colors, name = "Género") +
  labs(y = "Abundancia relativa", x = "Muestra") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "right")


##20_Control number of reads and plot them

readsumsdf=data.frame(nreads=sort(taxa_sums(physeq_clean),TRUE),sorted=1:ntaxa(physeq_clean),type="ASVs")
readsumsdf=rbind(readsumsdf,data.frame(nreads=sort(sample_sums(physeq_clean),TRUE),sorted=1:nsamples(physeq_clean),type="samples"))
title="total number of good quality reads in XX dataset"
physeq_clean.ReadDistribution<-ggplot(readsumsdf,aes(x=sorted,y=nreads))+geom_bar(stat="identity")+ggtitle(title)+scale_y_log10(labels=trans_format("log10",math_format(10^.x)))+facet_wrap(~type,1,scales="free")
physeq_clean.ReadDistribution
ggsave("physeq_clean.png",dpi=300)

##21_Construct tree of samples and adding phylogenetic tree to a phyloseq object
##best_option Construct ARB related tree

EXP.noposneg.SN1.clean3.def.tax.table<-tax_table(EXP.noposneg.SN1.clean3.def)
write.table(EXP.noposneg.SN1.clean3.def.tax.table,file="EXP.noposneg.SN1.clean3.def.tax.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)
refseq(EXP.noposneg.SN1.clean3.def)

EXP.noposneg.SN1.clean3.def.refseq<-refseq(EXP.noposneg.SN1.clean3.def)
write.table(EXP.noposneg.SN1.clean3.def.refseq,file="EXP.noposneg.SN1.clean3.def.refseq.txt",sep="\t",row.names=TRUE,col.names=TRUE)

#now open EXP.noposneg.SN1.clean3.def.refseq.txt file and replace "x" and "ASVx" by its corresponding format for "fasta format"
#save this new file as EXP.noposneg.SN1.clean3.def.refseq.fasta
#once saved, open new terminal window and work with SINA aligner and do ...
#$ cd /mnt/Dades/sv16SuntrimDada2BisDraftDef/
#$ conda activate sina 
#$ sina -i EXP.noposneg.SN1.clean3.def.refseq.fasta -r SILVA_138.1_SSURef_NR99_12_06_20_opt.arb -o EXP.noposneg.SN1.clean3.def.refseq.aligned.fasta

#then, use aligned fasta file to construct a tree as you will do it using ARB software, but using console script and fasttree tool
#$ fasttree -gamma -gtr -nt EXP.noposneg.SN1.clean3.def.refseq.aligned.fasta EXP.noposneg.SN1.clean3.def.refseq.aligned.tree


#now we come back to terminal screen of phyloseq and continue working

EXP.noposneg.SN1.clean3.def.refseq.tree<-read.tree("EXP.noposneg.SN1.clean3.def.refseq.aligned.tree")
EXP.noposneg.SN1.clean3.def.refseq.tree


##22_Merge tree to phyloseq object

EXP.noposneg.SN1.clean3.def.tree=merge_phyloseq(EXP.noposneg.SN1.clean3.def.refseq.tree)
EXP.noposneg.SN1.clean3.def

EXP.noposneg.SN1.clean3.def.tree

sample_data(EXP.noposneg.SN1.clean3.def.tree)

EXP.noposneg.SN1.clean3.def.tree.backup=EXP.noposneg.SN1.clean3.def.tree


##23_Summarize infos of the phyloseq object appart of the phyloseq object information of screen

summarize_phyloseq(EXP.noposneg.SN1.clean3.def.tree) #retrieve information on number of reads, sparsity, %singletons, sample variables ...

##24_Rarefaction curve. Figure must be saved as screen capture

EXP.noposneg.SN1.clean3.def.tree.raref.curve<-rarecurve(t(otu_table(EXP.noposneg.SN1.clean3.def.tree)),step=50,cex=0.5)


##25______Alpha Diversity______
##25.1_Write in screen first values of richness estimates of phyloseq object

head(estimate_richness(EXP.noposneg.SN1.clean3.def.tree))

##25.2_Write richness estimates of phyloseq object in a file

EXP.noposneg.SN1.clean3.def.tree.div.table<-estimate_richness(EXP.noposneg.SN1.clean3.def.tree,split=TRUE,measures=c("Observed","Chao1","Shannon","ACE","Simpson","InvSimpson","Fisher"))
EXP.noposneg.SN1.clean3.def.tree.div.table
write.table(EXP.noposneg.SN1.clean3.def.tree.DivTable,file="EXP.noposneg.SN1.clean3.def.tree.div.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)

##25.3_Generate diversity indices and visualize them as alpha-diversity boxplots

EXP.noposneg.SN1.clean3.def.tree.diversity<-plot_richness(EXP.noposneg.SN1.clean3.def.tree,x="group",measures=c("Observed","Chao1","Shannon"),color="group")+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=10)) 
EXP.noposneg.SN1.clean3.def.tree.diversity
ggsave("EXP.diversity.png",dpi=300) # dpi can also be at 700, and also can be saved as pdf or eps
EXP.noposneg.SN1.clean3.def.tree.diversity2<-plot_richness(EXP.noposneg.SN1.clean3.def.tree,x="group",measures=c("Observed","Chao1","Shannon"),color="group")+geom_boxplot(aes(fill="group"),alpha=0.2)+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=10))
EXP.noposneg.SN1.clean3.def.tree.diversity2
ggsave("EXP.diversity2.png",dpi=300) # dpi can also be at 700, and also can be saved as pdf or eps

####25.4_ANOVA analyses using Faith's Phylogenetic Diversity index. According to ISME LAB 2019. PENDING TO VERIFY AND DO
#safe in a dataframe the alpha diversity values computed
alpha_pd <- estimate_pd(EXP.noposnegTree)
#combine alpha diversity and metadata
data <- cbind(sample_data(EXP.noposnegTree), alpha_pd) 
#and compute ANOVA
EXP.noposnegTree.anova <- aov(PD ~ Group, data) 
#export as excel table
install.packages("xtable");library(xtable);packageVersion(xtable)
EXP.noposnegTree.anova.table <- xtable(EXP.noposnegTree.anova)
#to obtain xtable we need to do: print(xtable(EXP.noposnegTree.anova.table), include.rownames=T, include.colnames=T) but probably the exported table will be in latex format. Use any tool from the net to convert.

###26_ Save OTU and TAXONOMY Tables of the definitive phyloseq object

head(otu_table(EXP.noposneg.SN1.clean3.def.tree)) #or otu_table(EXP.noposneg.SN1.clean3.def.tree)[1:5,1:6] that will produce a 6 taxa (columns) and 5 samples (rows) table in screen

EXP.noposneg.SN1.clean3.def.tree.otu.table<-data.frame(otu_table(EXP.noposneg.SN1.clean3.def.tree))
write.table(EXP.noposneg.SN1.clean3.def.tree.otu.table,file="EXP.noposneg.SN1.clean3.def.tree.otu.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)
tax_table(EXP.noposneg.SN1.clean3.def.tree)[1:5,1:6] #or head(otu_table(EXP.noposneg.SN1.clean3.def.tree)) that will produce a XX taxa (columns) and 6 samples (rows) table in screen

EXP.noposneg.SN1.clean3.def.tree.tax.table<-data.frame(tax_table(EXP.noposneg.SN1.clean3.def.tree))
write.table(EXP.noposneg.SN1.clean3.def.tree.tax.table,file="EXP.noposneg.SN1.clean3.def.tree.tax.table.txt",sep="\t",row.names=TRUE,col.names=TRUE) #taxonomy table of all OTUs
EXP.noposneg.SN1.clean3.def.tree.otu.table<-data.frame(otu_table(EXP.noposneg.SN1.clean3.def.tree))
write.table(EXP.noposneg.SN1.clean3.def.tree.otu.table,file="EXP.noposneg.SN1.clean3.def.tree.otu.table.txt",sep="\t",row.names=TRUE,col.names=TRUE) #otu table of all OTUs

EXP.noposneg.SN1.clean3.def.tree.bis=EXP.noposneg.SN1.clean3.def.tree

##26_To generate Rank Abundance plots at "Phylum" (or other) level with the top 100 OTUs [1:100] or top 20 OTUs [1:20]. We can reproduce the same results by not doing relative abundance computations using the commented clusterData2 scripts 

EXP.noposneg.SN1.clean3.def.tree.RelAb<-transform_sample_counts(EXP.noposneg.SN1.clean3.def.tree,function(x) 100*x/sum(x))
EXP.noposneg.SN1.clean3.def.tree.RelAb

tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb)
otu_table(EXP.noposneg.SN1.clean3.def.tree.RelAb)[1:5,1:6]

clusterData=psmelt(EXP.noposneg.SN1.clean3.def.tree.RelAb)
clusterData=filter(clusterData,Abundance>0)

clusterAggPhylum=aggregate(Abundance ~ OTU+Phylum,data=clusterData,mean)
clusterAggPhylum=clusterAggPhylum[order(-clusterAggPhylum$Abundance),] #[1:100,] can be replaced by [1:20,] (generated with "detail" on file name) or removed to plot all OTUs. You can play with [1:XXX,]
clusterAggPhylumDetail=clusterAggPhylum[order(-clusterAggPhylum$Abundance),][1:20,] #[1:20,] can be replaced by [1:100,] to generate the clustering using top100 otus or removed to plot all OTUs. You can play with [1:XXX,]
EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.all<-ggplot(clusterAggPhylum,aes(x=reorder(OTU,-Abundance),y=Abundance))+geom_point(aes(color=Phylum),size=2)+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=6))+scale_y_log10(labels=trans_format("log10",math_format(10^.x))) #Replace "Phylum" by other taxonomy order
EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.all
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.all.png",dpi=300) #save as pdf and eps formats

EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.detail<-ggplot(clusterAggPhylumDetail,aes(x=reorder(OTU,-Abundance),y=Abundance))+geom_point(aes(color=Phylum),size=2)+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=8))+scale_y_log10(labels=trans_format("log10",math_format(10^.x))) #Replace "Phylum" by other taxonomy order
EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.detail
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.RankAbun.Phylum.detail.png",dpi=300) #save as pdf and eps formats


#clusterData2=psmelt(EXP.noposneg.SN1.clean3.def.tree)
#clusterData2=filter(clusterData2,Abundance>0)
#clusterAggPhylum2=aggregate(Abundance ~ OTU+Phylum,data=clusterData2,mean)
#clusterAggPhylum2Detail=clusterAggPhylum2[order(-clusterAggPhylum2$Abundance),][1:20,]

#In case to reproduce the same plots as above using other taxonomic rank (Family for example), MY ADVICE is to redo the clustering and aggregation as below:
#clusterData3=psmelt(EXP.noposneg.SN1.clean3.def.tree.RelAb)
#clusterData3=filter(clusterData3,Abundance>0)
#clusterAggFamily=aggregate(Abundance~OTU+Family,data=clusterData3,mean)
#clusterAggFamily=clusterAggFamily[order(-clusterAggFamily$Abundance),]
#clusterAggFamilyDetail=clusterAggFamily[order(-clusterAggFamily$Abundance),][1:20,]


##27_Barplots
##27.1_Barplot at phylum level of reads (with abundance considered)

EXP.noposneg.SN1.clean3.def.tree.RelAb.barplot.phylum<-plot_bar(EXP.noposneg.SN1.clean3.def.tree.RelAb,"group",fill="Phylum")+geom_bar(aes(),stat="identity",position="stack")+theme(axis.title.x=element_blank(),axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=10),legend.text=element_text(size=8),legend.title=element_text(size=8))
EXP.noposneg.SN1.clean3.def.tree.RelAb.barplot.phylum
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.barplot.phylum.png",dpi=300)

##27.2_Barplot at phylum level of reads (without relative abundance considered)

EXP.noposneg.SN1.clean3.def.tree.barplot.phylum<-plot_bar(EXP.noposneg.SN1.clean3.def.tree,"group",fill="Phylum")+geom_bar(aes(),stat="identity",position="stack")+theme(axis.title.x=element_blank(),axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=10),legend.text=element_text(size=8),legend.title=element_text(size=8))
EEXP.noposneg.SN1.clean3.def.tree.barplot.phylum
ggsave("EXP.noposneg.SN1.clean3.def.tree.barplot.phylum.png",dpi=300)

##27.3_Generate the top 20 first ASVs and generate a relative abundance barplot of them in each group merged at Family level. For other taxonomic levels replace "Family" by "Genus" in the above scripts

EXP.noposneg.SN1.clean3.def.tree.RelAb.bis=EXP.noposneg.SN1.clean3.def.tree.RelAb #to avoid extra columns in the original dataset, I propose to copy the phyloseq object to generate TOP20 barplot
top20otus=names(sort(taxa_sums(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis),TRUE)[1:20])
taxtab20=cbind(tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis),family20=NA)
taxtab20[top20otus,"family20"]<-as(tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis)[top20otus,"Family"],"character")
tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis)<-tax_table(taxtab20)
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged=merge_samples(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis,"group")
sample_data(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged)$group<-levels(sample_data(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis)$group)
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged=transform_sample_counts(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged,function(x) 100* x/sum(x))
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.barplot<-plot_bar(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged,fill="family20",title="top 20 otus @family")+ylab("percentage of retrieved sequences")+ylim(0,100)+xlab("group of subjects")+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=8))
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.barplot
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.barplot.png",dpi=300)

EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20=prune_taxa(top20otus,EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged) #to remove NAs from the graph
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.merged.barplot<-plot_bar(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20,fill="family20",title="top 20 otus @family")+ylab("percentage of retrieved sequences")+ylim(0,100)+xlab("group of subjects")+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=10))
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.merged.barplot
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.merged.top20.family.merged.barplot.png",dpi=300)


##27.4_Generate the top 20 first ASVs and generate a relative abundance barplot of them in each group merged at Genus level.
###THE GOOD ONE###

ps_betadiversity_filtered.genus.bis<-ps_betadiversity_filtered.genus
top20otus=names(sort(taxa_sums(ps_betadiversity_filtered.genus.Relab),TRUE)[1:20])
taxtab20=cbind(tax_table(ps_betadiversity_filtered.genus.bis),genus20=NA)
taxtab20[top20otus,"genus20"]<-as(tax_table(ps_betadiversity_filtered.genus.bis)[top20otus,"Genus"],"character")
tax_table(ps_betadiversity_filtered.genus.bis)<-tax_table(taxtab20)

ps_betadiversity_filtered.genus.bis.merged=merge_samples(ps_betadiversity_filtered.genus.bis,"Time")


ps_betadiversity_filtered.genus.bis.merged.Relab=transform_sample_counts(ps_betadiversity_filtered.genus.bis.merged,function(x) 100* x/sum(x))


ps_betadiversity_filtered.genus.bis.merged.Relab.top20.genus.barplot<-plot_bar(ps_betadiversity_filtered.genus.bis.merged.Relab,fill="genus20",title="top 20 otus @genus")+ylab("percentage of retrieved sequences")+ylim(0,100)+xlab("group of subjects")+theme(axis.text.x=element_text(hjust=1,vjust=1,angle=45,size=8))
ps_betadiversity_filtered.genus.bis.merged.Relab.top20.genus.barplot


ggsave("ps_betadiversity_filtered.genus.bis.merged.Relab.top20.genus.barplot.png",dpi=300)


##28______Beta Diversity______

##28.1_Generate PCoA plots without or with facet wrapping.

EXP.noposneg.SN1.clean3.def.tree.RelAb #read phyloseq object again to verify

ordPCoA=ordinate(EXP.noposneg.SN1.clean3.def.tree.RelAb,method="PCoA",distance="unifrac",weighted=TRUE)
EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA<-plot_ordination(EXP.noposneg.SN1.clean3.def.tree.RelAb,ordPCoA,color="group",title="XX project PCoA unifrac weighted")+scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))
EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.png",dpi=300)

EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.Phylum<-plot_ordination(EXP.noposneg.SN1.clean3.def.tree.RelAb,ordPCoA,type="taxa",color="Phylum",title="XX project PCoA unifrac weighted")
EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.Phylum
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.Phylum.png",dpi=300)

EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.facetFamily<-plot_ordination(EXP.noposneg.SN1.clean3.def.tree.RelAb,ordPCoA,type="taxa",color="Family",title="XX project PCoA unifrac weighted")+scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))+facet_wrap(~Family,3) #replace taxonomy level to do it at genus, class ...
EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.facetFamily
ggsave("EXP.noposneg.SN1.clean3.def.tree.RelAb.PCoA.facetFamily.png",dpi=300)


##28.2_Generate NMDS plots without or with facet wrapping. Similar way as for PCoA using Unifrac. Only depict here the ord command and the rest is as above.

ordNMDS=ordinate(EXP.noposneg.SN1.clean3.def.tree.RelAb,method="NMDS",distance="bray")


##28.3_Supported ordination method. Check page https://joey711.github.io/phyloseq/plot_ordination-examples.html and there the Supported ordination methods list and complete script.


##29_Write Taxonomy and OTU tables into a file

#TAX and OTU table for merged OTUs by group
EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged

sample_data(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged)

EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.otu.table<-data.frame(otu_table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged)
write.table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.otu.table,file="EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.otu.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)

EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.tax.table<-data.frame(tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged))
write.table(EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.tax.table,file="EXP.noposneg.SN1.clean3.def.tree.RelAb.bis.bis.merged.tax.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)

#TAX and OTU tables for original phyloseq object based on 100% relative abundances
EXP.noposneg.SN1.clean3.def.tree.RelAb
EXP.noposneg.SN1.clean3.def.tree.RelAb.otu.table<-data.frame(otu_table(EXP.noposneg.SN1.clean3.def.tree.RelAb)
write.table(EXP.noposneg.SN1.clean3.def.tree.RelAb.otu.table,file="EXP.noposneg.SN1.clean3.def.tree.RelAb.otu.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)

EXP.noposneg.SN1.clean3.def.tree.RelAb.tax.table<-data.frame(tax_table(EXP.noposneg.SN1.clean3.def.tree.RelAb))
write.table(EXP.noposneg.SN1.clean3.def.tree.RelAb.tax.table,file="EXP.noposneg.SN1.clean3.def.tree.RelAb.tax.table.txt",sep="\t",row.names=TRUE,col.names=TRUE)



##Option1_ OTU and TAX table for top 20 OTUs at genus level

otus.top20 <- names(sort(taxa_sums(EXP.noposneg.SN1.clean3.def.tree), decreasing=TRUE))[1:20]
EXP.noposneg.SN1.clean3.def.tree.top20 <- transform_sample_counts(EXP.noposneg.SN1.clean3.def.tree, function(x) 100*x/sum(x))
EXP.noposneg.SN1.clean3.def.tree.top20 <- prune_taxa(otus.top20, EXP.noposneg.SN1.clean3.def.tree.top20)
plot_bar(EXP.noposneg.SN1.clean3.def.tree.top20, x="group", fill="Family") #you can do facetting using other parameters.



##Option2_ OTU and TAX table for top 20 OTUs at genus level

#clusterData2=psmelt(EXP.noposneg.SN1.clean3.def.tree) #already done before
#clusterData2=filter(clusterData2,Abundance>0) #already done before
#clusterAggGenus2=aggregate(Abundance ~ OTU+Genus,data=clusterData2,mean) #already done before
#clusterAggGenus2Detail=clusterAggGenus2[order(-clusterAggPhylum2$Abundance),][1:20,] #already done before
head(clusterAggGenus2Detail,20)

top20otus.clusterAggGenus2<-data.frame(clusterAggGenus2)
allotus.clusterAggGenus2<-data.frame(clusterAggGenus2)
EXP.allotus.clusterAggGenus2<-data.frame(clusterAggGenus2)
write.table(EXP.allotus.clusterAggGenus2,file="EXP.allotus.clusterAggGenus2.txt",sep="\t",row.names=TRUE,col.names=TRUE)

EXP.top20otus.clusterAggGenus2Detail<-data.frame(clusterAggGenus2Detail,20)
write.table(EXP.top20otus.clusterAggGenus2Detail,file="EXP.top20otus.clusterAggGenus2Detail.txt",sep="\t",row.names=TRUE,col.names=TRUE)
>


#End of line and script and process

###############
permanova
