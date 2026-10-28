####PICRUST2 functional metagenome prediction####

#We will need a Linux with Phyton and conda environmnent installed#

#Installation of Picrust (better working on linux terminal rather than Rstudio (I got problems)):
# 1: wget https://github.com/picrust/picrust2/archive/v2.5.1.tar.gz
# 2: tar xvzf  v2.5.1.tar.gz
# 3: cd picrust2-2.5.1/ (to the directory where picrust2-2.5.1 is installed)
# 4: conda env create -f picrust2-env.yaml (create environment for picrust2)
# 5: conda activate picrust2
# 6: pip install --editable .
# 7: pytest (to test if works)


#Run Picrust2#

# First create a fasta file with sequences
ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered.refseq<-refseq(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered)
write.table(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered.refseq,file="ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered.refseq.txt",sep="\t",row.names=TRUE,col.names=TRUE)
# remember to convert to fasta format!!#

#Create biomformat file 
biomformat::write_biom(biomformat::make_biom(data = t(as.matrix(otu_table(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered)))), 
                       "/home/microbioma/Documentos/data_OvariosPoliquisticos/RawData_ANF/newfigures/expr.biom")

#go back to the terminal and run Picrust2 as follows:
picrust2_pipeline.py --study_fasta physeq_filtered_woSamplesinOnly1Visit_V2and3.genus.refseq.fasta --input expr.biom --output picrust2_out_stratified --processes 6 --stratified --remove_intermediate --verbose

#return to Rstudio, charge libraries and analyse Picrust2 results
library("data.table")
library("phyloseq")
library("ALDEx2")
library("tidyverse")
library("R.utils")

#get to the directory of picrust2_out_stratified files (setwd(XXX))

picrust2 = "picrust2_out_stratified3"
list.files(picrust2, recursive = TRUE)
p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

mapfile="description_mapfiles" #description_mapfiles need to be in the wd#
list.files(mapfile, recursive = TRUE)
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)


#run Aldex2 to ckeck for significances and effect sizes (data is compositional)
#mc.samples can be modified if processing is too slow
#Subset samples if needed by (where ps1a is phyloseq object already subset):
#p2EC1 = p2EC[,sample_names(ps1a)]
#p2KO1 = p2KO[,sample_names(ps1a)]
#p2PW1 = p2PW[,sample_names(ps1a)]

set.seed(12345)
system.time({
  aldex2_EC = aldex(p2EC, prova_RPwG_picrust, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)
})

hist(aldex2_EC$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")



system.time({
  aldex2_KO = aldex(p2KO, prova_RPwG_picrust, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)
})

hist(aldex2_KO$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")



system.time({
  aldex2_PW = aldex(p2PW, prova_RPwG_picrust, mc.samples = 500, test = "t", 
                    effect = TRUE, denom = "iqlr", verbose = TRUE)
})

hist(aldex2_PW$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")


##Ancombc and logfoldchange plot##

dataP2PW<-as.data.frame(t(p2PW))
metadata.v2v3<-as.data.frame(sample_data(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered))

physeq.picrust.v2v3<-phyloseq(otu_table(dataP2PW, taxa_are_rows = FALSE),sample_data(metadata.v2v3))

out <- ancombc(
phyloseq = physeq.picrust.v2v3,
formula = "RPWG",
p_adj_method = "fdr",
lib_cut = 1000,
group = "RPWG",
struc_zero = TRUE,
neg_lb = TRUE,
alpha = 0.05,
max_iter = 100,
conserve = TRUE,
global = TRUE # multi group comparison will be deactivated automatically
)
res_picrust_RPWG<-out$res
view(res_picrust)

###lot the lfc/effect size of the significant taxa###

ancombc_out_ASV_V23_RPWG_picrust <- data.frame(res_picrust_RPWG)

sigtab_RPWG_picrust<-ancombc_out_ASV_V23_RPWG_picrust[(ancombc_out_ASV_V23_RPWG_picrust$RPWGRPWG.4 < 0.05),]

view(sigtab_RPWG_picrust)

##we selected the significant taxa but we want to subset the taxa with >2 effect size##
sigtab_RPWG_picrust$PWnames<-rownames(sigtab_RPWG_picrust)
sigtab_RPWG_picrust<-subset(sigtab_RPWG_picrust, PWnames %in% c("ANAEROFRUCAT-PWY" , "ARGORNPROST-PWY","GLYCOCAT-PWY", "PRPP-PWY", "PWY-5100", "PWY-5121", "PWY-5910","PWY-6383" ,"PWY-6397","PWY-6470" , "PWY-6545", "PWY-7007"  ,"PWY-7392" ,"PWY-922" ,"PWY1G-0"))

pos_color <- "black"  # Color for positive values
neg_color <- "grey"   # Color for negative values

ancombcplot_picrust_RPWG_effectsize_ASV<-ggplot(sigtab_RPWG_picrust, aes(y = reorder(PWnames, RPWGRPWG.2), x = RPWGRPWG.2, fill = ifelse(RPWGRPWG.2 >= 0, "RPWG", "NPWG"))) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.2) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10)) +
  xlab("Effect size") + ylab("Significant pathways")+
  scale_fill_manual(values = c("RPWG" = pos_color, "NPWG" = neg_color) , name = "Postnatal weight gain") 
ancombcplot_picrust_RPWG_effectsize_ASV
ggsave("ancombcplot_picrust_RPWG_effectsize_ASV.png",dpi=300, height = 4, width = 10)


#change y=rownames(sigtab) for y=Genus in case we want genus names#
ggplot(sigtab, aes(y=rownames(sigtab), x=RPWGRPWG)) +
geom_vline(xintercept = 0.0, color = "gray", size = 0.2) +
geom_bar(stat = "identity") +
theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0, size = 10))+ ggtitle("NPWG vs RPWG (Visits 2 and 3)")+xlab("log2 fold change")+ylab("Predicted KEGG functional pathways")

ggsave("LFC_picrust_significant_genus_v2v3.png",dpi=300)

##alternative: volcano plot##

resOE_df_ordered <- ancombc_out_picrust_RPGW_V23[order(ancombc_out_picrust_RPGW_V23$RPWGRPWG.4), ] 

resOE_df_ordered$genelabels <- rownames(resOE_df_ordered) %in% rownames(resOE_df_ordered[1:3,])

View(resOE_df_ordered)

library(ggrepel)

ggplot(resOE_df_ordered) +
  geom_point(aes(x = RPWGRPWG.2, y = -log10(RPWGRPWG.4), colour = RPWGRPWG.5)) +
  geom_text_repel(aes(x = RPWGRPWG.2, y = -log10(RPWGRPWG.4), label = ifelse(genelabels == T, rownames(resOE_df_ordered),""))) +
  ggtitle("Predicted KEGG functional pathways") +
  xlab("Effect size") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +guides(color = guide_legend(title = "Differential abundance"))+ylim(0,2)


##Spearman's correlation heatmap with ASVs##

otutable.1<-data.frame(otu_table(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered))
dataP2PW.2<-subset(dataP2PW, colnames(dataP2PW) %in% c("ANAEROFRUCAT-PWY" , "ARGORNPROST-PWY","GLYCOCAT-PWY", "PRPP-PWY", "PWY-5100", "PWY-5121", "PWY-5910","PWY-6383" ,"PWY-6397","PWY-6470" , "PWY-6545", "PWY-7007"  ,"PWY-7392" ,"PWY-922" ,"PWY1G-0"))
#if not working, use this: 
#library(dplyr)
#dataP2PW.2 <- dataP2PW %>%
  #select("ANAEROFRUCAT-PWY", "ARGORNPROST-PWY", "GLYCOCAT-PWY", "PRPP-PWY", "PWY-5100", "PWY-5121", "PWY-5910", "PWY-6383", "PWY-6397", "PWY-6470", "PWY-6545", "PWY-7007", "PWY-7392", "PWY-922", "PWY1G-0")

cor_matrix <- cor(otutable.1, dataP2PW.2, method = "spearman")

p_values <- matrix(NA, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
for (i in 1:nrow(cor_matrix)) {
for (j in 1:ncol(cor_matrix)) {
result <- cor.test(otutable.1[, i], dataP2PW.2[, j], method = "spearman")
p_values[i, j] <- result$p.value
}
}

asterisk_matrix <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
asterisk_matrix[p_values < 0.05] <- "*"


taxonomy<-data.frame(tax_table(ps.rdp.noncontam.nopos.clean2.def.tree_V2V3_worepeated_filtered))
taxonomy$GenusSpeciesname<-with(taxonomy, ifelse(is.na(Genus), NA, 
                                                             ifelse(is.na(Species), paste(rownames(taxonomy), Genus, "spp."), 
                                                                    paste(rownames(taxonomy), paste(Genus, Species, sep = " ")))))
rownames(cor_matrix)<-taxonomy$GenusSpeciesname


png("heatmap_picrust_ASV.png", width = 1000, height = 600, res = 150)
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),  # Set color to NA to remove the default color scheme
         main = "Spearman's Correlation Heatmap",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_colors = list(corr = c("black")),
         annotation_legend = FALSE,
         display_numbers = asterisk_matrix)  # Use 'display_numbers' to include the asterisk matrix
dev.off()



