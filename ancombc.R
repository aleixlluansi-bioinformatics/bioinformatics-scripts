library(ANCOMBC)
out <- ancombc(
phyloseq = merged_physeq.2.genus,
formula = "Tipo_muestra",
p_adj_method = "fdr",
lib_cut = 0,
group = "Tipo_muestra",
struc_zero = TRUE,
neg_lb = TRUE,
alpha = 0.05,
max_iter = 100,
conserve = TRUE,
global = TRUE # multi group comparison will be deactivated automatically
)
res3<-out$res
view(res3)
rownames(res_ATB)<-tax_table(physeq_beta)$Genus
rownames(res)<-data.frame(tax_table(physeq_beta))$Genus
data.frame(tax_table(physeq_beta))$Genus



####plot the lfc/effect size of the significant taxa####

ancombc_out_ASV_BALvsANF <- data.frame(out$res)

  sigtab_AOM<-ancombc_out_ASV_BALvsANF[(ancombc_out_ASV_BALvsANF$Tipo_muestraBAL.3 < 0.05),]
sigtab_AOM = cbind(as(sigtab_AOM,"data.frame"), as(tax_table(merged_physeq.2.genus)[rownames(sigtab_AOM),],"matrix"))
sigtabgen_AOM<-subset(sigtab_AOM,!is.na(Genus))
sigtabgen_AOM$GenusSpeciesname<-with(sigtabgen_AOM, ifelse(is.na(Genus), NA, 
  ifelse(is.na(Species), paste(rownames(sigtabgen_AOM), Genus, "spp."), 
         paste(rownames(sigtabgen_AOM), paste(Genus, Species, sep = " ")))))
  
view(sigtabgen_AOM)

##we selected the significant taxa but we want to subset the taxa with >2 effect size##
sigtabgen_AOM <- subset(sigtabgen_AOM, Tipo_muestraBAL.2 > 1.9 | Tipo_muestraBAL.2 < -1.9)

#easy plot
ancombcplot_RPWG_effectsize_ASV<-ggplot(sigtabgen_RPWG, aes(y=Genus, x=RPWGRPWG.2)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.2) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0, size = 10))+ ggtitle("NPWG - RPWG (Visits 2 and 3)")+xlab("Effect size")

#now with colours
pos_color <- "#fb5858"  # Color for positive values
neg_color <- "darkorchid"   # Color for negative values

pos_color <- "black"  # Color for positive values
neg_color <- "grey"   # Color for negative values

ancombcplot_BQL_BALvsANF_genus<-ggplot(sigtabgen_AOM, aes(y = reorder(Genus, Tipo_muestraBAL.2), x = Tipo_muestraBAL.2, fill = ifelse(Tipo_muestraBAL.2 >= 0, "BAL", "ANF"))) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.2) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0, size = 10)) +
   xlab("Effect size") + ylab("Significant taxa")+
  scale_fill_manual(values = c("BAL" = pos_color, "ANF" = neg_color), labels= c("NPA", "BAL") , name = "Type of sample") 


ancombcplot_BQL_BALvsANF_genus
ggsave("ancombcplot_BQL_BALvsANF_genus.png",dpi=600, height = 4, width = 10)
