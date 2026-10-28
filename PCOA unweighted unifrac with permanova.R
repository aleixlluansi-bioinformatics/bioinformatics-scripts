###PCOA unweighted unifrac with permanova###


unifrac.dist <- UniFrac(alretosidibgi16S.noposneg.SN6.clean2.def.tree.RelAb,
weighted = FALSE,
normalized = TRUE,
parallel = FALSE,
fast = TRUE)

permanova <- adonis2(unifrac.dist ~ Group, data = metadata, permutations = 10000)
permanova
pcoaweigtedunifrac<-plot_ordination(alretosidibgi16S.noposneg.SN6.clean2.def.tree.RelAb,ordPCoA,color="Group2") +guides(color=guide_legend(title="Group"))+ggtitle("PCOA unweighted unifrac", "PERMANOVA: R2 = 0.053, p = 0.241")

