###PCOA weighted unifrac with permanova###
#first log transform read counts from phyloseq, then calculculate distance metrics and permanova test
ps.rdp.noncontam.nopos.clean2.def.tree_log<-transform_sample_counts(ps.rdp.noncontam.nopos.clean2.def.tree, function(x) log10(x+1))
dist.w <- phyloseq::distance(ps.rdp.noncontam.nopos.clean2.def.tree_log, method = "wunifrac") #wunifrac, unifrac, bray, jaccard
permanova <- adonis2(dist.w ~ RPWG, data = metadatanovaV23_2.2_URI_N, permutations = 10000, na.action = na.omit)
permanova
ordPCoA=ordinate(ps.rdp.noncontam.nopos.clean2.def.tree_log,method="PCoA",distance="unifrac",weighted=TRUE)
ordPCoA=ordinate(alretosidibgi16S.noposneg.SN6.clean2.def.tree.RelAb,method="PCoA",distance="unifrac",weighted=TRUE)


##NICE PLOTS##
bgb<-round(permanova$`Pr(>F)`[c(1)],3)
agb<-round(permanova$R2[c(1)],3)*100

plot_ordination(ps.betadiversity_woSamplesinOnly1Visit_V2and3,ordPCoA,color="RPWG")+
  guides(color=guide_legend(title="Rapid postnatal weight gain"))+
  ggtitle("Rapid postnatal weight gain (Visits 2 and 3)", "Weighted UniFrac")+
  theme_test()+theme(legend.title = element_blank()) +
  geom_vline(xintercept = 0,  linetype = "dashed", colour="grey") +
  geom_hline(yintercept = 0,  linetype = "dashed", colour = "grey")+
  ylim(-0.50,0.50)+ xlim(-1,1)+
  stat_ellipse(type="t", level= 0.85, alpha = 0.75)+
  geom_point()+coord_fixed()+ 
  annotate ("text",hjust = "left",size=3, x = -1, y = 0.45,parse=TRUE,  label=paste("R^2==", agb,"~'%'"))+ 
  annotate ("text",hjust = "left",size=3, x = -1, y = 0.40,parse=TRUE,  label =paste("italic(P)==", bgb,"~'*'"))

###Alternative plot with centroids and vectors
scrs <- ordPCoA$vectors
scrs <- cbind(as.data.frame(scrs), VAP = sample_data(subset_samples(phyloseq_RUN1_RUN2.tree.alldata_bql_t1and2_seleccio2_log,Tiempo=="T2"))$navm_confirm)
cent <- aggregate(cbind(Axis.1, Axis.2) ~ VAP , data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c('VAP','oNMDS1','oNMDS2')),
              by = 'VAP', sort = FALSE)

BD_WeigthedUnifrac_BALvsANF<-ggplot(scrs, aes(x = Axis.1, y = Axis.2, colour = VAP)) +
  geom_segment(data = segs,
               mapping = aes(xend = oNMDS1, yend = oNMDS2), alpha=0.5) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()+ scale_color_brewer(palette = "Set1")+  ggtitle("Weighted UniFrac") +
  theme_test() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey", alpha=0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey", alpha=0.2) +
  ylim(-0.30, 0.40) + xlim(-0.4, 0.3) +
  annotate("text", hjust = "left", size = 3, x = -0.4, y = 0.4, parse = TRUE, label = paste("R^2==", agb, "~'%'")) +
  annotate("text", hjust = "left", size = 3, x = -0.4, y = 0.35, parse = TRUE, label = paste("italic(P)==", bgb))
ggsave("BD_WeigthedUnifrac_BQL_NAVM.png", width = 8, height = 6, dpi = 300)
