###PCA Aitchison distances with permanova###

tse<-makeTreeSummarizedExperimentFromPhyloseq(ps2.species)
tse <- transformCounts(tse, method = "clr", pseudocount = 1)
clr_assay<-assays(tse)$clr
clr_assay <- t(clr_assay)
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")
euclidean_pcoa <- ecodist::pco(euclidean_dist)
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1],
                                pcoa2 = euclidean_pcoa$vectors[,2])
euclidean_patient_status_pcoa_df <- cbind(euclidean_pcoa_df,
                                          Group = colData(tse)$Group)



permanova <- adonis2(euclidean_dist ~ Group, data = metadata, permutations = 10000)
permanova


PCAaitchison<-ggplot(data = euclidean_patient_status_pcoa_df,
aes(x=pcoa1, y=pcoa2,
color = Group)) +
geom_point() +
labs(x = "PC1",
y = "PC2",
title = "PCA Aitchison distances", subtitle = "PERMANOVA: R2 = 0.034, p = 0.999") +
theme(title = element_text(size = 12)) +ggplot2::stat_ellipse(
ggplot2::aes(colour = Group))

###IMPROVED PLOTS###

bgb<-round(permanova$`Pr(>F)`[c(1)],3)
agb<-round(permanova$R2[c(1)],5)*100


betadiversity_aitchison_RPWG_visits2and3_wosamplesrepeated<-ggplot(data = euclidean_patient_status_pcoa_df,
aes(x=pcoa1, y=pcoa2,
color = Group)) +
geom_point() +
labs(x = "PC1",
y = "PC2",
title = "Rapid postnatal weight gain (Visits 2 and 3)", subtitle = "Aitchison distances") +
theme(title = element_text(size = 12))+
annotate ("text",hjust = "left",size=3, x = -0.02, y = 0.02,parse=TRUE,  label=paste("R^2==", agb,"~'%'"))+
annotate ("text",hjust = "left",size=3, x = -0.02, y = 0.018,parse=TRUE,  label =paste("italic(P)==", bgb,"~'*'"))+
geom_vline(xintercept = 0,  linetype = "dashed", colour="grey") +
geom_hline(yintercept = 0,  linetype = "dashed", colour = "grey")+ylim(-0.02,0.02)+ xlim(-0.02,0.02)+
stat_ellipse( alpha = 0.75)+geom_point()+coord_fixed()+theme(legend.title = element_blank())+
theme_test()+theme(legend.title = element_blank())
betadiversity_aitchison_RPWG_visits2and3_wosamplesrepeated



###PERMANOVA coefficients###

tse<-makeTreeSummarizedExperimentFromPhyloseq(taxa.phylum)
tse <- transformCounts(tse, method = "clr", pseudocount = 1)
tse <- transformCounts(tse, method = "relabundance")


rel_abund_assay <- assays(tse)$relabundance

# again transpose it to get taxa to columns
rel_abund_assay <- t(rel_abund_assay)

prova<-colData(tse)
prova$Group<-factor(prova$Group,levels=c("2BreadA1", "3breadB1", "1Bread0"))###Posar la refer?ncia com a ultim level

# then we can perform the method
permanova_cohort <- vegan::adonis(rel_abund_assay ~ Group,
                                  data = prova,
                                  permutations = 9999)

# we can obtain a the p value for our predictor:
print(paste0("Different different cohorts and variance of abundance ",
             "between samples, p-value: ", 
             as.data.frame(permanova_cohort$aov.tab)["Group", "Pr(>F)"]))




coef <- coefficients(permanova_cohort)["Group1",]

# Gets the highest coefficients
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

# Plots the coefficients
top_taxa_coeffient_plot <- ggplot(data.frame(x = top.coef,
                                             y = factor(names(top.coef),
                                                        unique(names(top.coef)))),
                                  aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="", y="", title="Coefficients 2BreadA1 at Phylum level") +
  theme_bw()

top_taxa_coeffient_plot



names <- rowData(tse)[names(top.coef), ][,"Phylum"]

# Adds new labels to the plot
top_taxa_coeffient_plot <- top_taxa_coeffient_plot +
  scale_y_discrete(labels = names) # Adds new labels
top_taxa_coeffient_plot