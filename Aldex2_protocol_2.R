###Aldex2 protocol 2###

#subset samples - remember it is only possible comparisons between two groups#
COVID2_rarefied10000_aggfamily_HCTLvsCOVID6<-subset_samples(COVID2_rarefied10000_aggfamily,Disease.level!="COVID 12M")

#make tse from phyloseq and keep only taxa with more than 10% prevalence#
tse_rarefied_aggfamily_COVID6vsCOVID12<-makeTreeSummarizedExperimentFromPhyloseq(COVID2_rarefied10000_aggfamily_COVID6vsCOVID12)
tse_rarefied_aggfamily_COVID6vsCOVID12.10prevalent<-subsetByPrevalentTaxa(tse_rarefied_aggfamily_COVID6vsCOVID12, detection = 0, prevalence = 0.1)

#perform aldex clr transformation, ttest and effect#
x <- aldex.clr(
  reads = assay(tse_rarefied_aggfamily_COVID6vsCOVID12.10prevalent),
  conds = colData(tse_rarefied_aggfamily_COVID6vsCOVID12.10prevalent)$Disease.level, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)
# calculates expected values of the Welch's t-test and Wilcoxon rank test on
# the data returned by aldex.clr
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)
# determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

rownames_to_column(aldex_out, "Family") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  select(Family, we.eBH, wi.eBH, effect, overlap) %>%
  kable()



#####GLM####


##Choose intercept variable (in the example we select "Remission at 6M")##
sample_data(provant_coses)$Cathegory_simple <- factor(sample_data(provant_coses)$Cathegory_simple, levels = c("Remission at 6M", "No remission at 6M", "Inactivity control"))
#check it#
levels(sample_data(provant_coses)$Cathegory_simple)

#create model and run glm test#
Counts<-data.frame(t(otu_table(provant_coses)))
Data <- data.frame(sample_data(provant_coses))
mm <- model.matrix(~Cathegory_simple, Data)
x <- aldex.clr(Counts, mm, mc.samples=128, denom="all")
glm.test <- aldex.glm(x, mm)

#run glm effect#
glm.effect <- aldex.glm.effect(x)

#create aldex output#
aldex_glm_out2 <- data.frame(glm.test, glm.effect)

#visualize all results#
view(aldex_glm_out2)

#visualize significant results#
rownames_to_column(aldex_glm_out2, "Genus") %>%
  filter(model.Cathegory_simpleNo.remission.at.6M.Pr...t.. <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  select(Genus, model.Cathegory_simpleNo.remission.at.6M.Pr...t.., Cathegory_simpleNo.remission.at.6M.effect,model.Cathegory_simpleNo.remission.at.6M.Pr...t...BH)


#plot of p values and size effect#
plot(aldex_glm_out2$Cathegory_simpleNo.remission.at.6M.effect, aldex_glm_out2$model.Cathegory_simpleNo.remission.at.6M.Pr...t.., log="y", cex=0.7, col=rgb(0,0,1,0.2),
     pch=19, xlab="Effect size", ylab="P value", main="Effect size plot")
points(aldex_glm_out2$Cathegory_simpleNo.remission.at.6M.effect, aldex_glm_out2$model.Cathegory_simpleNo.remission.at.6M.Pr...t...BH, cex=0.7, col=rgb(1,0,0,0.2),
       pch=19)
abline(h=0.05, lty=2, col="grey")
legend(15,1, legend=c("P value", "BH-adjusted"), pch=19, col=c("blue", "red"))


#create plot with taxonomy assigned#
taxonomia_dadesbusquets<-data.frame(phyloseq::tax_table(provant_coses))
aldex_glm_out2$Genus<-taxonomia_dadesbusquets$Genus
xsig=aldex_glm_out2[aldex_glm_out2$model.Cathegory_simpleNo.remission.at.6M.Pr...t..<0.05,]
xsig2=arrange(xsig,desc(abs(Cathegory_simpleNo.remission.at.6M.effect)))
ggplot(xsig2,aes(x=Genus,y=Cathegory_simpleNo.remission.at.6M.effect))+geom_bar(stat="identity")+coord_flip()+theme_classic()

