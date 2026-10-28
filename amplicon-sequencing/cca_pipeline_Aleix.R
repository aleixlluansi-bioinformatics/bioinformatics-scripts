prova1<- filter_taxa(alretosidibgi16S.noposneg.SN6.clean2.def.tree.RelAb, function(x) mean(x) > 0.01, TRUE)
otutable<-vegan_otu(prova1)
vegan_otu <- function(physeq) {
OTU <- otu_table(physeq)
if (taxa_are_rows(OTU)) {
OTU <- t(OTU)
}
return(as(OTU, "matrix"))
}
otutable<-vegan_otu(prova1)
sampledf <- data.frame(sample_data(prova1))


##per decidir si RDA o CCA cal mirar DCA axis lenght (<3 RDA, >4 CCA, entre 3 i 4 RDA o CCA)
dca <- decorana(otutable)
dca
summary(dca)



##assignar cada variable si factor, character, integer or numeric
##Si tenim factors i volem passar a numeric o integer: as.numeric(as.character(...)).
#Al tanto amb les comes com a decimals. Han de ser punts! Sinó no ens passara be a numeric/integer. Si tenim comes: as.numeric(gsub(",", ".", as.character(...)))
str(sampledf)
sampledf$Age<-as.integer(sampledf$Age)
sampledf$Weight<-as.integer(sampledf$Weight)
sampledf$IBS.SSS<-as.numeric(sampledf$IBS.SSS)
sampledf$Calprotectin<-as.numeric(sampledf$Calprotectin)
sampledf$PCR<-as.numeric(sampledf$PCR)
summary(sampledf)



###Adonis per mirar significancia de les variables
#podem seleccionar les variables
otutable.adonis <- vegan::adonis(otutable ~ Group+Time+Bread+Age+Sex+Weight+Smoker+Concomitant.medication+Mesalazine.treatment+Biologic.treatment+Mayo+RomaIV+Depression+Adherence.to.Mediterranean.diet+Calprotectin+PCR, data=sampledf)
#o posarles totes
otutable.adonis <- vegan::adonis(otutable ~ ., data=sampledf)
#En el meu cas selecciono perquè hi ha una variable factor amb només un grup (Ansietat - Normal) i per tant no em deixa fer adonis
#si tenim NAs i ens dona error (non-conformable arrows) podem fer el seguent:
k <- complete.cases (sampledf)
otutable.adonis <- vegan::adonis(otutable[k,] ~ Weight, data=sampledf[k,])
#si fem aixo, en propers passos tambe haurem dafegir [k,] despres de otutable i sampledf en les comandes (eval, fit1 i scrs<-cbinf)


#ara triem les variables significatives
bestEnvVariables<-rownames(otutable.adonis$aov.tab)[otutable.adonis$aov.tab$"Pr(>F)"<=0.05]
#podem triar significancia 
bestEnvVariables<-rownames(otutable.adonis$aov.tab)[otutable.adonis$aov.tab$"Pr(>F)"<=0.01]
#descartem els NA
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]

eval(parse(text=paste("ord1 <- cca(otutable ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=sampledf)",sep="")))

#fer l'anova cca et drià si hi ha alguna variable significativa
anova.cca(ord1, perm=9999)

fit1 <- envfit(ord1,sampledf[,bestEnvVariables], perm = 999, display = "lc", scaling = "sites")

#per variables numeric o integer
spp.scrs <- data.frame(scores(fit1, display = "vectors"))
pval <-  fit1$vectors$pvals
spp.scrs <- cbind(spp.scrs, Vector = rownames(spp.scrs), pval) # vector table
scrs <- as.data.frame(scores(ord1, display = "sites")) # sample table
scrs <- cbind(scrs, sampledf)
spp.scrs1<- subset(spp.scrs, pval<=0.05) #extracts relevant environment vectors from envifit


#per variables factor o character
spp.scrs2<-data.frame(scores(fit1, display = "factors"))
pval2<-fit1$factors$pvals
spp.scrs2 
spp.scrs2 <- cbind(spp.scrs2, Vector = rownames(spp.scrs2))
#amb variables factor no podem descartar les no significatives, per tant hem de mirar a goodness of fit o pval quines ho son i seleccionar-les
spp.scrs3<-subset(spp.scrs2, Vector%in% c("Smokerex","Smokerno","SmokerNO","Smokersi"))


#fer el plot

library(digest)
cca.plot <- ggplot(scrs) +
geom_point(mapping = aes(x = CCA1, y = CCA2, colour = Group), alpha = 0.8, size = 3) +
#coord_fixed() + ## need aspect ratio of 1!
geom_segment(data = spp.scrs3,
aes(x = 0, xend = CCA1*1.5, y = 0, yend = CCA2*1.5),
arrow = arrow(length = unit(0.2, "cm")), colour = "blue", size=0.5) +
geom_text(data = spp.scrs3, aes(x = CCA1*2, y = CCA2*2, label = Vector),
size = 3) +
geom_segment(data = spp.scrs,
aes(x = 0, xend = CCA1*1.5, y = 0, yend = CCA2*1.5),
arrow = arrow(length = unit(0.2, "cm")), colour = "blue", size=0.5) +
geom_text(data = spp.scrs, aes(x = CCA1*2, y = CCA2*2, label = Vector),
size = 3)+
theme_bw() +
theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=8)) +
#scale_color_manual(values = set_colors3, guide = guide_legend(ncol=1)) +
scale_color_manual(values =c("red", "blue", "green"), guide = guide_legend(ncol=1)) +
#xlim(-2,3)+
theme(legend.position="right")

cca.plot

ggsave("cca.plot.png",dpi=300)

