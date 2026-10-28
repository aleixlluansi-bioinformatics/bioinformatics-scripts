#####CODA4MICROBIOME##########

library(coda4microbiome)



####Without using Core 10 taxa###

#subset phyloseq if needed
taxa.genus.BreadAvsBreadB<-subset_samples(taxa.genus,Time!="0")

#Create otutable and change zeros for ones
otutable.taxa.genus.BreadAvsBreadB<-data.frame(otu_table(taxa.genus.BreadAvsBreadB))
otutable.genus.BreadAvsBreadB[otutable.genus.BreadAvsBreadB==0]=1

#Create sample data
sample_data_taxa.genus.BreadAvsBreadB <- data.frame(sample_data(taxa.genus.BreadAvsBreadB))


#Coda4microbiome function
coda_glmnet_genus1<-coda_glmnet(x=otutable.genus.BreadAvsBreadB,y=sample_data_taxa.genus.BreadAvsBreadB$Time)





###Using the Core10 taxa###

#subset phyloseq if needed
taxa.genus.BreadAvsBreadB<-subset_samples(taxa.genus,Time!="0")

#Create otutable transposed, add column Taxon, select Core10 taxa, transpose core10, convert character to numeric variables and change zeros for ones.
otutable.genus.BreadAvsBreadB<-data.frame(t(otu_table(taxa.genus.BreadAvsBreadB)))
Core10.taxa <- core_members(taxa.genus.BreadAvsBreadB, prevalence = 10/100, detection = 1/100)
otutable.genus.BreadAvsBreadB$Taxon<-rownames(otutable.genus.BreadAvsBreadB)
otutable.core10.genus.BreadAvsBreadB <- otutable.genus.BreadAvsBreadB[otutable.genus.BreadAvsBreadB$Taxon%in%Core10.taxa.genus,]
otutable.core10.genus.BreadAvsBreadB<- as.data.frame(x = t(otutable.core10.genus.BreadAvsBreadB), stringsAsFactors = FALSE)
colnames(otutable.core10.genus.BreadAvsBreadB) #Recomended to copy colanmes output to notebook (for example) and add commas between names
col.names<-c( "ASV7",  "ASV181","ASV21", "ASV128","ASV14", "ASV34", "ASV156","ASV97", "ASV63", "ASV27", "ASV327","ASV149","ASV174",              "ASV133","ASV45", "ASV165","ASV682","ASV219","ASV68", "ASV22", "ASV454","ASV17", "ASV1",  "ASV192","ASV499","ASV140",              "ASV57", "ASV144","ASV382","ASV297","ASV114","ASV312","ASV117","ASV317","ASV20", "ASV429","ASV278","ASV24", "ASV341",              "ASV52", "ASV215","ASV164","ASV483","ASV118","ASV79", "ASV83", "ASV231","ASV42", "ASV13", "ASV8",  "ASV18", "ASV47",              "ASV94", "ASV127","ASV419","ASV591","ASV12", "ASV15", "ASV9",  "ASV51", "ASV84", "ASV162","ASV40", "ASV6",  "ASV26",              "ASV169","ASV529","ASV33", "ASV101","ASV36", "ASV28", "ASV105","ASV5")
otutable.genus.BreadAvsBreadB$Taxon<-NULL
otutable.core10.genus.BreadAvsBreadB[col.names] <- sapply(otutable.core10.genus.BreadAvsBreadB[col.names],as.numeric)
otutable.core10.genus.BreadAvsBreadB[otutable.core10.genus.BreadAvsBreadB==0]=1
otutable.core10.genus.BreadAvsBreadB<-otutable.core10.genus.BreadAvsBreadB[-c(47), ]  #Must look at what number of column is the "Taxon" column to remove it

#Create sample data
sample_data_taxa.genus.BreadAvsBreadB <- data.frame(sample_data(taxa.genus.BreadAvsBreadB))


#Coda4microbiome function
coda_glmnet_genus1<-coda_glmnet(x=otutable.core10.genus.BreadAvsBreadB,y=sample_data_taxa.genus.BreadAvsBreadB$Time)









