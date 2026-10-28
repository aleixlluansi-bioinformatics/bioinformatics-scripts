### Estimate the optimal number of clusters by using silhouette

# Make the hierarchical clustering using the bray dissimilarity matrix previously calculated
HC_bray = hclust(dist.w23_2.nou2, method="average")

# Calculate silhouette index for each agrupation. The optimal number of clusters is 11.
#pdf("silhouette_number_groups.pdf", height=5, width=7)
plot(2:15, sapply(2:15, function(i) {
  mean(silhouette(cutree(HC_bray, i), dmatrix=as.matrix(dist.w23_2.nou2))[,"sil_width"]) }),
  xlab="Number of clusters", ylab="Average Silhouette", type="b", pch=20)
#dev.off()

### Add the clustering to the phyloseq object
clusterCut_bray <- cutree(HC_bray, 2)  
table<-table(clusterCut_bray, metadatanovaV23_2.2$RPWG4)
all.equal(names(clusterCut_bray), rownames(sample_data(ps.rdp.noncontam.nopos.clean2.def.tree_V23_filtered2_Genus)))
sample_data(ps.rdp.noncontam.nopos.clean2.def.tree_V23_filtered2_Genus)$cluster_bray<-clusterCut_bray
sample_data(ps.rdp.noncontam.nopos.clean2.def.tree_V23_filtered2_Genus)$cluster_bray<-as.factor(sample_data(ps.rdp.noncontam.nopos.clean2.def.tree_V23_filtered2_Genus)$cluster_bray)
plot_ordination(ps.rdp.noncontam.nopos.clean2.def.tree_V23_filtered2_Genus,ordPCoA,color="cluster_bray")+
  stat_ellipse(type="t", level= 0.85, alpha = 0.75)
