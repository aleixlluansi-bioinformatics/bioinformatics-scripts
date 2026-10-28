# ============================================================
# Tutorial on drawing an NMDS plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(vegan);library(ggplot2)
abund_table<-abund_table(alretosidibgi16S.noposneg.SN6.clean2.def.tree)
meta_table<-sample_data(alretosidibgi16S.noposneg.SN6.clean2.def.tree)
#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),meta_table)
head(grouping_info)
#Get MDS stats
sol<-metaMDS(abund_table,distance = "bray", k = 2, trymax = 50)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Group=as.factor(grouping_info[,1]),Time=as.factor(grouping_info[,2]),Sex=as.factor(grouping_info[,5]))
head(sol)
plot.new()
ord<-ordiellipse(sol, as.factor(grouping_info[,1]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()
#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level 
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
+ {
+   theta <- (0:npoints) * 2 * pi/npoints
+   Circle <- cbind(cos(theta), sin(theta))
+   t(center + scale * t(Circle %*% chol(cov)))
+ }


#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Group)){
+   if(g!="" && (g %in% names(ord))){
+  
+     df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Group==g,],
+                                                      veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
+                                   ,Group=g))
+   }
+ }

# head(df_ell)
head(df_ell)

#Generate mean values from NMDS plot grouped on Countries
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Group),mean)
NMDS.mean

#Now do the actual plotting
p<-ggplot(data=NMDS,aes(x,y,colour=Group))
p
p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
p
p<-p+geom_point(aes(shape=Time))+theme_bw()
p
pdf("NMDS.pdf")
print(p)
dev.off()
