

############################
# Z scoring normalization
############################
Znorm <- t( apply(expression, 1 , function(x) scale( x , center = T, scale = T)) ) # Znorm
colnames(Znorm) <- colnames(expression)

a <- t(Znorm)
a[1:3,1:3]
# NMNAT1      NAXE        KMO
# GPX_CN002_HuGene2.CEL  0.03115502 0.4325783 -0.9127581
# GPX_CN003_HuGene2.CEL  2.12713308 0.9187591 -0.6874679
# GPX_CN004_HuGene2.CEL -2.14317557 0.6721825  0.6900182


library(ggpubr)
library(factoextra)

###############################
# Plot of the Silhouette method
###############################
fviz_nbclust(a3[,1:29], kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

#########################
# K means clustering plot
#########################
res.km <- kmeans(scale(a), 2, nstart = 25)
res.km$cluster
head(res.km, 2)

test <- fviz_cluster(res.km, data = a, # (ncol(a3)-6)
           #  palette = c("#2E9FDF", "#00AFBB"), 
             geom = "point",
             ellipse.type = "convex",
             show.clust.cent = FALSE,
             ggtheme = theme_bw()
)
test # it will make a picture.

