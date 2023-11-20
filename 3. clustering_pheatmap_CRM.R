
############################
# adding clinical data
############################

# clinical data...
run <- read.csv("data/survival_SIRT1.csv") # data from E-MTAB-7264
head(run) # os, meaning overall survival...
# X time os grade subtype
# 1 GPX_CN002_HuGene2.CEL   NA  0    G2      E1
# 2 GPX_CN003_HuGene2.CEL  140  0    G2      E1
# 3 GPX_CN004_HuGene2.CEL   NA NA  <NA>    <NA>
# 4 GPX_CN005_HuGene2.CEL   NA NA  <NA>    <NA>
# 5 GPX_CN006_HuGene2.CEL   NA NA  <NA>    <NA>
# 6 GPX_CN007_HuGene2.CEL   NA NA  <NA>    <NA>

colnames(run)[1] <- "sample"
library(tibble)
a <- as.data.frame(a)
a2 <- rownames_to_column(a)
a3 <- dplyr::left_join(a2, run, by=c("rowname"="sample"))
head(a3)
a3 <- column_to_rownames(a3, "rowname")

a3$Grade <- a3$grade
a3$grade <- NULL
a3$Subtype <- a3$subtype
a3$subtype <- NULL

###############################
# adding clustering information
###############################
res.km$cluster
a4 <- cbind(a3, res.km$cluster) # as the samples order matched, simple cbind was conducted. 
a4$cluster <- a4$'res.km$cluster' 
a4$'res.km$cluster' <- NULL
a4$Cluster <- ifelse(a4$cluster==1, 2, 1) # changed the number of clusters
o <- order(a4[,ncol(a4)]) # order the last column (clutering column)
a4 <- a4[o,] # order the matrix according to the order of the last column

#################################
# heatmap generation
#################################

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(seriation)
library(dendsort)

a5 <- t(a4[,1:29]) # 29 - the number of genes of the geneset.
dim(a5) # 29 164

a4$Cluster <- as.character(a4$Cluster)
a4$Subtype[is.na(a4$Subtype)] <- "N/A"
a4$Grade[is.na(a4$Grade)] <- "N/A"
mydf <- dplyr::select(a4, Subtype, Grade, Cluster)

my_color = list(
  Grade = c(benign="#00ff7f", G1="#00ffff", G2="#007fff", G3="#0000ff", dedifferentiated="#7f00ff", "N/A"="#b7b7b7"),
  Subtype = c(E1="#ffff00", E2="#ff0000", "N/A"="#b7b7b7"),
  Cluster = c("1"="#007FFF", "2"="#FF0000")
)
col<- colorRampPalette(c("#007FFF", "grey", "red"))(50)

quantile_breaks <- function(xs, n = 50) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(a5, n = 51)


callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(
  mat               = a5,
  color             = col, # bluered(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
#  gaps_row=c(17),
  gaps_col=c(132),
#  cutree_rows=2,
  show_colnames     = FALSE,
  # show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  #  main              = "Quantile Color Scale",
  cluster_rows = T, cluster_cols = F,
  clustering_method = "complete",
  annotation_col = mydf,
  annotation_colors = my_color,
  clustering_callback = callback
#  labels_row = as.expression(newnames)
)


