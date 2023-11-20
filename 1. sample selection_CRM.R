
#################################################
# uploading matching data between official gene name and affy number of the NAD pathway
#################################################

# official gene names of NAD pathway
nad <- read.csv("data/nadgene_all.csv") 
head(nad)

# gene kegg go_all salvage denovo
# 1 ACMSD    0      1       0      1
# 2 AFMID    0      1       0      1
# 3 ASPDH    0      1       0      0
# 4  HAAO    0      1       0      1
# 5  IDO1    0      1       0      1
# 6  IDO2    0      1       0      1

# matching data between official gene name and affy number
nadaffy <- read.csv("data/nadgene_affy.csv") 
head(nadaffy)

# gene     affy
# 1 NNMT 16731464
# 2 NNMT 16731465
# 3 NNMT 16731462
# 4 NNMT  3349858
# 5 NNMT 16731463
# 6 NNMT  3349859

library(dplyr)
nadaffy <- dplyr::left_join(nad, nadaffy, by="gene")
head(nadaffy)

# gene kegg go_all salvage denovo     affy
# 1 ACMSD    0      1       0      1  2507187
# 2 ACMSD    0      1       0      1  8045378
# 3 ACMSD    0      1       0      1 16885897
# 4 ACMSD    0      1       0      1  2507188
# 5 ACMSD    0      1       0      1  8045377
# 6 ACMSD    0      1       0      1  2507185


a <- subset(nadaffy, go_all==1) # go_all, meaning all genes of GO NAD pathway


#################################################
# bringing reference data
#################################################

exp_palmieri <- readRDS("data/exp_palmieri.Rds") # mRNA data from E-MTAB-7264
exp_palmieri[1:3,1:3]

# GPX_CN002_HuGene2.CEL GPX_CN003_HuGene2.CEL GPX_CN004_HuGene2.CEL
# 16650001              2.909333              2.695944              2.672546
# 16650003              3.063851              3.375709              3.265057
# 16650005              5.493629              5.103422              5.517001

# select data only for NAD pathway
b <- subset(exp_palmieri, rownames(exp_palmieri) %in% a$affy)
dim(b) # [1]  29 164


#################################################
# changing affy number to official gene name
#################################################
library(tibble)
b <- as.data.frame(b)
b1 <- rownames_to_column(b)
nadaffy2 <- dplyr::select(nadaffy, gene, affy)
nadaffy2$affy <- as.character(nadaffy2$affy)
b2 <- dplyr::left_join(b1, nadaffy2, by=c("rowname"="affy"))
b2$gene <- make.names(b2$gene, unique=T)
b2 <- column_to_rownames(b2, "gene")
head(b2)
b2$rowname <- NULL
expression <- b2
nrow(expression) # 29

head(expression)[1:5, 1:5]
# GPX_CN002_HuGene2.CEL GPX_CN003_HuGene2.CEL GPX_CN004_HuGene2.CEL GPX_CN005_HuGene2.CEL GPX_CN006_HuGene2.CEL
# NMNAT1              6.113283              6.655344              5.550958              6.451859              6.501470
# NADK                7.381978              7.131490              7.007270              7.311928              7.252173
# NT5C1A              4.843266              5.078110              4.495855              4.927131              4.896955
# NMNAT2              5.198799              5.020943              6.370806              5.042517              5.651878
# NAMPT               8.880582              9.028991             10.528145              9.012381              9.239500

