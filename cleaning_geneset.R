setwd("~/Desktop/Madhesh/Analysis")
library(readr)

g1 <- read_csv("rawdata/genesets/geneset.txt")
g2 <- read_csv("rawdata/genesets/geneset-1.txt")
g3 <- read_csv("rawdata/genesets/geneset-2.txt")
g4 <- read_csv("rawdata/genesets/geneset-3.txt")
g5 <- read_csv("rawdata/genesets/geneset-4.txt")
g6 <- read_csv("rawdata/genesets/geneset-5.txt")
g7 <- read_csv("rawdata/genesets/geneset-6.txt")
g8 <- read_csv("rawdata/genesets/geneset-7.txt")
g9 <- read_csv("rawdata/genesets/geneset-8.txt")

g1 <- g1[-1,]
g2 <- g2[-1,]
g3 <- g3[-1,]
g4 <- g4[-1,]
g5 <- g5[-1,]
g6 <- g6[-1,]
g7 <- g7[-1,]
g8 <- g8[-1,]
g9 <- g9[-1,]

g1$Group <- rep(colnames(g1))
g2$Group <- rep(colnames(g2))
g3$Group <- rep(colnames(g3))
g4$Group <- rep(colnames(g4))
g5$Group <- rep(colnames(g5))
g6$Group <- rep(colnames(g6))
g7$Group <- rep(colnames(g7))
g8$Group <- rep(colnames(g8))
g9$Group <- rep(colnames(g9))

colnames(g1) <- c("Gene", "Group")
colnames(g2) <- c("Gene", "Group")
colnames(g3) <- c("Gene", "Group")
colnames(g4) <- c("Gene", "Group")
colnames(g5) <- c("Gene", "Group")
colnames(g6) <- c("Gene", "Group")
colnames(g7) <- c("Gene", "Group")
colnames(g8) <- c("Gene", "Group")
colnames(g9) <- c("Gene", "Group")

genes <- rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9)
table(genes$Group)
#write.csv(genes, "cleaned/GENESETS.csv", row.names = F)

#### Finding probe annotation for genesets #####
genes <- read.csv("cleaned/GENESETS.csv")
annotation <- read_delim("rawdata/GPL17586-45144.txt", 
                         "\t", escape_double = FALSE, comment = "#", 
                         trim_ws = TRUE)
annotation <- dplyr::filter(annotation, `locus type` == "Coding")
table(annotation$category)
v1 <- unique(as.character(genes$Gene))
annotation2 <- annotation[grep(paste(v1, collapse="|"), annotation$gene_assignment),]
annotation3 <- annotation2[,c("ID", "gene_assignment")]

library(dplyr)
library(tidyr)
require(reshape)

annotation3 <- transform(annotation3, gene_assignment = colsplit(gene_assignment, split = "\\///"))

before <- data.frame(
  attr = c(1, 30 ,4 ,6 ), 
  type = c('foo_and_bar', 'foo_and_bar_2')
)

annotation3 %>%
  separate(type, c("foo", "bar"), "_and_")