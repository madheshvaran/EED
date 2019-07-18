setwd("C:/Users/Madheshvaran/Desktop/Analysis")
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(qgraph)
library(readr)

###### Correlation of 51 genes between High and Low LM#################

participant <- read_excel("rawdata/participant.xlsx")
str(participant)
participant <- participant[, c(1, 3)]
participant$category <- participant$`%L`
participant <- participant[order(participant$category),]
participant <- rbind (head(participant, 30), tail(participant, 30))
participant[1:30, "category"] <- c(rep("Low_LM"))
participant[31:60, "category"] <- c(rep("High_LM"))
participant$`%L` <- NULL
str(participant)
table(participant$category)

G51 <- read_excel("rawdata/genelist_51.xlsx")
str(G51)
G51 <- G51[, c(2,3)]
colnames(G51) <- c("Gene", "Probe")

transcriptome_all <- read_excel("rawdata/Transcriptome_all60.xlsx")
str(transcriptome_all)
colnames(transcriptome_all)[1] <- "Probe"
transcriptome <- melt(transcriptome_all)
transcriptome <- filter(transcriptome, Probe %in% G51$Probe)
transcriptome <- merge(G51, transcriptome)
colnames(transcriptome)[3] <- "ID"
transcriptome <- merge(transcriptome, participant)

transcriptome <- dcast(transcriptome, ID+category ~ Gene)

LMhigh <- dplyr::filter(transcriptome, category == "High_LM")
LMlow <- dplyr::filter(transcriptome, category == "Low_LM")
LMhigh <- LMhigh[,-c(1,2)]
LMlow <- LMlow[,-c(1,2)]

LMhighcor <- cor(LMhigh)
qgraph(LMhighcor, 
       vsize=3,
       labels=TRUE, 
       theme ="classic",
       legend.cex=0.5,
       legend=TRUE,
       #layout="groups", 
       layout="spring",
       #layout="circle",
       legend.mode="nodeNames",
       nodeNames=T,
       minimum=0.2,
       cut=0.6,
       posCol= "darkgreen", 
       negCol="red")

LMlowcor <- cor(LMlow)
qgraph(LMlowcor, 
       vsize=3,
       labels=TRUE, 
       theme ="classic",
       legend.cex=0.5,
       legend=TRUE,
       #layout="groups", 
       layout="spring",
       #layout="circle",
       legend.mode="nodeNames",
       nodeNames=T,
       minimum=0.2,
       cut=0.6,
       posCol= "darkgreen", 
       negCol="red")

###### Correlation of TLR4 with interesting gene lists ####################

genes <- read.csv("cleaned/GENES_PROBES.csv")
genes$Probe <- paste(genes$Probe, 1, sep = ".")
TLR4 <- dplyr::filter(genes, Probe == "TC09000601.hg.1")
TLR4$Gene <- NULL
TLR4$Group <- NULL

### first let us consider intestinal gene sets ####

table(genes$Group)
genes1 <- dplyr::filter(genes, Group %in% c("GO_INTESTINAL_EPITHELIAL_CELL_DEVELOPMENT","GO_INTESTINAL_EPITHELIAL_CELL_DIFFERENTIATION") )
id <- transcriptome_all[,c(1)]
genes1 <- id[sample(nrow(id), 50), ]
str(genes1)
genes1 <- genes1[order(genes1$Group),]
genes1 <- rbind(genes1, TLR4)

order_i_want <- unique(as.character(genes1$Probe))

# pick out these genes from transcriptome_all

toanalyze <- dplyr::filter(transcriptome_all, Probe %in% genes1$Probe)
toanalyze <- melt(toanalyze)
str(toanalyze)
toanalyze <- toanalyze[order(factor(toanalyze$Probe, levels = order_i_want)),]
str(toanalyze)
toanalyze$Probe <- factor(toanalyze$Probe, levels = order_i_want)
str(genes1)
colnames(toanalyze) <- c("Probe", "ID", "value")
toanalyze <- merge(toanalyze, participant, by = "ID")
toanalyze <- dcast(toanalyze, ID+category ~ Probe)
TAhighLM <- dplyr::filter(toanalyze, category == "High_LM")
TAlowLM <- dplyr::filter(toanalyze, category == "Low_LM")

GROUPS <- as.character(genes1$Group)

library(Hmisc)
TAhighLM <- rcorr(as.matrix((TAhighLM[,-c(1:2)])))
str(TAhighLM)
GROUPS <- data.frame(Probe = rownames(TAhighLM))
GROUPS <- merge(GROUPS, genes1)
GROUPS$Probe <- factor(GROUPS$Probe, levels = rownames(TAhighLM))
GROUPS <- GROUPS[order(GROUPS$Probe),]
str(GROUPS)
TAlowLM <- rcorr(as.matrix((TAlowLM[,-c(1:2)])))
str(TAlowLM)
TAhighLM1 <- melt(TAhighLM$r)
TAhighLM2 <- melt(TAhighLM$P)
colnames(TAhighLM1)[3] <- "R_value"
colnames(TAhighLM2)[3] <- "P_value"

TAhighLM1 <- TAhighLM1[which(TAhighLM1$Var1 == "TC09000601.hg.1"),]
TAhighLM1 <- TAhighLM1[which(TAhighLM1$Var2 != "TC09000601.hg.1"),]
TAhighLM2 <- TAhighLM2[which(TAhighLM2$Var1 == "TC09000601.hg.1"),]
TAhighLM2 <- TAhighLM2[which(TAhighLM2$Var2 != "TC09000601.hg.1"),]
TAhighLM1$Var1 <- NULL
TAhighLM2$Var1 <- NULL
TAhighLM <- merge(TAhighLM1,TAhighLM2, by = "Var2")
TAhighLM$category <- rep("LM_high")

TAlowLM1 <- melt(TAlowLM$r)
TAlowLM2 <- melt(TAlowLM$P)
colnames(TAlowLM1)[3] <- "R_value"
colnames(TAlowLM2)[3] <- "P_value"

TAlowLM1 <- TAlowLM1[which(TAlowLM1$Var1 == "TC09000601.hg.1"),]
TAlowLM1 <- TAlowLM1[which(TAlowLM1$Var2 != "TC09000601.hg.1"),]
TAlowLM2 <- TAlowLM2[which(TAlowLM2$Var1 == "TC09000601.hg.1"),]
TAlowLM2 <- TAlowLM2[which(TAlowLM2$Var2 != "TC09000601.hg.1"),]
TAlowLM1$Var1 <- NULL
TAlowLM2$Var1 <- NULL
TAlowLM <- merge(TAlowLM1,TAlowLM2, by = "Var2")
TAlowLM$category <- rep("LM_low")

TAlowLM$P_value <- p.adjust(TAlowLM$P_value, "fdr", nrow(TAlowLM))
TAhighLM$P_value <- p.adjust(TAhighLM$P_value, "fdr", nrow(TAhighLM))

TAlowLM[which(TAlowLM$P_value < 0.3), "P_value"] <- rep("Significant")
TAlowLM[which(TAlowLM$P_value != "Significant"), "P_value"] <- rep("Not-Significant")
colnames(TAlowLM) <- c("Probe", "R_lowLM", "P_lowLM")
TAlowLM[4] <- NULL

TAhighLM[which(TAhighLM$P_value < 0.3), "P_value"] <- rep("Significant")
TAhighLM[which(TAhighLM$P_value != "Significant"), "P_value"] <- rep("Not-Significant")
colnames(TAhighLM) <- c("Probe", "R_highLM", "P_highLM")
TAhighLM[4] <- NULL

TA <- merge(TAlowLM, TAhighLM)
ggplot(TA, aes(R_highLM, R_lowLM, shape = P_highLM, colour = P_highLM))+
        geom_point()

TA <- merge(TA, genes, by.x = "Probe", by.y = "Probe", all.x = T)
TA$Group <- NULL
TA <- TA[!duplicated(TA),]

x=TA[which(TA$P_highLM == "Significant"),]
write.csv(x, "Intestinal Gene Sets.csv", row.names = F)

################ Continue the same code for the remaining gene sets ###########################

##### Correlation network for significantly TLR-4 correlated genes in the genesets ############

######## Lets start with Apoptotic Signalling Gene sets ###########

genes <- read.csv("cleaned/GENES_PROBES.csv")
genes$Probe <- paste(genes$Probe, 1, sep = ".")
TLR4 <- dplyr::filter(genes, Probe == "TC09000601.hg.1")
TLR4 <- TLR4[,-c(3)]

G <- read.csv("results/Data/Apoptosis.csv")
str(G)
G <- G[, c(6,1)]
G <- rbind(G,TLR4)

transcriptome <- melt(transcriptome_all)
transcriptome <- filter(transcriptome, Probe %in% G51$Probe)
transcriptome <- merge(G, transcriptome)
colnames(transcriptome)[3] <- "ID"
transcriptome <- merge(transcriptome, participant)
transcriptome <- transcriptome[,-c(2)]

transcriptome <- dcast(transcriptome, ID+category ~ Gene, mean)

LMhigh <- dplyr::filter(transcriptome, category == "High_LM")
LMlow <- dplyr::filter(transcriptome, category == "Low_LM")
LMhigh <- LMhigh[,-c(1,2)]
LMlow <- LMlow[,-c(1,2)]

LMhighcor <- cor(LMhigh)
write.csv(LMhighcor,"Apoptosis Signalling High.csv")
LMlowcor <- cor(LMlow)
write.csv(LMlowcor,"Apoptosis Signalling Low.csv")

qgraph(LMhighcor, 
       vsize=3,
       labels=TRUE, 
       theme ="classic",
       legend.cex=0.5,
       legend=TRUE,
       #layout="groups", 
       layout="spring",
       #layout="circle",
       legend.mode="nodeNames",
       nodeNames=T,
       minimum=0.2,
       cut=0.6,
       posCol= "darkgreen", 
       negCol="red")

LMlowcor <- cor(LMlow)
qgraph(LMlowcor, 
       vsize=3,
       labels=TRUE, 
       theme ="classic",
       legend.cex=0.5,
       legend=TRUE,
       #layout="groups", 
       layout="spring",
       #layout="circle",
       legend.mode="nodeNames",
       nodeNames=T,
       minimum=0.2,
       cut=0.6,
       posCol= "darkgreen", 
       negCol="red")

############ Similarly continue for other gene sets ############