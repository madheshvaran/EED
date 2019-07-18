library(biomaRt)

genes_df <- read.csv("cleaned/GENESETS.csv")
genes_v <- unique(as.character(genes$Gene))

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
genes <- getBM(attributes = c("affy_hta_2_0", "ensembl_gene_id","ensembl_transcript_id","external_gene_name","entrezgene", "gene_biotype","transcript_biotype", "go_id","strand","chromosome_name","transcript_length","description"),
               filters = "external_gene_name",
               values = genes_v,
               mart = ensembl)

table(genes$gene_biotype)
genes <- dplyr::filter(genes, gene_biotype == "protein_coding")
table(genes$transcript_biotype)
genes <- dplyr::filter(genes, transcript_biotype == "protein_coding")
colnames(genes)
genes <- dplyr::select(genes, affy_hta_2_0, ensembl_gene_id, external_gene_name)
genes <- dplyr::filter(genes, affy_hta_2_0 != "")
genes <- genes[!duplicated(genes),]
genes <- dplyr::select(genes, affy_hta_2_0, external_gene_name)
genes <- genes[!duplicated(genes),]
colnames(genes) <- c("Probe", "Gene")
genes_merged <- merge(genes, genes_df)
length(unique(genes_merged$Gene))
#write.csv(genes_merged, "cleaned/GENES_PROBES.csv", row.names = F)
