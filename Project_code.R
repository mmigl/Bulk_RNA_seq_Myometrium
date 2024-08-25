install.packages("ggfortify")
BiocManager::install('EnhancedVolcano')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(ggfortify)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(tibble)
library(tidyr)
library(clusterProfiler)

data <- read.csv("C:/Users/mehak/Downloads/bulk_rna_seq_data.csv")
View(data)
rownames(data) <- data$Sample
data <- data[-1]
View(data)

ensembl_ids_to_gene_symbols <- function(data) {
  names <- rownames(data)
  symbols <- mapIds(org.Hs.eg.db, keys = names,
                    column = 'SYMBOL', keytype = 'ENSEMBL')
  symbols[is.na(symbols)] <- "Unknown"
  unique_symbols <- make.unique(symbols)
  rownames(data) <- unique_symbols
  keep <- !grepl("^Unknown", rownames(data))
  data <- data[keep,]
  return(data)
}
processed_data <- ensembl_ids_to_gene_symbols(data)
View(processed_data)

sample_name <- names(processed_data)
conditions <- factor(c(rep("WT", ncol(processed_data) / 2), rep("G44D", ncol(processed_data) / 2)))
conditions <- relevel(conditions, ref = "WT")
colData <- data.frame(row.names = colnames(processed_data), condition = conditions)
#all(colnames(data) %in% rownames(colData))
#all(colnames(data) == rownames(colData)) 
dds <- DESeqDataSetFromMatrix(countData = processed_data,
                              colData = colData,
                              design = ~condition)
dds <- DESeq(dds)

#Quality control and assessment 

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab("PC1") +
  ylab("PC2") + 
  coord_fixed() +
  geom_text(aes(label = sample_name), vjust = 2, hjust = 0.5, check_overlap = TRUE, show.legend = FALSE) 

vst_mat <- assay(vsd)
vst_cor <- cor(vst_mat)
library(pheatmap)
pheatmap(vst_cor)

#DE analysis

resultsNames(dds)
#contrast <-  c("condition", "G44D", "WT")
#res <- results(dds, alpha = 0.05)
#res_G44D_vs_WT <- results(dds, contrast = contrast, alpha = 0.05)
res <- results(dds, contrast = c("condition", "G44D", "WT"))
res %>%
  as.data.frame() %>%
  arrange((padj), desc(log2FoldChange)) %>%
  head(n=30)

significant_genes<- res %>%
  as.data.frame() %>%
  filter(padj <=0.01, abs(log2FoldChange) >= 2) %>% 
  rownames()
significant_genes

lfc_res <- lfcShrink(dds, coef = "condition_G44D_vs_WT", res = res)
lfc_res

lfc_res_n <- lfcShrink(dds, contrast = c("condition", "G44D", "WT"), res = res, alpha = 0.05, type = "normal")
lfc_res_n

BiocManager::install("ashr")
lfc_res_a <- lfcShrink(dds, contrast = c("condition", "G44D", "WT"), res = res, alpha = 0.05, type = "ashr")
lfc_res_a

plotDispEsts(dds)
plotMA(res)
plotMA(lfc_res)
plotMA(lfc_res_n)
plotMA(lfc_res_a)  #This plot is better aligned with the pca plot

#The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DESeq2 results",
                subtitle = "Differential expression")
EnhancedVolcano(lfc_res_a,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "DESeq2 results",
                subtitle = "Differential expression")

padj.cutoff <- 0.01
lfc.cutoff <- 2
lfc_res_a_table <- lfc_res_a %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
lfc_res_a_table
sig_de <- lfc_res_a_table %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_de

#Visualizations

coldata_table <- colData %>% 
  rownames_to_column(var="samples") %>% 
  as_tibble()
View(coldata_table)

normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_table <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
normalized_counts_table

d <- plotCounts(dds, gene="MED12", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MED12") +
  theme(plot.title = element_text(hjust = 0.5))

top20_sig_de_genes <- sig_de %>% 
  arrange(padj) %>% 
  pull(gene) %>% 	
  head(n=20) 	
top20_sig_de_genes
top20_normalized_counts <- normalized_counts_table %>%
  filter(gene %in% top20_sig_de_genes)
top20_normalized_counts

plotting_data <- top20_normalized_counts %>%
  gather(key = "samples", value = "normalized_counts", -gene)
View(plotting_data)

plotting_data <- inner_join(coldata_table, plotting_data)
View(plotting_data)

ggplot(plotting_data) +
  geom_point(aes(x = gene, y = normalized_counts, color = condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

heatmap_data <- top20_normalized_counts %>% column_to_rownames(var = "gene")
heatmap_data
pheatmap(heatmap_data, 
         cluster_rows = TRUE,
         show_rownames = TRUE,
         annotation = colData, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#Over-representation analysis

install.packages("msigdbr")
library(msigdbr)
background_genes<- res %>% 
  as.data.frame() %>% 
  filter(baseMean != 0) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pull(gene)
background_genes
entrez_ids_background <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids_background
entrez_ids_de <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids_de$ENTREZID

m_df <- msigdbr(species = "Homo sapiens")

h_database <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
h_database
table(h_database$gs_name)

ora_h <- enricher(entrez_ids_de$ENTREZID, TERM2GENE=m_t2g, 
               universe = entrez_ids_background$ENTREZID)
ora_h
barplot(ora_h)

ego <- enrichGO(gene = entrez_ids_de$ENTREZID, 
                universe = entrez_ids_background$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
go_data <- data.frame(ego)
View(go_data)
dotplot(ego, showCategory=20)
barplot(ego, showCategory = 10)

kegg <- enrichKEGG(gene = entrez_ids_de$ENTREZID,
                   universe = entrez_ids_background$ENTREZID,
                   organism = 'hsa',
                   keyType = "kegg",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)
kegg_data <- data.frame(kegg)
View(kegg_data)
dotplot(kegg, showCategory=20)
barplot(kegg, showCategory = 15)

#Gene set enrichment analysis

res_df<- res %>% 
  as.data.frame() %>% 
  filter(baseMean != 0) %>%
  tibble::rownames_to_column(var = "gene")
res_df
ranked_res_df<- res_df %>% 
  mutate(signed_rank_stats = sign(log2FoldChange) * -log10(pvalue)) %>%
  left_join(entrez_ids_background, by= c("gene" = "SYMBOL")) %>%
  arrange(desc(signed_rank_stats))
ranked_res_df
gene_list<- ranked_res_df$signed_rank_stats
gene_list
names(gene_list)<- ranked_res_df$ENTREZID
gene_list
gsea_h <- GSEA(gene_list, TERM2GENE=m_t2g)
gsea_h@result
gsea_h@result %>% View()
p1 <- gseaplot(gsea_h, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
              by = "runningScore", title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

p2 <- gseaplot(gsea_h, geneSetID = "HALLMARK_E2F_TARGETS", 
               by = "runningScore", title = "HALLMARK_E2F_TARGETS")
p1/p2
dotplot(gsea_h)