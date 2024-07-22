rm(list = ls())

require(ggplot2)
require(pheatmap)

files <- list.files(path= ".", pattern="*.out.tab")
# using perl to manpulate file names by trimming file extension
labs <- paste("", gsub(".out.tab", "", files, perl=TRUE), sep="")
cov <- list()
for (i in labs) {
  filepath <- file.path(paste(i,".out.tab",sep=""))
  tmp_read_per_gene_counts <- read.table(filepath,sep = "\t", header=F, stringsAsFactors=FALSE)
  cov[[i]] <- tmp_read_per_gene_counts[c(1,2)] # pick gene id and unstranded
  colnames(cov[[i]]) <- c("ENSEMBL_GeneID", i)
}

## construct one data frame from list of data.frames using reduce function
df <-Reduce(function(x,y) merge(x = x, y = y, by ="ENSEMBL_GeneID"), cov)
#df <- df[-c(1:5),]
dim(df)
tmpdf<-df[grep("PAR", df$ENSEMBL_GeneID),]
df<-df[!df$ENSEMBL_GeneID%in%tmpdf$ENSEMBL_GeneID,]
dim(df)


rownames(df) <- df$ENSEMBL_GeneID
df$ENSEMBL_GeneID<-NULL
df_t <- t(df)
rownames(df_t)<-gsub('AReadsPerGene','',rownames(df_t))
df_t <- data.frame(df_t)
df_t$Group <- gsub('.*_','',rownames(df_t))
table(df_t$Group)

df_t$Group <- factor(df_t$Group, levels = c("V02.0", "V02a", "V05.0", "V05.a", "V05a"),
                     labels = c("V02.0", "V02a", "V05.0", "V05.0", "V05a"))
table(df_t$Group)
sum(table(df_t$Group))

mapping <- df_t[c('N_ambiguous','N_noFeature','N_multimapping','N_unmapped')]
names(mapping) <- c('Ambiguous','No feature','Multimapping','Unmapped')

sapply(mapping, quantile)

require(tidyverse)

# Convert wide to long format
df_long <- mapping %>%
  rownames_to_column(var = "Gene") %>%  # Move row names to a column
  pivot_longer(cols = -Gene,            # Select all columns except 'Gene'
               names_to = "Sample",     # New column for former column names
               values_to = "Value")     # New column for former cell values

# Print the long format data frame
print(df_long)

# Create boxplots
ggplot(df_long, aes(x = Sample, y = Value, fill=Sample)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill=NA),
        axis.text  = element_text(size = 12, color='black'),
        axis.text.x = element_text(angle = -45, hjust = -0.1),
        axis.title  = element_text(size = 14),
        legend.position = "none"
  )+  labs(x = "",y = "Number of reads")

df_t$N_ambiguous<-NULL
df_t$N_noFeature<-NULL
df_t$N_multimapping<-NULL
df_t$N_unmapped<-NULL
group_df <- df_t["Group"]
group_df$Group2<-group_df$Group
group_df$Group2[group_df$Group%in%c('V02.0','V05.0')]<-"V02.0 and V05.0"
group_df$Group2[group_df$Group%in%c('V02a','V05a')]<-"V02a and V05a"

require(dplyr)
group_df <- group_df %>%
  mutate(Group2 = case_when(
    Group %in% c("V02.0", "V05.0") ~ "V02 and V05",
    Group %in% c("V02a", "V05a") ~ "V02a and V05a"
    ))

df_t$Group<-NULL

n_top_genes = 50
total_gene_counts <- data.frame(reads=colSums(df_t))
total_gene_counts$gene_read_counts <- total_gene_counts$reads
total_gene_counts <- total_gene_counts[order(total_gene_counts$reads, decreasing = TRUE),]
top_genes <- total_gene_counts[1:n_top_genes, ]

top_n_genes_df <- df_t[colnames(df_t)%in%rownames(top_genes)]

pheatmap(log2(t(top_n_genes_df)),
                  treeheight_row = 0,
                 treeheight_col = 0,
         annotation_col = group_df[c("Group")],
         show_colnames = FALSE)

############ DEA #######

require(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = t(df_t),
                              colData = group_df,
                              design = ~ Group2) # 
dds

dds <- DESeq(dds)

res <- results(dds)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Print summary of results
summary(res)

plotMA(res)

res_df <- as.data.frame(res)
res_df$log2FoldChange <- as.numeric(res_df$log2FoldChange)
res_df$padj <- as.numeric(res_df$padj)
res_df <- res_df[!is.na(res_df$padj),]
#res_df <- res_df[res_df$padj<=0.05,]
#res_df <- res_df[abs(res_df$log2FoldChange)>1.5,]

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed")


plotDispEsts(dds, ylim = c(1e-6, 1e2) )

# Convert the DESeq2 results to a data frame
res_df <- as.data.frame(res)

# Remove rows with NA values in the log2FoldChange column
res_df <- res_df[!is.na(res_df$log2FoldChange), ]

# Select the top 10 genes with the highest positive log2 fold change
top10_pos <- res_df %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

# Select the top 10 genes with the highest negative log2 fold change
top10_neg <- res_df %>%
  arrange(log2FoldChange) %>%
  head(10)

# Combine the top 10 positive and top 10 negative genes
top_genes <- bind_rows(top10_pos, top10_neg)

# Order the combined data by log2 fold change
top_genes <- top_genes %>%
  arrange(log2FoldChange)

# Create a column for the rank of each gene based on log2 fold change
top_genes$rank <- seq_len(nrow(top_genes))
top_genes$gene <- rownames(top_genes)

# Determine vertical positioning of labels based on log2 fold change
top_genes$label_y <- ifelse(top_genes$log2FoldChange > 0, top_genes$log2FoldChange + 0.2, top_genes$log2FoldChange - 0.2)

# Create the waterfall plot using ggplot2
ggplot(top_genes, aes(x = rank, y = log2FoldChange, fill = log2FoldChange > 0, label = gene)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1", direction = -1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill=NA),
        axis.text  = element_text(size = 12),
        legend.position = "none"
  )+
  labs(x = "",
       y = "Log2 Fold Change") +
  geom_text(data = top_genes,
            aes(y = 0, label = gene, angle = ifelse(log2FoldChange > 0, -90, 90)),
            hjust = -0.1,
            size = 3.5) 



# Regularized log transformation for visualization
rld <- vst(dds, blind = FALSE)

# Select top 25 genes by smallest adjusted p-value
topgenes <- head(rownames(resOrdered), 50)

# Extract normalized counts for these genes
mat <- assay(rld)[topgenes, ]

# Plot heatmap
pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = TRUE, annotation_col = group_df['Group2'],
         treeheight_row = 0,
         treeheight_col=0,
         show_colnames = 0)



vsd_matrix <- assay(rld)
# Calculate Bray-Curtis distance
bc_dist <- vegan::vegdist(t(vsd_matrix), method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bc_dist, eig = TRUE)

# Create a data frame with PCoA results
pcoa_df <- data.frame(Sample = rownames(pcoa_result$points),
                      PCoA1 = pcoa_result$points[, 1],
                      PCoA2 = pcoa_result$points[, 2])

# Merge with metadata
group_df$Sample<-rownames(group_df)
pcoa_df <- merge(pcoa_df, group_df, by = "Sample")

# Plot PCoA
p2 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group2)) +
  geom_point(size = 3) +
  theme(
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill=NA),
        axis.text  = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right"
  )+  labs(title = "",
       x = "PCoA1",
       y = "PCoA2")+
  stat_ellipse(aes(group = Group2), level = 0.95) +
  scale_color_brewer(palette = "Set1", direction = -1)
p2

p3<-ggplot(pcoa_df[pcoa_df$Group!="Others",], aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 3) +
  theme(
        panel.background = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(fill=NA),
        axis.text  = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right")+  labs(title = "",
           x = "PCoA1",
           y = "PCoA2")+
  stat_ellipse(aes(group = Group), level = 0.95) +
  scale_color_brewer(palette = "Set2", direction = -1)
  
p3
cowplot::theme_cowplot(p3,p2)
gridExtra::grid.arrange(p3,p2, nrow=1)

########## Gene set analysis #########################

require(fgsea)

# Example: Extract gene symbols and their corresponding statistics
gene_symbols <- rownames(res)
logFC <- res$log2FoldChange
pvalues <- res$padj  # or res$pvalue depending on what you have

# Create a data frame with gene symbols, log2FoldChange, and p-values
gene_stats <- data.frame(symbol = gene_symbols, logFC = logFC, pvalue = pvalues)

# Filter genes based on a significance threshold (e.g., adjusted p-value < 0.05)
gene_stats <- gene_stats[gene_stats$pvalue < 0.05, ]

gene_stats$GeneSymbol <- gsub("\\.\\d+$", "", gene_stats$symbol)

gene_stats <- na.omit(gene_stats)
gene_stats2 <- gene_stats
# Rank genes based on log2FoldChange (negative values should be ranked lower)
gene_rank <- order(gene_stats$logFC, decreasing = TRUE)

# Prepare gene list and associated statistics for GSEA
gene_list <- gene_stats$symbol[gene_rank]
gene_stats <- gene_stats$logFC[gene_rank]  # Assuming log2FoldChange as the statistic


# Ensure gene_stats is named with gene symbols
names(gene_stats) <- gene_list
gene_ids <- gsub("\\.\\d+$", "", names(gene_stats))
names(gene_stats) <- gene_ids

# Load necessary libraries
pathways_reactome <- msigdbr::msigdbr(species = "Homo sapiens")
#pathways_reactome <- msigdbr::msigdbr(species = "Homo sapiens")

gene_stats2$Entez <- pathways_reactome$human_gene_symbol[match(gene_stats2$GeneSymbol, pathways_reactome$human_ensembl_gene)]
gene_stats2 <- gene_stats2[abs(gene_stats2$logFC)>=1.5,]

write.csv(na.omit(gene_stats2), 'tmpdf.csv', quote = FALSE)

# Perform GSEA using fgsea
gsea_result <- fgsea(pathways = pathways_reactome,
                     stats = gene_stats)

# View the top enriched pathways
head(gsea_result)



