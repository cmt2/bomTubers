library(DESeq2)
library(tximport)
library(tidyverse)
library(ggthemes)
setwd("~/Desktop/trans_diff")
#### load in data and process stuff ####

#conditions <- read.table("samples_described.txt",sep="\t", header = TRUE,
#                         row.names = 1)
#files <- file.path(paste0("isoform_results/RSEM.isoforms.results.", 
#                   rownames(conditions)))
#names(files) <- rownames(conditions)
#
#txi.rsem <- tximport(files, txIn = TRUE, txOut = TRUE,  
#                     geneIdCol =  "transcript_id", 
#                     txIdCol = "transcript_id",
#                     lengthCol = "effective_length",
#                     abundanceCol = "FPKM",
#                     countsCol = "expected_count")
#
#dds <- DESeqDataSetFromTximport(txi.rsem,
#                                colData = conditions,
#                                design = ~ condition)
#save(dds, file = "dds.RData")

load("dds.RData")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$condition, levels = c("rhi", "roo", "SAM", "tub"))

##### Make a beautiful PCA #####

vst <- vst(dds, blind = FALSE)
pcaData_vst <- plotPCA(vst, returnData = TRUE)
percentVar <- round(100 * attr(pcaData_vst, "percentVar"))
myColors <- c("tub" = "#a6cee3", "roo" = "#1f78b4",
              "rhi" = "#b2df8a", "SAM" = "#33a02c")
samples_to_individuals <- data.frame(
  sample = c("tub_2","roo_2","rhi_2",
             "SAM_3","SAM_1","SAM_2",
             "tub_1","roo_3","rhi_1",
             "rhi_3","roo_1","tub_3"),
  individual = c("A", "A", "A", "A",
                 "B", "C", "D", "D",
                 "D", "E", "E", "E"))
pcaData_vst <- inner_join(pcaData_vst, 
                          samples_to_individuals, 
                          by = c("name"="sample"))
pdf("pca_2.pdf", width = 5, height = 5, useDingbats = FALSE)
ggplot(pcaData_vst, aes(x = PC1, 
                        y = PC2, 
                        color = group)) +
  geom_point(size = 2) +
  geom_text(aes(label = individual),size = 3,
            color = "black",
            #check_overlap = TRUE,
            nudge_x = 4, nudge_y = 0,
            show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_equal() +
  scale_color_manual(values = myColors,
                     name = "Tissue Type",
                     breaks = c("SAM","rhi",
                                "roo","tub"),
                     labels = c("Aerial Shoot",
                                "Rhizome",
                                "Fibrous Root",
                                "Tuberous Root") ) +
  ggtitle("VST transformed transcript counts") +
  theme_few() +
  guides(color = guide_legend(override.aes = 
                                list(shape = 20))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0, 0),
        legend.justification = c(0,0),
        legend.key.size = unit(.2, "cm"),
        legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        legend.text.align = 0,
        legend.title.align = 0,
        #legend.box.background = element_rect(colour = "black"),
        #panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

dev.off()
#####

##### Run Differential Gene Expression #####
dds <- DESeq(dds)

results(dds, contrast = c("condition","tub","roo"),
         alpha  = 0.01, lfcThreshold = 2) %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) %>%
  filter(log2FoldChange > 2 | log2FoldChange < -2) %>%
  filter(padj < 0.01) -> res
#####

### read in annotated transcriptome and subset for only selected genes ###
annotated <- read.csv(file = "/Users/carrietribble/Documents/Berkeley_IB/research_projects/transcriptomes/trans_annotate/trinotate_annotation_report.csv",
                      na.strings = ".", stringsAsFactors = FALSE)

inner_join(x = annotated, y = res, by = "transcript_id") %>%
  distinct(transcript_id, gene_ontology_blast, .keep_all = T) -> res_ann
write.csv(res_ann, "res.csv")

res_ann %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(abs_log2FoldChange) %>%
  tail(n = 10) -> res_ann_strict
write.csv(res_ann_strict, "res_strict.csv")

### Plot expression of those identified genes 
pdf("logfoldchange.pdf")
results(dds, contrast = c("condition", "tub","roo"),
        alpha  = 0.01, lfcThreshold = 2) %>%
  as.data.frame() %>%
  mutate(color = ifelse(.$padj < 0.01 &
                          abs(.$log2FoldChange) > 2, "darkred","black")) -> res

plot(res$log2FoldChange, pch = 20, col = alpha("black",0.3),
     ylab = "Log2 Fold Change", xlab = NA,
     xaxt = 'n',
     main = "Logfold Change in Transcript Counts\n Between Fibrous and Tuberous Roots")
points(x =  which(res$color == "darkred"),
       y =  res$log2FoldChange[which(res$color == "darkred")], 
       col  = "darkred", pch = 20, cex =  1.2)
dev.off()




