#### search for 'candidate' genes #####
library(tidyverse)
library(DESeq2)
library(knitr)
library(kableExtra)
# read in annotated transcriptome
annotated <- read.csv(file = "~/trans_diff/data/trinotate_annotation_report.csv",
                      na.strings = ".", stringsAsFactors = FALSE)

# read in and analyze diff expression data
load("~/trans_diff/data/dds.RData")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition,levels = c("rhi", "roo", "SAM", "tub"))
dds <- DESeq(dds)
results(dds, contrast = c("condition","tub","roo"),
        alpha  = 0.01, lfcThreshold = 2) %>%
  as.data.frame() %>%
  mutate(transcript_id = rownames(.) ) -> diff_expression


##### groups/ cellular processes #####
# Starch biosynthesis pathway ("GO:0019252") [logical, also found upregulated in most studies]
# 168
starch <- annotated[grep("GO:0019252",annotated$gene_ontology_blast),]
# Hormones cytokinin, abscisic acid, and auxin [Noh et al 2010]
# cty = 356; ab = 1357; auxin = 716
cytokinin <- annotated[grep("GO:0009735",annotated$gene_ontology_blast),]
abscisic_acid <- annotated[grep("GO:0009737",annotated$gene_ontology_blast),]
auxin <- annotated[grep("GO:0009733",annotated$gene_ontology_blast),]
# MADS Box: IbMADS1, IbMADS3, IbMADS4, IbMADS10, IbMADS79, IbAGL17, IbAGL20 and SRD1 [Noh et al 2010]
# 116
MADS <- annotated[grep("MADS-box",annotated$sprot_Top_BLASTX_hit),] 
# KNOX genes [Noh et al 2010]
# 18
knox <- annotated[grep("KNOX",annotated$Pfam),]

### Root thickening 
# (lack of) gibberellin https://academic.oup.com/aob/article/110/2/373/2769213
# 295
gibb <- annotated[grep("GO:0009739", annotated$gene_ontology_blast),]
# * Expansins 
# 117
expansin <- annotated[grep("Expansin",annotated$sprot_Top_BLASTX_hit),]
# * lignin biosynthesis (downregulated)
# 325
lignin <- annotated[grep("GO:0009809", annotated$gene_ontology_blast),]
# * 14-3-3 genes
# * 14-3-3 genes
# 62
fourteen <- annotated[grep("14-3-3-like protein", annotated$sprot_Top_BLASTX_hit),]
# * calcium‐dependent protein kinase [manihot]
# 180
cdpk <- annotated[grep("Calcium-dependent protein kinase", annotated$sprot_Top_BLASTX_hit),]

# lipoxygenase
lipo <- annotated[grep("lipoxygenase", annotated$gene_ontology_blast),]



### identify diff exp levels for these processes
processes <- list(starch,cytokinin,abscisic_acid,
                  auxin,MADS,knox,gibb,expansin,
                  lignin,fourteen,cdpk, lipo)

processes_res <- list()
for (i in 1:length(processes)){
  processes_res[[i]] <- inner_join(x = diff_expression, y = processes[[i]], by = "transcript_id")
}
names(processes_res) <- c("starch","cytokinin","abscisic_acid","auxin",
                          "MADS","knox","gibb","expansin","lignin","fourteen","cdpk", "lipo")
log2FoldChange <- NA
pvalue <- NA
process <- NA
transcript_id <- NA
for (i in 1:length(processes_res)){
  log2FoldChange <- append(log2FoldChange, processes_res[[i]]$log2FoldChange)
  pvalue <- append(pvalue, processes_res[[i]]$padj)
  transcript_id <- append(transcript_id, processes_res[[i]]$transcript_id)
  process <- append(process, rep(names(processes_res)[i], 
                                     times = length(processes_res[[i]]$log2FoldChange)))
}

processes_all <- data.frame(log2FoldChange = append(log2FoldChange, 
                                                    diff_expression$log2FoldChange),
                            pvalue = append(pvalue, 
                                            diff_expression$pvalue),
                            transcript_id = append(transcript_id,
                                                   diff_expression$transcript_id),
                            process = c(process,
                                        rep("all", times = nrow(diff_expression))),
                            sig = replace_na(append(pvalue,diff_expression$pvalue) < 0.01, FALSE))[-1,]


# test for differences in mann whitney U test statistc 
pdf("~/Desktop/test_ws_two_tailed.pdf", 
    width = 15, height = 7)
par(mfrow = c(1,2))
for (i in 1:length(processes_res)) {
  #means <- NA
  #sds <- NA
  w_stat <- NA
  print(names(processes_res)[i])
  
  for (j in 1:10000) {
  s <- sample_n(diff_expression, nrow(processes_res[[i]]))
  s$abs_log2foldchange <- abs(s$log2FoldChange)
  #means[j] <- mean(s$abs_log2foldchange)
  #sds[j] <- sd(s$abs_log2foldchange)
  w_stat[j] <- wilcox.test(s$abs_log2foldchange,
                           abs(diff_expression$log2FoldChange))$statistic
  }
  w_stat <- w_stat / (nrow(processes_res[[i]]) * length(diff_expression$log2FoldChange))
  w <- wilcox.test(abs(processes_res[[i]]$log2FoldChange),
                  abs(diff_expression$log2FoldChange))
  #tot_mean_values <- c(means, mean(abs(processes_res[[i]]$log2FoldChange)))
  #xlims_means <- c((min(tot_mean_values)-(sum(abs(range(tot_mean_values)))/4)), 
  #                 (max(tot_mean_values)+(sum(abs(range(tot_mean_values)))/4)))
  #hist(means,
  #     main=paste0(names(processes_res)[i], ": P-Value = ", round(w$p.value, digits = 4)), 
  #     xlab= "Mean Log2FoldChange",
  #     col="grey", 
  #     border = NA,
  #     las=1, 
  #     xlim = xlims_means,
  #     breaks = 50, 
  #     prob = TRUE)
  #abline(v = mean(abs(processes_res[[i]]$log2FoldChange)),
  #       col = "red")
  #     
  #tot_sd_values <- c(sds, sd(abs(processes_res[[i]]$log2FoldChange)))
  #xlims_sds <- c((min(tot_sd_values)-(sum(abs(range(tot_sd_values)))/4)), 
  #                 (max(tot_sd_values)+(sum(abs(range(tot_sd_values)))/4)))
  #
  #hist(sds, 
  #     main=paste0(names(processes_res)[i], ": P-Value = ", round(w$p.value, digits = 4)), 
  #     xlab= "SD Log2FoldChange", 
  #     col="grey", 
  #     border = NA,
  #     las=1, 
  #     xlim = xlims_sds,
  #     breaks = 50, 
  #     prob = TRUE)
  #abline(v = sd(abs(processes_res[[i]]$log2FoldChange)),
  #       col = "red")
  #
  
  #tot_w_values <- c(w_stat, 
  #                  w$statistic/ (nrow(processes_res[[i]]) * 
  #                                length(diff_expression$log2FoldChange)))
  #xlims_ws <- c((min(tot_w_values)-(sum(abs(range(tot_w_values)))/4)), 
  #                (max(tot_w_values)+(sum(abs(range(tot_w_values)))/4)))
  my_hist <- hist(w_stat, breaks = 50, plot = F)
  my_color <-  ifelse(my_hist$breaks < quantile(w_stat, 0.025) | my_hist$breaks > quantile(w_stat, 0.975), "grey", "dodgerblue" )
  
  hist(w_stat, 
       main=paste0(names(processes_res)[i], ": P-Value = ", round(w$p.value, digits = 4)), 
       xlab= "W Statistic Effect Size", 
       col= my_color, 
       border = NA,
       las=1, 
       xlim = c(0,1),
       breaks = 50, 
       prob = TRUE)
  abline(v = w$statistic/ (nrow(processes_res[[i]]) * length(diff_expression$log2FoldChange)),
         col = "red")
}
dev.off()

# test for differences in number of significant genes 
processes_all %>%
  group_split(process) -> processes_groups

pdf("~/Desktop/test_sig_genes_col.pdf", 
    width = 15, height = 7)
par(mfrow = c(1,2))
for (i in 1:length(processes_groups)) {
  if (processes_groups[[i]]$process[1] != "all") {
    df <- processes_groups[[i]]
    r <- nrow(df) 
    n_sigs <- NA
    sig_in_group <- sum(df$sig)
    for (j in 1:10000) {
      s <- sample_n(processes_groups[[1]], r)
      n_sigs[j] <- sum(s$sig)
    }
    
    x_values <- c(max(n_sigs), sig_in_group)
    breaks <- seq(from = -0.5, to = max(x_values + 0.5))
    
    my_hist <- hist(n_sigs, breaks = breaks, plot = F)
    
    table <- sort(table(n_sigs)/length(n_sigs), decreasing = TRUE)
    cs_table <- as.table(table[1:min(which((cumsum(table) >= 0.95) == TRUE))])
    credible_set <- as.integer(names(cs_table))
    
    my_color <- my_hist$breaks
    my_color[my_color %in% c(credible_set - 0.5)] <- "dodgerblue"
    my_color[my_color != "dodgerblue"] <- "grey"

    hist(n_sigs, 
         breaks = breaks,
         xaxt='n',
         col = my_color,
         xlab = "Number of Significant Isoforms",
         main = paste0(unique(df$process)))
    text(x = sig_in_group, 
         y = my_hist$counts[sig_in_group + 1], 
         labels = "*", 
         col = "red",
         cex = 4)
    axis(side = 1, 
         at = (0:max(x_values)), 
         labels = (0:max(x_values)),
         tick = F,
         )
  }
}
dev.off()


# reorder processes by significance
processes_all$process <- factor(processes_all$process, 
                                levels = c("all", "expansin", "lignin", "MADS",
                                           "starch", "fourteen", "lipo", "abscisic_acid", 
                                           "auxin", "cdpk", "cytokinin", "gibb", "knox"))

#reorder within processes to plot non-sig first 
processes_all <- processes_all[order(processes_all$process, processes_all$sig),]

#### plot results ####

png("~/trans_diff/figures/bomTubCandidates.png", height = 5, width = 10, units = "in", res = 500)
processes_all %>%
  ggplot(aes(x=process, y = log2FoldChange)) + 
  geom_jitter(position=position_jitter(0.4),
              aes(color = sig, shape = sig, size = process)) +
  geom_hline(yintercept = median(processes_all[processes_all$process == "all", "log2FoldChange"]), 
             color = "black", linetype = "dashed") +
  geom_violin(fill = "beige", draw_quantiles = 0.5, trim = T) +
  scale_color_manual(values = c("TRUE" = "red", 
                                "FALSE" = "#bbbbbb60"),
                     labels = c("Not Signficant","Significant")) +
  scale_shape_manual(values = c("TRUE" = 3,
                                "FALSE" = 20), 
                     labels = c("Not Signficant","Significant")) +
  scale_size_manual(values = c("all" = 0.5,
                               "expansin" = 2.5, 
                               "lignin" = 2.5,
                               "MADS" = 2.5,
                               "starch" = 2.5,
                               "fourteen" = 2.5,
                               "lipo" = 2.5,
                               "abscisic_acid" = 2.5,
                               "auxin" = 2.5,
                               "cdpk" = 2.5,
                               "cytokinin" = 2.5,
                               "gibb" = 2.5,
                               "knox" = 2.5)) +
  scale_x_discrete(name = NULL, 
                   labels = c("all" = "All Genes", 
                              "expansin" = "Expansins", 
                              "lignin" = "Lignin",
                              "MADS" = "MADS Box",
                              "starch" = expression(paste("Starch"^"*")),
                              "fourteen" = "14-3-3",
                              "lipo" = "Lipoxygenases",
                              "abscisic_acid" = "Abscisic Acid",
                              "auxin" = "Auxin",
                              "cdpk" = "CDPK",
                              "cytokinin" = "Cytokinin",
                              "gibb" = "Gibberellin",
                              "knox" = "Knox")) +
  scale_y_continuous(name = "Absolute Value of Log2 Fold Change in Expression") +
  theme_few() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 12,
                                   color = c("black", 
                                             rep("red", 4),
                                             rep("dodgerblue", 2),
                                             rep("black", 6))),
        legend.title = element_blank(), legend.position = c(1,-25)) 
dev.off()


# make table of groups sign. genes 
processes_all[which(processes_all$sig & processes_all$process != "all"), 
              c("log2FoldChange","transcript_id","process")] %>%
  inner_join(annotated) -> process_ann

transcript <- character()
ann_name <- character()
names <- character()
process <- character()
log2FoldChange <- numeric()
cell_components <- character()
mol_functions <- character()
bio_process <- character()

for (i in 1:length(process_ann$gene_ontology_blast)){
  transcript[i] <- process_ann$transcript_id[i]
  ann_name[i] <- strsplit(strsplit(strsplit(strsplit(process_ann$sprot_Top_BLASTP_hit[i], "\\^")[[1]][6], "=")[[1]][2], ";")[[1]][1], "\\{")[[1]][1]
  names[i] <- paste0(transcript[i], ": ", ann_name[i])
  log2FoldChange[i] <- round(process_ann$log2FoldChange[i], digits = 2)
  process[i] <- as.character(process_ann$process[i])
  c <- character()
  m <- character()
  b <- character() 
  if (!is.na(process_ann$gene_ontology_blast[i])){
    l <- unlist(strsplit(process_ann$gene_ontology_blast[i],  "`"))
    for (j in 1:length(l)){
      s <- unlist(strsplit(l[j], "\\^"))
      if (s[2] == "cellular_component"){
        c <- append(c, s[3])
      } else if (s[2]  == "molecular_function"){
        m <- append(m, s[3])
      } else if (s[2] == "biological_process"){
        b <- append(b, s[3])
      }
    }
    cell_components[i] <- paste(c, collapse = ", ")
    mol_functions[i] <- paste(m, collapse = ", ")
    bio_process[i] <- paste(b, collapse = ", ")
  } else {
    cell_components[i] <- NA
    mol_functions[i] <- NA
    bio_process[i] <- NA
  }
}

data.frame(process = recode(process, cdpk = "Calcium-Dependent Protein Kinases",
                            lignin = "Lignin Biosynthesis",
                            gibb = "Response to Gibberellin",
                            MADS = "MADS-Box Genes",
                            abscisic_acid = "Abscisic Acid Signaling",
                            cytokinin = "Cytokinin Signaling",
                            starch = "Starch Biosynthesis"),
           names = names,
           log2FoldChange = log2FoldChange,
           cell_components = cell_components, 
           mol_functions = mol_functions,
           bio_process = bio_process) %>%
  arrange(process,-log2FoldChange) %>%
  kable("latex", longtable = T, booktabs = T, align=c(rep('l',2),'c',rep('l',3)),
        col.names = c("Process Group","Transcript ID: Annotated Name (SPROT)",
                      "Log2 Fold Change",
                      "Cellular Components", "Molecular Functions", 
                      "Biological Processes")) %>%
  add_header_above(c("","","", "Gene Ontology" = 3),
                   bold = T) %>%
  column_spec(1, width = "10em", bold = T) %>%
  kable_styling(latex_options = c("striped") ) %>%
  row_spec(0, bold=TRUE) %>%
  column_spec(2, width = "15em") %>%
  column_spec(3, width = "5em") %>%
  column_spec(4, width = "30em") %>%
  column_spec(5, width = "30em") %>%
  column_spec(6, width = "30em") 


  
##### non-group candidates #####
####

# rice root thickening
# qRT9 /  OsbHLH120 https://academic.oup.com/jxb/article/66/9/2723/679410
# blasted sequence to transcriptome - 3 hits, 2 found in diff expr df, none sign.
#diff_expression[which(diff_expression$transcript_id == "TRINITY_DN122858_c1_g1_i1"),]
#diff_expression[which(diff_expression$transcript_id == "TRINITY_DN111446_c0_g1_i1"),]

rice_AA_blast <- c("TRINITY_DN122858_c1_g1_i1","TRINITY_DN111446_c0_g1_i1","TRINITY_DN109913_c0_g1_i1","TRINITY_DN122858_c0_g2_i1","TRINITY_DN121077_c10_g1_i3","TRINITY_DN121077_c10_g1_i1","TRINITY_DN122858_c1_g1_i2","TRINITY_DN112828_c0_g1_i1","TRINITY_DN112828_c0_g2_i1","TRINITY_DN113155_c0_g1_i1","TRINITY_DN116814_c5_g2_i1","TRINITY_DN115984_c0_g1_i2","TRINITY_DN115984_c0_g1_i1","TRINITY_DN114144_c0_g1_i2","TRINITY_DN118606_c2_g2_i9","TRINITY_DN118606_c2_g2_i1","TRINITY_DN118606_c2_g2_i7","TRINITY_DN121793_c5_g2_i6","TRINITY_DN121793_c5_g2_i4","TRINITY_DN121793_c5_g2_i5","TRINITY_DN121793_c5_g2_i7","TRINITY_DN114964_c0_g1_i2","TRINITY_DN117607_c0_g1_i1","TRINITY_DN114964_c0_g1_i1","TRINITY_DN129122_c1_g1_i2","TRINITY_DN129122_c1_g1_i1")
rice_AA <- diff_expression[diff_expression$transcript_id %in% rice_AA_blast,]

# Ipomea: 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2837253/pdf/erp399.pdf
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-460
# https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-018-1307-4
# * IDD5 Protein indeterminate-domain 5, regulates starch synthesis 
# find 3 genes, 2 in DE, none significant
#idd5 <- annotated[grep("IDD5",annotated$sprot_Top_BLASTX_hit),] 
#idd5_de <- inner_join(x = diff_expression, y = idd5, by = "transcript_id")
#finds 63, none sig
idd5_AA_blast <- c("TRINITY_DN122462_c4_g3_i8","TRINITY_DN122462_c4_g3_i3","TRINITY_DN124922_c5_g3_i1","TRINITY_DN119660_c4_g1_i4","TRINITY_DN124922_c5_g3_i7","TRINITY_DN124922_c2_g1_i1","TRINITY_DN122462_c4_g3_i2","TRINITY_DN122462_c4_g3_i7","TRINITY_DN124922_c5_g3_i3","TRINITY_DN122462_c4_g3_i5","TRINITY_DN124922_c5_g3_i2","TRINITY_DN124922_c5_g3_i6","TRINITY_DN119867_c3_g1_i1","TRINITY_DN119660_c4_g2_i1","TRINITY_DN122462_c4_g3_i10","TRINITY_DN122462_c4_g3_i4","TRINITY_DN119660_c4_g1_i1","TRINITY_DN124922_c4_g1_i2","TRINITY_DN119867_c3_g1_i2","TRINITY_DN124922_c1_g1_i1","TRINITY_DN124922_c5_g2_i4","TRINITY_DN124922_c5_g1_i1","TRINITY_DN126974_c0_g6_i1","TRINITY_DN124922_c5_g1_i2","TRINITY_DN126974_c0_g1_i1","TRINITY_DN124922_c5_g2_i3","TRINITY_DN119867_c1_g1_i1","TRINITY_DN124922_c5_g2_i2","TRINITY_DN126974_c0_g6_i2","TRINITY_DN124922_c5_g2_i5","TRINITY_DN126974_c0_g2_i4","TRINITY_DN126974_c0_g5_i1","TRINITY_DN126974_c0_g1_i2","TRINITY_DN126974_c0_g6_i4","TRINITY_DN124922_c5_g2_i1","TRINITY_DN126974_c0_g2_i1","TRINITY_DN126974_c0_g2_i2","TRINITY_DN126974_c0_g6_i3","TRINITY_DN126974_c0_g5_i2","TRINITY_DN122462_c4_g4_i2","TRINITY_DN122462_c4_g8_i1","TRINITY_DN126974_c0_g2_i3","TRINITY_DN124922_c5_g6_i1","TRINITY_DN119867_c2_g1_i1","TRINITY_DN124922_c5_g4_i1","TRINITY_DN126399_c3_g3_i1","TRINITY_DN122462_c4_g5_i1","TRINITY_DN124922_c1_g2_i1","TRINITY_DN125234_c0_g1_i3","TRINITY_DN125234_c0_g1_i2","TRINITY_DN124629_c1_g1_i1","TRINITY_DN119550_c2_g6_i2","TRINITY_DN106353_c0_g1_i1","TRINITY_DN119550_c2_g1_i2","TRINITY_DN119550_c2_g2_i2","TRINITY_DN124922_c5_g7_i1","TRINITY_DN124922_c5_g7_i2","TRINITY_DN119550_c2_g3_i1","TRINITY_DN119550_c2_g1_i1","TRINITY_DN119550_c2_g3_i2","TRINITY_DN119867_c0_g1_i1","TRINITY_DN119550_c2_g2_i1","TRINITY_DN122462_c4_g9_i1","TRINITY_DN124922_c5_g5_i1","TRINITY_DN150478_c0_g1_i1","TRINITY_DN124922_c5_g5_i3","TRINITY_DN124922_c4_g1_i1","TRINITY_DN124922_c5_g3_i5","TRINITY_DN116737_c0_g1_i2","TRINITY_DN116737_c0_g1_i1")
idd5_AA <- diff_expression[diff_expression$transcript_id %in% idd5_AA_blast,]

# * WOX4 
# one result, not significant 
#WOX4 <- annotated[grep("WOX4",annotated$sprot_Top_BLASTX_hit),]
#WOX4_de <- inner_join(x = diff_expression, y = WOX4, by = "transcript_id")
#finds 31, none sig
wox4_AA_blast <- c("TRINITY_DN120806_c0_g3_i2","TRINITY_DN120806_c0_g1_i3","TRINITY_DN120806_c0_g1_i5","TRINITY_DN120806_c0_g5_i1","TRINITY_DN117102_c3_g2_i1","TRINITY_DN122482_c0_g1_i1","TRINITY_DN123697_c2_g2_i9","TRINITY_DN123697_c2_g2_i3","TRINITY_DN123697_c2_g2_i7","TRINITY_DN119023_c2_g1_i11","TRINITY_DN119023_c2_g1_i6","TRINITY_DN119023_c2_g1_i2","TRINITY_DN119023_c2_g1_i4","TRINITY_DN119023_c2_g1_i3","TRINITY_DN123697_c2_g2_i10","TRINITY_DN119023_c2_g1_i1","TRINITY_DN119023_c2_g1_i7","TRINITY_DN123697_c2_g2_i2","TRINITY_DN123697_c2_g2_i8","TRINITY_DN120214_c1_g1_i5","TRINITY_DN120214_c1_g1_i3","TRINITY_DN120214_c1_g1_i7","TRINITY_DN120214_c1_g1_i2","TRINITY_DN120214_c1_g1_i9","TRINITY_DN120214_c1_g1_i6","TRINITY_DN122355_c7_g1_i1","TRINITY_DN122355_c7_g1_i2","TRINITY_DN120214_c1_g1_i4","TRINITY_DN123697_c2_g2_i1","TRINITY_DN123697_c2_g2_i11","TRINITY_DN123017_c0_g1_i1")
wox4_AA <- diff_expression[diff_expression$transcript_id %in% wox4_AA_blast,]

# Manihot esculenta
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1399-3054.2010.01389.x
# * Sulfite reductase
# find 17 genes, non significant
#sulf <- annotated[grep("Sulfite reductase", annotated$sprot_Top_BLASTX_hit),]
#sulf_de <- inner_join(x = diff_expression, y = sulf, by = "transcript_id")
# fine 9 genes, none sig
sulf_AA_blast <- c("TRINITY_DN126636_c0_g2_i2","TRINITY_DN126636_c0_g2_i1","TRINITY_DN127936_c6_g1_i1","TRINITY_DN126636_c0_g1_i2","TRINITY_DN126636_c0_g1_i1","TRINITY_DN127936_c6_g1_i2","TRINITY_DN123380_c1_g2_i2","TRINITY_DN123380_c1_g2_i3","TRINITY_DN123380_c1_g2_i1","TRINITY_DN6600_c0_g1_i1 l")
sulf_AA <- diff_expression[diff_expression$transcript_id %in% sulf_AA_blast,]


# * ft genes
# Potato 
#StSP6A [FT-like gene]*
#StFDL1 [FD-like gene]
#Onion
#AcFT1 [FT-like gene] *
#AcFT4 [FT-like gene]
# first, search for genes annotated as FT
# search produces 4, only 3 in DE dataset, none significant
ft_search_ann <- annotated[grep("FLOWERING LOCUS T", annotated$sprot_Top_BLASTX_hit), ]
ft_search_ann_de <- inner_join(x = diff_expression, y = ft_search_ann, by = "transcript_id")

# blast StSP6A and pick results with E value greater that 0.010
# find 27 genes, 22 in DE, 1 significant 
pot_blast <- c("TRINITY_DN120193_c1_g1_i3",
               "TRINITY_DN120193_c1_g1_i1",
               "TRINITY_DN120856_c1_g2_i2",
               "TRINITY_DN124927_c0_g2_i1",
               "TRINITY_DN124927_c0_g1_i1",
               "TRINITY_DN124927_c0_g1_i7",
               "TRINITY_DN124927_c0_g1_i6",
               "TRINITY_DN124927_c0_g1_i2",
               "TRINITY_DN124927_c0_g1_i9",
               "TRINITY_DN124927_c0_g1_i4",
               "TRINITY_DN127729_c5_g3_i3",
               "TRINITY_DN111334_c0_g1_i1",
               "TRINITY_DN129076_c1_g1_i1",
               "TRINITY_DN45900_c0_g1_i1",
               "TRINITY_DN127729_c5_g3_i2",
               "TRINITY_DN124927_c0_g1_i8",
               "TRINITY_DN120193_c1_g1_i2",
               "TRINITY_DN124927_c0_g1_i3",
               "TRINITY_DN120856_c1_g2_i1",
               "TRINITY_DN120856_c1_g1_i3",
               "TRINITY_DN120856_c1_g1_i11",
               "TRINITY_DN127729_c5_g3_i1",
               "TRINITY_DN127729_c5_g9_i1",
               "TRINITY_DN117914_c5_g2_i1",
               "TRINITY_DN120193_c1_g2_i1",
               "TRINITY_DN10791_c0_g1_i1",
               "TRINITY_DN56951_c0_g1_i1")
ft_blast_potato <- diff_expression[diff_expression$transcript_id %in% pot_blast,]

ft1_blast_onion_names <- c("TRINITY_DN124927_c0_g2_i1",
                 "TRINITY_DN120193_c1_g1_i1",
                 "TRINITY_DN120193_c1_g1_i3",
                 "TRINITY_DN120856_c1_g2_i2",
                 "TRINITY_DN124927_c0_g1_i1",
                 "TRINITY_DN124927_c0_g1_i7",
                 "TRINITY_DN124927_c0_g1_i4",
                 "TRINITY_DN124927_c0_g1_i6",
                 "TRINITY_DN124927_c0_g1_i2",
                 "TRINITY_DN124927_c0_g1_i9",
                 "TRINITY_DN127729_c5_g3_i3",
                 "TRINITY_DN111334_c0_g1_i1",
                 "TRINITY_DN129076_c1_g1_i1",
                 "TRINITY_DN45900_c0_g1_i1",
                 "TRINITY_DN127729_c5_g3_i2",
                 "TRINITY_DN124927_c0_g1_i8",
                 "TRINITY_DN124927_c0_g1_i3",
                 "TRINITY_DN120193_c1_g1_i2",
                 "TRINITY_DN120856_c1_g2_i1",
                 "TRINITY_DN127729_c5_g3_i1",
                 "TRINITY_DN120856_c1_g1_i3",
                 "TRINITY_DN120856_c1_g1_i11",
                 "TRINITY_DN127729_c5_g9_i1",
                 "TRINITY_DN120193_c1_g2_i1",
                 "TRINITY_DN117914_c5_g2_i1",
                 "TRINITY_DN10791_c0_g1_i1",
                 "TRINITY_DN56951_c0_g1_i1",
                 "TRINITY_DN121353_c0_g1_i2",
                 "TRINITY_DN121353_c0_g1_i1",
                 "TRINITY_DN121353_c0_g1_i7",
                 "TRINITY_DN121353_c0_g1_i15",
                 "TRINITY_DN127729_c5_g1_i2",
                 "TRINITY_DN127729_c5_g1_i1",
                 "TRINITY_DN121353_c0_g1_i8",
                 "TRINITY_DN121353_c0_g1_i10")
ft1_blast_onion <- 
  diff_expression[diff_expression$transcript_id %in% ft1_blast_onion_names,]

ft4_blast_onion_names <- c("TRINITY_DN120856_c1_g2_i2",
                           "TRINITY_DN124927_c0_g2_i1",
                           "TRINITY_DN120193_c1_g1_i3",
                           "TRINITY_DN120193_c1_g1_i1",
                           "TRINITY_DN127729_c5_g3_i3",
                           "TRINITY_DN124927_c0_g1_i1",
                           "TRINITY_DN124927_c0_g1_i7",
                           "TRINITY_DN124927_c0_g1_i6",
                           "TRINITY_DN124927_c0_g1_i4",
                           "TRINITY_DN124927_c0_g1_i2",
                           "TRINITY_DN124927_c0_g1_i9",
                           "TRINITY_DN111334_c0_g1_i1",
                           "TRINITY_DN129076_c1_g1_i1",
                           "TRINITY_DN45900_c0_g1_i1",
                           "TRINITY_DN127729_c5_g3_i2",
                           "TRINITY_DN117914_c5_g2_i1",
                           "TRINITY_DN124927_c0_g1_i8",
                           "TRINITY_DN127729_c5_g3_i1",
                           "TRINITY_DN120856_c1_g2_i1",
                           "TRINITY_DN127729_c5_g9_i1",
                           "TRINITY_DN124927_c0_g1_i3",
                           "TRINITY_DN120193_c1_g1_i2",
                           "TRINITY_DN120856_c1_g1_i3",
                           "TRINITY_DN120856_c1_g1_i11",
                           "TRINITY_DN10791_c0_g1_i1",
                           "TRINITY_DN120193_c1_g2_i1",
                           "TRINITY_DN56951_c0_g1_i1",
                           "TRINITY_DN118335_c0_g1_i4",
                           "TRINITY_DN118335_c0_g1_i1",
                           "TRINITY_DN127729_c5_g1_i2",
                           "TRINITY_DN127729_c5_g1_i1")

ft4_blast_onion <- 
  diff_expression[diff_expression$transcript_id %in% ft4_blast_onion_names,]

all_fts <- unique(c(ft4_blast_onion_names, ft1_blast_onion_names, pot_blast))

# summarize the non-group candidate genes

gene_names <- c("OsbHLH120","IDD5","WOX4","Sulfite reductase","FT-like")
taxon <- c("Oryza sativa","Ipomea batatas","Ipomea batatas",
           "Manihot esculenta","Solanum tuberosum, Allium cepa")
copies <- c(22,63,31,9,37)
de <- c("1*","0","0","0","1")

data.frame(gene_names = gene_names,
           taxon = taxon,
           copies = copies,
           de = de) %>%
  kable("latex",longtable = F, booktabs = T,
        col.names = c("Gene Names","Original Taxon",
                      "Isoforms Found", "Number DE")) %>%
  row_spec(0, bold = TRUE, italic = FALSE) 


  
