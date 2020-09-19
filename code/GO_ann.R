# pull out all cellular components and all molecular functions from GO annotations
library(knitr)
library(kableExtra)
library(tidyverse)
setwd("~/Documents/trans_diff")

# identify go categories that are overrepresented in the deg genes compared to full transcriptome
go_degs <- read.csv("res.csv", stringsAsFactors = FALSE)$gene_ontology_blast
go_all <- read.csv(file = "~/Documents/trans_diff/trans_annotate/trinotate_annotation_report.csv",
                   na.strings = ".", stringsAsFactors = FALSE)

parse_go_terms <- function(terms) {
  cell_components <- character()
  mol_functions <- character()
  bio_process <- character()
  for (i in 1:length(terms)){
    if (!is.na(terms[i])){
      
      l <- unlist(strsplit(terms[i],  "`"))
      
      for (j in 1:length(l)){
        
        s <- unlist(strsplit(l[j], "\\^"))
        
        if (s[2] == "cellular_component"){
          cell_components <- append(cell_components, s[3])
        } else if (s[2]  == "molecular_function"){
          mol_functions <- append(mol_functions, s[3])
        } else if (s[2] == "biological_process"){
          bio_process <- append(bio_process, s[3])
        }
      }
    }
  }
  return(list(cell_components = cell_components,
              mol_functions = mol_functions,
              bio_process =  bio_process))
}

go_terms_degs <- parse_go_terms(terms = go_degs)
go_terms_all <- parse_go_terms(terms = go_all$gene_ontology_blast)

##### bio process #####
degs_go_table_bio <- data.frame(table(go_terms_degs$bio_process))
colnames(degs_go_table_bio) <- c("bio_process", "deg_count")
all_go_table_bio <- data.frame(table(go_terms_all$bio_process))
colnames(all_go_table_bio) <- c("bio_process", "all_transcript_count")
bio_process_combined <- left_join(all_go_table_bio, degs_go_table_bio)
bio_process_combined$deg_count[is.na(bio_process_combined$deg_count)] <- 0
bio_process_combined$ref <- bio_process_combined$all_transcript_count - 
                            bio_process_combined$deg_count

# test if distributions are different
test_bio <- fisher.test(bio_process_combined[,c("deg_count","ref")], simulate.p.value = T, B = 10000)

# find individual GO categories that are enriched 
calc_GO_p_value <- function(deg_count, full_count, 
                            deg_size = length(go_degs), 
                            full_size = length(go_all$gene_ontology_blast), 
                            num_compare ) {
  p <- pbinom( deg_count, size = deg_size, prob = (full_count - deg_count) / (full_size - deg_size), lower.tail = T)
  p_0 <- dbinom( 0, size = deg_size, prob = (full_count -deg_count) / (full_size - deg_size))
  p_adj <- (1 - (p - p_0)) * num_compare
  return(p_adj)
}
  
p_values <- numeric()
GO_terms <- character()
for (i in 1:nrow(bio_process_combined)) {
  if (bio_process_combined$deg_count[i] > 0) {
    p_values[i] <- calc_GO_p_value(deg_count = bio_process_combined$deg_count[i],
                                   full_count = bio_process_combined$all_transcript_count[i],
                                   num_compare = length(go_terms_degs$bio_process))
    GO_terms[i] <- as.character(bio_process_combined$bio_process[i])
  }
}

bio_processes_sig <- data.frame(GO_bio = GO_terms[which(p_values < 0.05)],
                                p_values = p_values[which(p_values < 0.05)],
                                count_in_degs = bio_process_combined$deg_count[which(p_values < 0.05)],
                                count_in_full = bio_process_combined$all_transcript_count[which(p_values < 0.05)])

##### mol function #####
degs_go_table_mol <- data.frame(table(go_terms_degs$mol_functions))
colnames(degs_go_table_mol) <- c("mol_function", "deg_count")
all_go_table_mol <- data.frame(table(go_terms_all$mol_functions))
colnames(all_go_table_mol) <- c("mol_function", "all_transcript_count")
mol_function_combined <- left_join(all_go_table_mol, degs_go_table_mol)
mol_function_combined$deg_count[is.na(mol_function_combined$deg_count)] <- 0
mol_function_combined$ref <- mol_function_combined$all_transcript_count - 
  mol_function_combined$deg_count

# test if distributions are different
test_mol <- fisher.test(mol_function_combined[,c("deg_count","ref")], simulate.p.value = T, B = 10000)

# find individual GO categories that are enriched 
p_values <- numeric()
GO_terms <- character()
for (i in 1:nrow(mol_function_combined)) {
  if (mol_function_combined$deg_count[i] > 0) {
    p_values[i] <- calc_GO_p_value(deg_count = mol_function_combined$deg_count[i],
                                   full_count = mol_function_combined$all_transcript_count[i],
                                   num_compare = length(go_terms_degs$mol_functions))
    GO_terms[i] <- as.character(mol_function_combined$mol_function[i])
  }
}

mol_functions_sig <- data.frame(GO_mol = GO_terms[which(p_values < 0.05)],
                                p_values = p_values[which(p_values < 0.05)],
                                count_in_degs = mol_function_combined$deg_count[which(p_values < 0.05)],
                                count_in_full = mol_function_combined$all_transcript_count[which(p_values < 0.05)])


##### cellular components #####
degs_go_table_cel <- data.frame(table(go_terms_degs$cell_components))
colnames(degs_go_table_cel) <- c("cell_comp", "deg_count")
all_go_table_cel <- data.frame(table(go_terms_all$cell_components))
colnames(all_go_table_cel) <- c("cell_comp", "all_transcript_count")
cell_comp_combined <- left_join(all_go_table_cel, degs_go_table_cel)
cell_comp_combined$deg_count[is.na(cell_comp_combined$deg_count)] <- 0
cell_comp_combined$ref <- cell_comp_combined$all_transcript_count - 
  cell_comp_combined$deg_count

# test if distributions are different
test_cel <- fisher.test(cell_comp_combined[,c("deg_count","ref")], simulate.p.value = T, B = 10000)

# find individual GO categories that are enriched 

p_values <- numeric()
GO_terms <- character()
for (i in 1:nrow(cell_comp_combined)) {
  if (cell_comp_combined$deg_count[i] > 0) {
    p_values[i] <- calc_GO_p_value(deg_count = cell_comp_combined$deg_count[i],
                                   full_count = cell_comp_combined$all_transcript_count[i],
                                   num_compare = length(go_terms_degs$cell_components))
    GO_terms[i] <- as.character(cell_comp_combined$cell_comp[i])
  }
}

cell_comp_sig <- data.frame(GO_cell = GO_terms[which(p_values < 0.05)],
                            p_values = p_values[which(p_values < 0.05)],
                            count_in_degs = cell_comp_combined$deg_count[which(p_values < 0.05)],
                            count_in_full = cell_comp_combined$all_transcript_count[which(p_values < 0.05)])

##### GO anntoations of the top 10 genes #####
GO_strict <- read.csv("res_strict.csv", stringsAsFactors = FALSE)

transcript <- character()
ann_name <- character()
names <- character()
log2FoldChange <- numeric()
cell_components <- character()
mol_functions <- character()
bio_process <- character()

for (i in 1:length(GO_strict$gene_ontology_blast)){
  transcript[i] <- GO_strict$transcript_id[i]
  ann_name[i] <- strsplit(strsplit(strsplit(strsplit(GO_strict$sprot_Top_BLASTP_hit[i], "\\^")[[1]][6], "=")[[1]][2], ";")[[1]][1], "\\{")[[1]][1]
  names[i] <- paste0(transcript[i], ": ", ann_name[i])
  log2FoldChange[i] <- round(GO_strict$log2FoldChange[i], digits = 2)
  c <- character()
  m <- character()
  b <- character() 
  if (!is.na(GO_strict$gene_ontology_blast[i])){
    l <- unlist(strsplit(GO_strict$gene_ontology_blast[i],  "`"))
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

data.frame(names = names,
           log2FoldChange = log2FoldChange,
           cell_components = cell_components, 
           mol_functions = mol_functions,
           bio_process = bio_process) %>%
  arrange(-row_number()) %>%
  kable("latex", longtable = F, booktabs = T, align=c('l','c',rep('l',3)),
        col.names = c("Transcript ID: Annotated Name (SPROT)",
                      "Log2 Fold Change",
                      "Cellular Components", "Molecular Functions", 
                      "Biological Processes")) %>%
  add_header_above(c("","", "Gene Ontology" = 3),
                   bold = T) %>%
  column_spec(1, width = "15em", bold = T) %>%
  kable_styling(latex_options = c("striped") ) %>%
  row_spec(0, bold=TRUE) %>%
  column_spec(2, width = "5em") %>%
  column_spec(3, width = "25em") %>%
  column_spec(4, width = "25em") %>%
  column_spec(5, width = "25em") 

