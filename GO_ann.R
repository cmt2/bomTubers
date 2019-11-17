# pull out all cellular components and all molecular functions from GO annotations
library(knitr)
library(kableExtra)
library(tidyverse)
setwd("~/Desktop/trans_diff")
GO <- read.csv("res.csv", stringsAsFactors = FALSE)$gene_ontology_blast

cell_components <- character()
mol_functions <- character()
bio_process <- character()

for (i in 1:length(GO)){
  if (!is.na(GO[i])){
    
    l <- unlist(strsplit(GO[i],  "`"))
    
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

top_GO <- data.frame(
top_cell_comp = round(sort(table(cell_components)/length(mol_functions), 
                           decreasing = TRUE), digits = 3)[1:10],
top_mol_fun = round(sort(table(mol_functions)/length(mol_functions), 
                         decreasing = TRUE), digits = 3)[1:10],
top_bio_pro = round(sort(table(bio_process)/length(mol_functions), 
                         decreasing = TRUE), digits = 3)[1:10])

kable(top_GO, "latex", longtable = F, booktabs = T, 
      col.names = c(rep(c("Annotation", "Frequency"),3))) %>%
  add_header_above((c("Cellular Components" = 2, 
                      "Molecular Functions" = 2,
                      "Biological Processes" = 2)),bold = T) %>%
  row_spec(0, bold=TRUE) %>%
  kable_styling(latex_options = c("repeat_header"))

### Do the same but for the top 10 genes, keeping gene name intact 
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

