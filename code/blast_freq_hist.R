perc_bins <- seq(1:10)*10
num_genes <- rev(c(6217,2210,1368,1102,1142,1209,1344,1623,1592,504))
gene_freq <- num_genes/sum(num_genes) *100

xx <- barplot(num_genes, names.arg = perc_bins, 
        xlab = "Percent sequence length identity",
        ylab = "Number of genes in bin", ylim = c(0,max(num_genes)*1.25))
text(x = xx, y = num_genes + 50 , label = paste0(round(gene_freq, digits = 2), "%"), 
     pos = 3, cex = 0.6, col = "red")