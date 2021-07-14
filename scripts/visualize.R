# load the libraries
if (!require("plotscale")) install.packages('scripts/plotscale_0.1.6.tar.gz', repos = NULL, type="source")
library(yaml)
library(hash)
library(mygene)
library(EnhancedVolcano)
library(plotscale)

# ====================== load parameters in config file ======================

# passing the params from command line
args <- commandArgs(TRUE)
norm.path <- args[1]
dea.path <- args[2]
out.path <- args[3]

# load the config file
if (length(args) > 3) {  # this script is used in  visualize_test.rules
    yaml.file <- yaml.load_file(args[4])
} else {  # this script is used in visualize.rules
    yaml.file <- yaml.load_file('configs/config_main.yaml')
}

# extract the information from the yaml file
controls <- yaml.file$CONTROL
treats <- yaml.file$TREAT
dea.tool <- yaml.file$DEATOOL

# check the number of comparisons
num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

convert.id2symbol <- function(gene.id) {
  gene.symbol <- gene.id  # initialize the gene symbol with the gene id

  # it may happen that no symbol can be found for any id. In that case, "queryMany" will throw an error
  # so "try" is used here to take care of that error
  try({
    gene.symbol.all <- queryMany(gene.id, scopes = 'ensembl.gene', fields = 'symbol')

    h <- hash()
    for (i in 1:nrow(gene.symbol.all)) {
      query <- gene.symbol.all$query[i]
      symbol <- gene.symbol.all$symbol[i]
      if (has.key(query, h)) {  # if there's duplicate for the same query
        h[[query]] <- paste(hash::values(h, keys = query), symbol, sep = ', ')
      } else {
        if (is.na(symbol)) {  # if there's no hit for the query, keep the original id
          h[[query]] <- query
        } else {
          h[[query]] <- symbol
        }
      }
    }
    
    for (i in c(1:length(gene.symbol))) {
      gene.symbol[i] <- h[[gene.id[i]]]
    }
  })

  return(gene.symbol)
}

# function to plot volcano plot and heatmap
plot.volcano.heatmap <- function(name.control, name.treat) {
  file.dea.table <- paste(dea.path, "/dea_", name.control, "_", name.treat, ".tsv", sep = "")
  norm.control <- paste(norm.path, "/", name.control, "_gene_norm.tsv", sep = "")  # normalized table of control
  norm.treat <- paste(norm.path, "/", name.treat, "_gene_norm.tsv", sep = "")  # normalized table of treat

  dea.table <- read.table(file.dea.table, header = TRUE, row.names = 1)
  # sort the dea table: ascending of FDR then descending of absolute valued of logFC
  if (dea.tool == 'edgeR') {
    dea.table <- dea.table[order(dea.table$FDR, -abs(dea.table$logFC), decreasing = FALSE), ]  
  } else if (dea.tool == 'DESeq2') {
    dea.table <- dea.table[order(dea.table$padj, -abs(dea.table$log2FoldChange), decreasing = FALSE), ]
  }

  gene.id.dea <- row.names(dea.table)

  # convert the gene id to gene symbol, and keep the ID if no symbol can be found for that ID
  gene.dea <- convert.id2symbol(gene.id.dea)

  # volcano plot
  dea.table.volcano <- dea.table  # for better volcano plot, 0 FDRs/padj will be changed to a very low value
  
  if (dea.tool == 'edgeR') {
    # change the 0 FDR to a low value (100 times smaller than the minumum non-zero value)
    FDR <- dea.table.volcano$FDR
    FDR.min.non.zero <- min(FDR[FDR>0])
    dea.table.volcano$FDR[FDR==0] <- FDR.min.non.zero/100
    ## define the range of x-axis and y-axis
    log2FC_lim <- max(abs(dea.table.volcano$logFC))
    FDR_lim <- -log10(min(dea.table.volcano$FDR)) # NAs already removed from dea.table

    fig.volcano <- EnhancedVolcano(dea.table.volcano, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
                                   FCcutoff = 1, xlim = c(-log2FC_lim-1, log2FC_lim+1), ylim = c(0, FDR_lim), title = NULL, subtitle = NULL)
  } else if (dea.tool == 'DESeq2') {
    # change the 0 padj to a low value (100 times smaller than the minumum non-zero value)
    padj <- dea.table.volcano$padj
    padj.min.non.zero <- min(padj[padj>0])
    dea.table.volcano$padj[padj==0] <- padj.min.non.zero/100
    ## define the range of x-axis and y-axis
    log2FC_lim <- max(abs(dea.table.volcano$log2FoldChange))
    padj_lim <- -log10(min(dea.table.volcano$padj)) # NAs already removed from dea.table

    fig.volcano <- EnhancedVolcano(dea.table.volcano, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'log2FoldChange', y = 'padj', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
                                   FCcutoff = 1, xlim = c(-log2FC_lim-1, log2FC_lim+1), ylim = c(0, padj_lim), title = NULL, subtitle = NULL)
  }
  
  # ## if the range of x-axis or y-axis gets crazy, you may also manually define it to show a subset of the genes (but to be noted, the genes out of range will not be shown)
  # if (dea.tool == 'edgeR') {
  #   fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
  #                                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5, title = NULL, subtitle = NULL)
  # } else if (dea.tool == 'DESeq2') {
  #   fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'log2FoldChange', y = 'padj', pCutoff = 10e-5, col = c("grey30", "orange2", "royalblue", "red2"),
  #                                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5, title = NULL, subtitle = NULL)
  # }
  
  as.pdf(fig.volcano, width = 9, height = 6, scaled = TRUE, file = file.path(out.path, paste('volcano_plot_', name.control, '_', name.treat, '.pdf', sep = '')))

  # heatmap
  norm.table.control <- read.table(norm.control, header = TRUE, row.names = 1)
  norm.table.treat <- read.table(norm.treat, header = TRUE, row.names = 1)

  num.control <- dim(norm.table.control)[2]
  num.treat <- dim(norm.table.treat)[2]

  norm.table <- cbind(norm.table.control, norm.table.treat)
  groups <- c(name.control, name.treat)

  # instead using all genes, only use the top 20 genes in dea.table
  index.deg <- which(row.names(norm.table) %in% gene.id.dea[1:20])
  norm.table.deg <- norm.table[index.deg,]

  gene.id.norm.table <- rownames(norm.table.deg)

  # convert the gene id to gene symbol, and keep the ID if no symbol can be found for that ID
  gene.norm.table <- convert.id2symbol(gene.id.norm.table)

  palette <- c("#999999", "#377EB8")
  palette.group <- c(rep(palette[1], num.control), rep(palette[2], num.treat))

  ## draw heatmap
  pdf(file = file.path(out.path, paste('heatmap_', name.control, '_', name.treat, '.pdf', sep = '')), width = 15, height = 12, title = 'Heatmap using the top features')
  heatmap(as.matrix(norm.table.deg), ColSideColors = palette.group, margins = c(9,5.5), labRow = gene.norm.table, cexRow = 1.9, cexCol = 1.9)
  legend("topleft", title = 'Group', legend=groups, text.font = 15,
         col = palette, fill = palette, cex=1.8)
  dev.off()
}

# the main function
for (ith.comparison in c(1:num.comparison)) {
  name.control <- controls[ith.comparison]
  name.treat <- treats[ith.comparison]
  plot.volcano.heatmap(name.control, name.treat)
}
