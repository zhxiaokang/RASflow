# load the libraries
if (!require("plotscale")) install.packages('scripts/plotscale_0.1.6.tar.gz', repos = NULL, type="source")
library(yaml)
library(mygene)
library(EnhancedVolcano)
library(plotscale)

# ====================== load parameters in config file ======================

# load the config file
yaml.file <- yaml.load_file('configs/config_main.yaml')

# extract the information from the yaml file
controls <- yaml.file$CONTROL
treats <- yaml.file$TREAT

# passing the params from command line
args <- commandArgs(TRUE)
norm.path <- args[1]
dea.path <- args[2]
out.path <- args[3]

# check the number of comparisons
num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

# function to plot volcano plot and heatmap
plot.volcano.heatmap <- function(name.control, name.treat) {
  file.dea.table <- paste(dea.path, "/dea_", name.control, "_", name.treat, ".tsv", sep = "")
  norm.control <- paste(norm.path, "/", name.control, "_gene_norm.tsv", sep = "")  # normalized table of control
  norm.treat <- paste(norm.path, "/", name.treat, "_gene_norm.tsv", sep = "")  # normalized table of treat

  dea.table <- read.table(file.dea.table, header = TRUE, row.names = 1)
  # sort the dea table: ascending of FDR then descending of absolute valued of logFC
  dea.table <- dea.table[order(dea.table$FDR, -abs(dea.table$logFC), decreasing = FALSE), ]

  gene.id.dea <- row.names(dea.table)

  gene.symbol.dea <- queryMany(gene.id.dea, scopes = 'ensembl.gene', fields = 'symbol')$symbol

  # if can't find a symbol for the id, then keep the id as it is
  gene.dea <- gene.symbol.dea
  for (i in c(1:length(gene.dea))) {
    if (is.na(gene.dea[i])) {
      gene.dea[i] <- gene.id.dea[i]
    }
  }

  # volcano plot
  fig.volcano <- EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5,
                                 FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5, title = 'Volcano plot for DEA', subtitle = NULL)
  as.pdf(fig.volcano, width = 8, height = 5, scaled = TRUE, file = file.path(out.path, paste('volcano_plot_', name.control, '_', name.treat, '.pdf', sep = '')))

  # heatmap
  norm.table.control <- read.table(norm.control, header = TRUE, row.names = 1)
  norm.table.treat <- read.table(norm.treat, header = TRUE, row.names = 1)

  num.control <- dim(norm.table.control)[2]
  num.treat <- dim(norm.table.treat)[2]

  norm.table <- cbind(norm.table.control, norm.table.treat)
  groups <- c(name.control, name.treat)

  # instead using all genes, only use the top 50 genes in dea.table
  index.deg <- which(row.names(norm.table) %in% gene.id.dea[1:50])
  norm.table.deg <- norm.table[index.deg,]

  gene.id.norm.table <- rownames(norm.table.deg)

  gene.symbol.norm.table <- queryMany(gene.id.norm.table, scopes = 'ensembl.gene', fields = 'symbol')$symbol

  # if can't find a symbol for the id, then keep the id as it is
  gene.norm.table <- gene.symbol.norm.table
  for (i in c(1:length(gene.norm.table))) {
    if (is.na(gene.norm.table[i])) {
      gene.norm.table[i] <- gene.id.norm.table[i]
    }
  }

  palette <- c("#4DAF4A", "#377EB8")
  palette.group <- c(rep(palette[1], num.control), rep(palette[2], num.treat))

  ## draw heatmap
  pdf(file = file.path(out.path, paste('heatmap_', name.control, '_', name.treat, '.pdf', sep = '')), width = 15, height = 12, title = 'Heatmap using the top features')
  heatmap(as.matrix(norm.table.deg), ColSideColors = palette.group, margins = c(8,5), labRow = gene.norm.table)
  legend("topleft", title = 'Group', legend=groups, text.font = 2,
         col = palette, fill = palette, cex=0.8)
  dev.off()
}

# the main function
for (ith.comparison in c(1:num.comparison)) {
  name.control <- controls[ith.comparison]
  name.treat <- treats[ith.comparison]
  plot.volcano.heatmap(name.control, name.treat)
}
