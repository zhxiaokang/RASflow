# load the libraries
if (!require("plotscale")) install.packages('scripts/plotscale_0.1.6.tar.gz', repos = NULL, type="source")
library(yaml)
library(mygene)
library(EnhancedVolcano)
library(plotscale)

# ====================== load parameters in config file ======================

# load the config file
yaml.file <- yaml.load_file('configs/config_visualize.yaml')

# extract the information from the yaml file
file.dea.table <- yaml.file$DEAFILE
file.deg.table <- yaml.file$DEGFILE
norm.control <- yaml.file$NORMFILES$control
norm.treat <- yaml.file$NORMFILES$treat
name.control <- yaml.file$NAME_CONTROL
name.treat <- yaml.file$NAME_TREAT

outpath.volcano <- dirname(file.dea.table)
outpath.heatmap <- dirname(norm.control)

dea.table <- read.table(file.dea.table, header = TRUE, row.names = 1)
deg.table <- read.table(file.deg.table, header = TRUE, row.names = 1)

gene.id.dea <- row.names(dea.table)
gene.id.deg <- row.names(deg.table)

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
as.pdf(fig.volcano, width = 8, height = 5, scaled = TRUE, file = file.path(outpath.volcano, paste('volcano_plot_', name.control, '_', name.treat, '.pdf', sep = '')))

# heatmap
norm.table.control <- read.table(norm.control, header = TRUE, row.names = 1)
norm.table.treat <- read.table(norm.treat, header = TRUE, row.names = 1)

num.control <- dim(norm.table.control)[2]
num.treat <- dim(norm.table.treat)[2]

norm.table <- cbind(norm.table.control, norm.table.treat)
groups <- c(name.control, name.treat)

# instead using all genes, only use the top 50 degs if there are more than 50
if (length(gene.id.deg) > 50) {
  index.deg <- which(row.names(norm.table) %in% gene.id.deg[1:50])
} else {
  index.deg <- which(row.names(norm.table) %in% gene.id.deg)
}
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
palette.group <- as.factor(c(rep(name.control, num.control), rep(name.treat, num.treat)))

# replace the group name with color code
for (i in c(1:length(levels(palette.group)))) {
  levels(palette.group)[levels(palette.group) == levels(palette.group)[i]] <- palette[i]
}

## draw heatmap
pdf(file = file.path(outpath.heatmap, paste('heatmap_', name.control, '_', name.treat, '.pdf', sep = '')), width = 15, height = 12, title = 'Heatmap using the top features')
heatmap(as.matrix(norm.table.deg), ColSideColors = as.character(palette.group), margins = c(8,5), labRow = gene.norm.table)
legend("topleft", title = 'Group', legend=groups, text.font = 2,
       col=as.character(levels(palette.group)), fill=as.character(levels(palette.group)), cex=0.8)
dev.off()
