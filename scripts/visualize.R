# load the libraries
if (!require("plotscale")) install.packages('scripts/plotscale_0.1.6.tar.gz', repos = NULL, type="source")
library(yaml)
library(mygene)
library(EnhancedVolcano)
library(RColorBrewer)
library(data.table)
library(plotscale)

# ====================== load parameters in config file ======================

# load the config file
yaml.file <- yaml.load_file('configs/config_visualize.yaml')

# extract the information from the yaml file
dea.table <- yaml.file$DEAFILE
deg.table <- yaml.file$DEGFILE
count.path <- yaml.file$COUNTPATH  # the folder where all the normalized count tables are stored
output.path <- yaml.file$OUTPUTPATH
meta.file <- yaml.file$METAFILE

# extract the metadata
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
group.all <- meta.data$group
groups <- levels(group.all)

dea.table <- read.csv(dea.table, header = TRUE, row.names = 1)
deg.table <- read.csv(deg.table, header = TRUE, row.names = 1)

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
as.pdf(fig.volcano, width = 8, height = 5, scaled = TRUE, file = paste(output.path, '/volcano_plot.pdf', sep = ''))

# heatmap
## collect all the files
files <- file.path(count.path, paste(groups, "_norm.csv", sep = ''))
norm.tables <- lapply(files, fread)
norm.table.temp <- do.call(merge, norm.tables)
norm.table.temp <- as.data.frame(norm.table.temp)
norm.table <- norm.table.temp[, -1]
row.names(norm.table) <- unlist(norm.table.temp[, 1])

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

# replace the rownames
rownames(norm.table.deg) <- gene.norm.table

palette <- brewer.pal(n = length(levels(group.all)), name = "Set1")
palette.group <- group.all
for (i in c(1:length(levels(palette.group)))) {
  levels(palette.group)[levels(palette.group) == levels(palette.group)[i]] <- palette[i]
}
## draw heatmap
pdf(file = paste(output.path, '/heatmap.pdf', sep = ''), width = 15, height = 12, title = 'Heatmap using the top features')
heatmap(as.matrix(norm.table.deg), ColSideColors = as.character(palette.group), margins = c(8,5))
legend("topleft", title = 'Group', legend=groups, text.font = 2,
       col=as.character(levels(palette.group)), fill=as.character(levels(palette.group)), cex=0.8)
dev.off()
