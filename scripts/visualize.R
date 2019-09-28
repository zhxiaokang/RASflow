# load the libraries
library(yaml)
library(mygene)
library(EnhancedVolcano)
library(RColorBrewer)
library(data.table)

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

gene.id.norm <- row.names(norm.table)
gene.id.dea <- row.names(dea.table)
gene.id.deg <- row.names(deg.table)

gene.symbol.norm <- queryMany(gene.id.norm, scopes = 'ensembl.gene', fields = 'symbol')$symbol
gene.symbol.dea <- queryMany(gene.id.dea, scopes = 'ensembl.gene', fields = 'symbol')$symbol
gene.symbol.deg <- queryMany(gene.id.deg, scopes = 'ensembl.gene', fields = 'symbol')$symbol

# if can't find a symbol for the id, then keep the id as it is
gene.dea <- gene.symbol.dea
for (i in c(1:length(gene.dea))) {
  if (is.na(gene.dea[i])) {
    gene.dea[i] <- gene.id.dea[i]
  }
}

# volcano plot
pdf(file = paste(output.path, '/volcano_plot.pdf'), width = 6, height = 3.5)
EnhancedVolcano(dea.table, lab = gene.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5,
                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5)
dev.off()

# heatmap
## collect all the files
files <- file.path(count.path, paste(groups, "_norm.csv", sep = ''))
norm.tables <- lapply(files, fread)
norm.table.temp <- do.call(merge, norm.tables)
norm.table.temp <- as.data.frame(norm.table.temp)
norm.table <- norm.table.temp[, -1]
row.names(norm.table) <- unlist(norm.table.temp[, 1])

index.deg <- which(row.names(norm.table) %in% gene.id.deg)
norm.table.deg <- norm.table[index.deg,]
palette <- brewer.pal(n = length(levels(group.all)), name = "Set1")
palette.group <- group.all
for (i in c(1:length(levels(palette.group)))) {
  levels(palette.group)[levels(palette.group) == levels(palette.group)[i]] <- palette[i]
}
## draw heatmap
pdf(file = paste(output.path, '/heatmap.pdf'), width = 6, height = 3.5)
heatmap(as.matrix(norm.table.deg), ColSideColors = as.character(palette.group))
dev.off()