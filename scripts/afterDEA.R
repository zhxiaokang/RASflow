# load the libraries
library(mygene)
library(EnhancedVolcano)
library(pheatmap)

norm.table <- read.csv('oil_control_oil_low_norm.csv', header = TRUE, row.names = 1)
dea.table <- read.csv('dea_oil_control_oil_low.csv', header = TRUE, row.names = 1)
deg.table <- read.csv('deg_oil_control_oil_low.csv', header = TRUE, row.names = 1)

gene.id.norm <- row.names(norm.table)
gene.id.dea <- row.names(dea.table)
gene.id.deg <- row.names(deg.table)

# gene.symbol.norm <- queryMany(gene.id.norm, scopes = 'ensembl.gene', fields = 'symbol')$symbol
# gene.symbol.dea <- queryMany(gene.id.dea, scopes = 'ensembl.gene', fields = 'symbol')$symbol
# gene.symbol.deg <- queryMany(gene.id.deg, scopes = 'ensembl.gene', fields = 'symbol')$symbol

# volcano plot
EnhancedVolcano(dea.table, lab = gene.id.dea, xlab = bquote(~Log[2]~ "fold change"), x = 'logFC', y = 'FDR', pCutoff = 10e-5,
                FCcutoff = 1, xlim = c(-5, 5), ylim = c(0, 10), transcriptPointSize = 1.5)

# heatmap
index.deg <- which(row.names(norm.table) %in% gene.id.deg)
norm.table.deg <- norm.table[index.deg,]
heatmap(as.matrix(norm.table.deg))

