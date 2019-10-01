library(yaml)
library(edgeR)

# Use edgeR to do DEA

# load the config file
yaml.file <- yaml.load_file('configs/config_dea_genome.yaml')

# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name of this analysis
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
cpm.threshold <- yaml.file$FILTER$cpm
pair.test <- yaml.file$PAIR
meta.file <- yaml.file$METAFILE
output.path <- file.path(yaml.file$OUTPUTPATH, project, "genome/dea")

# extract the metadata
meta.data <- read.csv(meta.file, header = TRUE, sep = '\t')
group.all <- meta.data$group
subject.all <- meta.data$subject

num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

num.comparison <- num.control

DEA <- function(control, treat) {
  count.control <- read.table(paste(output.path, '/countGroup/', control, '_count.tsv', sep = ''), header = TRUE, row.names = 1)
  count.treat <- read.table(paste(output.path, '/countGroup/', treat, '_count.tsv', sep = ''), header = TRUE, row.names = 1)
  count.table <- cbind(count.control, count.treat)  # merge the control and treat tables together
  
  # number of samples in control and treat groups (should be the same if it's a pair test)
  num.sample.control <- ncol(count.control)
  num.sample.treat <- ncol(count.treat)
  
  # save gene list in gene.list for extracting gene names later
  gene.list <- rownames(count.table)
  
  # get the sample id
  samples <- colnames(count.table)
 
  # normalize the two groups and save the normalized count table
  y.control <- DGEList(counts = count.control, genes = gene.list)
  y.treat <- DGEList(counts = count.treat, genes = gene.list)

  y.control <- calcNormFactors(y.control, method="TMM")
  count.table.control.norm <- cpm(y.control)
  write.csv(count.table.control.norm, paste(output.path, '/countGroup/', control, '_norm.csv', sep = ''))

  y.treat <- calcNormFactors(y.treat, method="TMM")
  count.table.treat.norm <- cpm(y.treat)
  write.csv(count.table.treat.norm, paste(output.path, '/countGroup/', treat, '_norm.csv', sep = ''))

  # Put the data into a DGEList object
  y <- DGEList(counts = count.table, genes = gene.list)
  
  # Filtering
  if (filter.need) {
    countsPerMillion <- cpm(y)
    countCheck <- countsPerMillion > cpm.threshold
    keep <- which(rowSums(countCheck) > 1)
    y <- y[keep, ]
  }
  
  # Normalization
  y <- calcNormFactors(y, method="TMM")
  
  # Do DEA !!!
  # define the group
  subject <- factor(subject.all[c(which(group.all == control), which(group.all == treat))])
  group <- factor(group.all[c(which(group.all == control), which(group.all == treat))])
  
  y$samples$subject <- subject
  y$samples$group <- group
  
  # The design matrix
  if (pair.test) {
    design <- model.matrix(~subject+group)
  } else {
    design <- model.matrix(~group)
  }
  rownames(design) <- colnames(y)
  
  # Estimating the dispersion
  
  # estimate the NB dispersion for the dataset
  y <- estimateDisp(y, design, robust = TRUE)
  
  # Differential expression
  
  # determine differentially expressed genes
  # fit genewise glms
  fit <- glmFit(y, design)
  
  # conduct likelihood ratio tests for tumour vs normal tissue differences and show the top genes
  lrt <- glmLRT(fit)
  
  # the DEA result for all the genes
  # dea <- lrt$table
  toptag <- topTags(lrt, n = nrow(y$genes), p.value = 1)
  dea <- toptag$table  # just to add one more column of FDR
  
  # differentially expressed genes
  toptag <- topTags(lrt, n = nrow(y$genes), p.value = 0.05)
  deg <- toptag$table
  
  # save the DEA result and DEGs to files
  write.csv(dea, paste(output.path, '/DEA/dea_', control, '_', treat, '.csv', sep = ''), row.names = F)
  write.csv(deg, paste(output.path, '/DEA/deg_', control, '_', treat, '.csv', sep = ''), row.names = F) 
}

# the main function
for (ith.comparison in c(1:num.comparison)) {
  control <- controls[ith.comparison]
  treat <- treats[ith.comparison]
  DEA(control, treat)
}
