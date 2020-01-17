# combine quant of samples to group

library(biomaRt)
library(yaml)
library(tximport)

# ====================== define some functions ======================

# remove the version in the transcript ID
remove_version <- function(files) {  # input files (file names with directory) are output from Salmon
  for (i in c(1:length(files))) {
    quant.file <- files[i]
    quant.table <- read.table(quant.file, header = TRUE, stringsAsFactors = FALSE)
    trans.id.version <- quant.table$Name
    trans.id <- rep('ID', length(trans.id.version))
    for (j in c(1:length(trans.id.version))) {
      trans.id[j] <- strsplit(trans.id.version[j], ".", fixed = TRUE)[[1]][1]
    }
    quant.table$Name <- trans.id
    write.table(quant.table, file.path(input.path, samples[i], "quant_noVersion.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  return(trans.id)
}

# ====================== load parameters in config file ======================
# load the config file
yaml.file <- yaml.load_file('configs/config_main.yaml')

# extract the information from the yaml file
project <- yaml.file$PROJECT  # project name
input.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/quant")
gene.level <- yaml.file$GENE_LEVEL
controls <- yaml.file$CONTROL  # all groups used as control
treats <- yaml.file$TREAT  # all groups used as treat, should correspond to control
filter.need <- yaml.file$FILTER$yesOrNo
cpm.threshold <- yaml.file$FILTER$cpm
pair.test <- yaml.file$PAIR
meta.file <- yaml.file$METAFILE
dataset <- yaml.file$EnsemblDataSet
output.path <- file.path(yaml.file$FINALOUTPUT, project, "trans/dea")

num.control <- length(controls)  # number of comparisons that the user wants to do
num.treat <- length(treats)  # should equals to num.control

if (num.control != num.treat) {
  message("Error: Control groups don't mathch with treat groups!")
  message("Please check config_dea.yaml")
  quit(save = 'no')
}

# extract the metadata
meta.data <- read.table(meta.file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
samples <- meta.data$sample
group.all <- meta.data$group
subject.all <- meta.data$subject

# ====================== prepare the quant files ======================

# the original quant files from Salmon
files <- file.path(input.path, samples, "quant.sf")
names(files) <- samples

# get the normalized abundance tables for all samples (scaled using the average transcript length over samples and then the library size (lengthScaledTPM))

# ====================== prepare the tx2gene table ======================
if (gene.level) {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = dataset)
    datasets <- listDatasets(ensembl)

    attributes <- listAttributes(mart = ensembl)

    # remove the version to get more matches
    trans.id <- remove_version(files)

    files.noVersion <- file.path(input.path, samples, "quant_noVersion.sf")
    names(files.noVersion) <- samples

    tx2gene <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id'),
                    filters = 'ensembl_transcript_id', values = trans.id, mart = ensembl)
    # save tx2gene
    output.file.tx2gene <- file.path(output.path, "countGroup", 'tx2gene.RData')
    save(tx2gene, file = output.file.tx2gene)
}
# ====================== get raw and normalized abundance tables ======================

trans.matrix <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no")
trans.count <- trans.matrix$counts

trans.matrix.tpm <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
trans.count.tpm <- trans.matrix.tpm$counts

if (gene.level) {
    gene.matrix <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
    gene.count <- gene.matrix$counts

    gene.matrix.tpm <- tximport(files.noVersion, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    gene.count.tpm <- gene.matrix.tpm$counts
}

# ====================== combine the samples to groups ======================

for (group in unique(group.all)) {
  index.group = which(group.all == group)  # all the index of this group
  samples.group = samples[index.group]  # the samples belonging to this group
  
  group.count.trans <- trans.count[, index.group]
  group.count.trans.tpm <- trans.count.tpm[, index.group]

  if (gene.level) {
    group.count.gene <- gene.count[, index.group]
    group.count.gene.tpm <- gene.count.tpm[, index.group]
  }

  # write to files
  output.file.trans <- file.path(output.path, "countGroup", paste(group, "_trans_abundance.tsv", sep = ""))
  write.table(group.count.trans, output.file.trans, sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

  output.file.trans.norm <- file.path(output.path, "countGroup", paste(group, "_trans_norm.tsv", sep = ""))
  write.table(group.count.trans.tpm, output.file.trans.norm, sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

  if (gene.level) {
    output.file.gene <- file.path(output.path, "countGroup", paste(group, "_gene_abundance.tsv", sep = ""))
    write.table(group.count.gene, output.file.gene, sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

    output.file.gene.norm <- file.path(output.path, "countGroup", paste(group, "_gene_norm.tsv", sep = ""))
    write.table(group.count.gene.tpm, output.file.gene.norm, sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
}


