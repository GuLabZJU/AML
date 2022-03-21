

basedir <- "/Users/yinjie/Project/TARGET_AML/HTseq_count"

gdc_file <- read.csv("/Users/yinjie/Project/TARGET_AML/HTseq_count/gdc_manifest.2020-01-14.txt",sep="\t",header = T)


listfiles <- c()
for (i in 1: length(gdc_file[,1])){
  listfiles[i] <- paste(basedir, gdc_file[i,1],gdc_file[i,2], sep="/")
}

# tophat.all <- list.files(path = rawfiles, 
#                          pattern = "_all_counts.txt",
#                          all.files = TRUE, 
#                          recursive = FALSE, 
#                          ignore.case = FALSE, 
#                          include.dirs = FALSE)

# we choose the 'all' series
myfiles <- listfiles
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# for (i in 1:length(myfiles) ) {
#   infile = paste("htseq_counts", myfiles[i], sep = "/")
#   DT[[myfiles[i]]] <- read.table(infile, header = F)
#   cnts <- gsub("(.*)_all_mdup_counts.txt", "\\1", myfiles[i])
#   colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
# }

# for (i in length(myfiles) ) {
#   
#   DT[[myfiles[i]]] <- read.table(gzfile(myfiles[i]), header = F)
#   cnts <- gsub("\\/Users\\/yinjie\\/Project\\/TARGET_AML\\/HTseq_count\\/(.*)\\/(.*).htseq.counts.gz", "\\2", myfiles[i])
#   colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
# }

for (i in 1:length(myfiles) ) {
  DT[[myfiles[i]]] <- read.table(gzfile(myfiles[i]), header = F)
  cnts <- gsub("\\/Users\\/yinjie\\/Project\\/TARGET_AML\\/HTseq_count\\/(.*)\\/(.*).htseq.counts.gz", "\\2", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}


# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}




# remove top two rows
data.all.summary <- data[grep("^ENS", data$ID, perl=TRUE, invert=TRUE), ]
rownames(data.all.summary) <- data.all.summary$ID
data.all.summary <- data.all.summary[,-1]




data.all <- data[grep("^ENS", data$ID, perl=TRUE, invert=FALSE), ]
rownames(data.all) <- data.all$ID
data.all <- data.all[,-1]


# final merged table
head(data.all, 3)

# summary counts by type
t(data.all.summary)

write.table(data.all,
            file = "TARGET_AML_HTSeq_count.csv",
            append = FALSE,
            quote = FALSE,
            sep = " ",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = TRUE,
            col.names = TRUE,
            qmethod = c("escape", "double"))

write.table(data.all.summary,
            file = "hTARGET_AML_HTSeq_count-summary.csv",
            append = FALSE,
            quote = FALSE,
            sep = " ",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = TRUE,
            col.names = TRUE,
            qmethod = c("escape", "double"))
save.image("TARGET_AML_HTSeq-Count.RData")

head(data)

gene_name <- gsub("(.*)\\.\\d+","\\1",rownames(data.all))
rownames(data.all) <- gene_name

length(unique(gene_name))

#retrieve gene length
library(GenomicFeatures)
#txdb <- makeTxDbFromGFF('Homo_sapiens.GRCh38.98.gtf')
txdb <- makeTxDbFromGFF('/Users/yinjie/Database/Homo_sapiens.GRCh38.98.gtf')

## get gene information
all.genes <- genes(txdb)
## import your list of gene names
my.genes <-  gene_name

library(data.table)
library(tibble)
library(dplyr)

data.all.tb <- as_tibble(data.all)
data.all.tb <- mutate(data.all.tb,ID=rownames(data.all))

#filter unexisted gene id
data.all.sub <- filter(data.all.tb, ID %in% all.genes$gene_id)
data.all.sub.df <- data.all.sub %>% as.data.frame()


rownames(data.all.sub.df) <- data.all.sub.df[,188]
data.all.sub.df <- data.all.sub.df[,-188]

## get the length of each of those genes
my.genes.lengths <- width(all.genes[data.all.sub$ID])
gene_length <- as_tibble(my.genes.lengths)
gene_length <- mutate(gene_length,ID=data.all.sub$ID)
colnames(gene_length) <- c("Length","ID")
gene_length <- as.data.frame(gene_length)
rownames(gene_length) <- gene_length$ID

source("c2tpm.R")

data_normalized <- c2tpm(data.all.sub.df,gene_length$Length)

data_normalized <- as.data.frame(data_normalized)

data_normalized.cd <- data_normalized[c("ENSG00000123374","ENSG00000111640"),]
#CDK2
boxplot(log2(data_normalized["ENSG00000123374",]+1))
#GAPDH
boxplot(data_normalized["ENSG00000111640",])


write.table(data_normalized.cd,"AML_CDK2vsGAPDH.tsv",quote=FALSE,sep="\t",col.names = T,row.names = T)


save.image("TARGET_AML_HTSeq-Count.Normalized.RData")
