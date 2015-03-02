setwd("/Users/elijahklowe/Desktop/collaborations/Wei/")
samples <- c("handrdnfgfr15hpf1_GGACCC_L006","handrdnfgfr15hpf2_CCTCGG_L006","handrlaz15hpf1_AAGGGA_L006","handrlaz15hpf2_TTCAGC_L006","handrnoggin15hpf1_CCTTCA_L006","handrnoggin15hpf2_AAGACG_L006")

read.sample <- function(sample.name) {
  file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"))
}
sample.1 <- read.sample(samples[1])
head(sample.1)
nrow(sample.1)

#Read the second sample
sample.2 <- read.sample(samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(sample.1) == nrow(sample.2)
all(sample.1$gene == sample.2$gene)

#Now let's combine them all into one dataset
all.data <- sample.1
all.data <- cbind(sample.1, sample.2$count)
for (c in 3:length(samples)) {
  temp.data <- read.sample(samples[c])
  all.data <- cbind(all.data, temp.data$count)
}

#We now have a data frame with all the data in it:
head(all.data)

colnames(all.data)[2:ncol(all.data)] <- samples

all.data <- all.data[1:(nrow(all.data)-5),]

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
library("DESeq")

#Remove the first column
raw.deseq.data <- all.data[,2:ncol(all.data)]
#Set row names to the gene names
rownames(raw.deseq.data) <- all.data$gene

head(raw.deseq.data)

#Create metadata
handr.design <- data.frame(
  row.names=samples,
  condition=c(rep("fgfr", 2), rep("lacz", 2), rep("noggin", 2)),
  libType=rep("single-end", 6)
)
#Double check it...
handr.design

cds = newCountDataSet(raw.deseq.data,handr.design$condition)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds=estimateDispersions(cds)

plotDispEsts(cds)

res1=nbinomTest(cds, "lacz","fgfr")
plotMA(res1, main="FGFR")
res2=nbinomTest(cds, "lacz","noggin")
plotMA(res2, main="noggin")
res3=nbinomTest(cds, "fgfr","noggin")
plotMA(res3, main="bFGFR vs noggin")
res1Sig=res1[res1$padj < 0.15,]
res1.5Sig=res1[res1$padj < 0.25,]
res2Sig=res2[res2$padj < 0.15,]
res3Sig=res3[res3$padj < 0.15,]
head(res1Sig[order(res1Sig$pval),])
head(res2Sig[order(res2Sig$pval),])
head(res3Sig[order(res3Sig$pval),])
