setwd("/Users/elijahklowe/Desktop/collaborations/Wei/DESeq2_round1/")
hand_samples <- c("handrdnfgfr15hpf1_GGACCC_L006","handrdnfgfr15hpf2_CCTCGG_L006","handrlaz15hpf1_AAGGGA_L006","handrlaz15hpf2_TTCAGC_L006","handrnoggin15hpf1_CCTTCA_L006","handrnoggin15hpf2_AAGACG_L006")

read.sample <- function(sample.name) {
  file.name <- paste(sample.name, "_htseq_counts.txt", sep="")
  result <- read.delim(file.name, col.names=c("gene", "count"), sep="\t", colClasses=c("character", "numeric"))
}
handr.1 <- read.sample(hand_samples[1])
head(handr.1)
nrow(handr.1)

#Read the second sample
handr.2 <- read.sample(hand_samples[2])

#Let's make sure the first and second samples have the same number of rows and the same genes in each row
nrow(handr.1) == nrow(handr.2)
all(handr.1$gene == handr.2$gene)

#Now let's combine them all into one dataset
all_hand.data <- handr.1
all_hand.data <- cbind(handr.1, handr.2$count)
for (c in 3:length(hand_samples)) {
  temp.data <- read.sample(hand_samples[c])
  all_hand.data <- cbind(all_hand.data, temp.data$count)
}

#We now have a data frame with all the data in it:
head(all_hand.data)

colnames(all_hand.data)[2:ncol(all_hand.data)] <- hand_samples

all_hand.data <- all_hand.data[1:(nrow(all_hand.data)-5),]

library("DESeq2")

#Remove the first column
raw.hand_deseq.data <- all_hand.data[,2:ncol(all_hand.data)]
#Set row names to the gene names
rownames(raw.hand_deseq.data) <- all_hand.data$gene

head(raw.hand_deseq.data)

#Create metadata
handr.design <- data.frame(
  row.names=hand_samples,
  batch=rep(c("sample1","sample2"),3),
  condition=c(rep("fgfr", 2), rep("lacz", 2), rep("noggin", 2)),
  libType=rep("single-end", 6)
)
#Double check it...
handr.design
handr.design$condition <- relevel(handr.design$condition, ref="lacz")

#DESeq data object
(hand_deseq.data <- DESeqDataSetFromMatrix(countData = raw.hand_deseq.data,
                                      colData = handr.design,
                                      design = ~ batch + condition))

hand_deseq.data <- hand_deseq.data[ rowSums( counts(hand_deseq.data) ) > 0 , ]
hand_dds <- DESeq(hand_deseq.data)
plotDispEsts(hand_dds)
handr_rld <- rlog(hand_dds)
plotPCA( handr_rld, intgroup = c("batch","condition"))

fgfr_lacZ.res <- results( hand_dds, contrast = c("condition", "fgfr", "lacz"),
                          pAdjustMethod = "fdr")
noggin_lacZ.res <- results( hand_dds, contrast = c("condition", "noggin", "lacz"),
                          pAdjustMethod = "fdr")
fgfr_noggin.res <- results( hand_dds, contrast = c("condition", "fgfr", "noggin"),
                          pAdjustMethod = "fdr")

summary(fgfr_lacZ.res)
summary(noggin_lacZ.res)
summary(fgfr_noggin.res)

fgfr_lacZ.resOrdered <- fgfr_lacZ.res[order(fgfr_lacZ.res$padj),]
noggin_lacZ.resOrdered <- noggin_lacZ.res[order(noggin_lacZ.res$padj),]
fgfr_noggin.resOrdered <- fgfr_noggin.res[order(fgfr_noggin.res$padj),]

fgfr_lacZ.resSig <- subset(fgfr_lacZ.resOrdered, padj < 0.1)
noggin_lacZ.resSig <- subset(noggin_lacZ.resOrdered, padj < 0.1)
fgfr_noggin.resSig <- subset(fgfr_noggin.resOrdered, padj < 0.1)

# output=c("fgfr_lacZ.resOrdered","noggin_lacZ.resOrdered","fgfr_noggin.resOrdered",
#         "fgfr_lacZ.resSig","noggin_lacZ.resSig","fgfr_noggin.resSig")
# 
# for (c in output) {
#   file.name <- paste(c, "_results.csv", sep="")
#   temp <-noquote(c)
#   write.csv(as.data.frame(temp),
#           file=file.name)
# }

write.csv(as.data.frame(fgfr_lacZ.resOrdered),file="15hpf_fgfr_lacZ.resOrdered_results.csv")
write.csv(as.data.frame(noggin_lacZ.resOrdered),file="15hpf_noggin_lacZ.resOrdered_results.csv")
write.csv(as.data.frame(fgfr_noggin.resOrdered),file="15hpf_fgfr_noggin.resOrdered_results.csv")
write.csv(as.data.frame(fgfr_lacZ.resSig),file="15hpf_fgfr_lacZ.resSig_results.csv")
write.csv(as.data.frame(noggin_lacZ.resSig),file="15hpf_noggin_lacZ.resSig_results.csv")
write.csv(as.data.frame(fgfr_noggin.resSig),file="15hpf_fgfr_noggin.resSig_results.csv")

#HTML report for FGFr vs LacZ
des2Report <- HTMLReport(shortName = "15_hpr_fgfr_vs_lacz",
                         title = "DESeq2 RNA-seq analysis of FGFr vs LacZ at 15hpf",
                         reportDirectory = "./reports")
publish(hand_dds,des2Report, pvalueCutoff=0.1, contrast = c("condition","fgfr", "lacz"),
        factor = colData(hand_deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)

#HTML report for FGFr vs Noggin
des2Report <- HTMLReport(shortName = "15_hpr_fgfr_vs_noggin",
                         title = "DESeq2 RNA-seq analysis of FGFr vs Noggin at 15hpf",
                         reportDirectory = "./reports")
publish(hand_dds,des2Report, pvalueCutoff=0.1, contrast = c("condition","fgfr", "noggin"),
        factor = colData(hand_deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)

#HTML report for Noggin vs FGFr
des2Report <- HTMLReport(shortName = "15_hpr_noggin_vs_lacz",
                         title = "DESeq2 RNA-seq analysis of Noggin vs LacZ at 15hpf",
                         reportDirectory = "./reports")
publish(hand_dds,des2Report, pvalueCutoff=0.10, contrast = c("condition","noggin", "lacz"),
        factor = colData(hand_deseq.data)$condition,
        reportDir="./reports")
finish(des2Report)
