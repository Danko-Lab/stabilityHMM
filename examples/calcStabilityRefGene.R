##
## Computes the stability of RefSeq annotations on chr21

## Sanity check.
refGene <- read.table("refGene.hg19.chr21.bed.gz")
hg19 <- "~/storage/data/hg19/hg19.2bit"
seqLen= 1000

seqs <- collectSequences(hg19, refGene, seq.length = seqLen)
mData <- prepareData(seqs)
gene_Scores <- unstableScore(mData)

hist(gene_Scores)

## Get outliers.
cbind(refGene, gene_Scores)[gene_Scores > 0.8,]


