##
## Tunes length and cutoff parameters of the stability model 
## using empirical data in k562 and gm12878.
##

require(stabilityHMM)
require(twoBit)
require(rqhmm)
require(dREG)
require(vioplot)

hg19 <- "/local/storage/data/hg19/hg19.2bit" ## hg19

###
## Read in sequences of known stability.
path<-"/local/storage/projects/validateStabilityHMM/data/"
setwd("/local/storage/projects/validateStabilityHMM")
cell<-"k562"#"gm12878"
stable   <- rbind(read.table(paste(path,"tss_SS_",cell,"_plus.bed",sep="")), 
			read.table(paste(path,"tss_SS_",cell,"_minus.bed",sep="")),
			read.table(paste(path,"tss_US_",cell,"_plus.bed",sep="")),
			read.table(paste(path,"tss_SU_",cell,"_minus.bed",sep="")))

unstable   <- rbind(read.table(paste(path,"tss_UU_",cell,"_plus.bed",sep="")), 
                        read.table(paste(path,"tss_UU_",cell,"_minus.bed",sep="")),
                        read.table(paste(path,"tss_SU_",cell,"_plus.bed",sep="")),
                        read.table(paste(path,"tss_US_",cell,"_minus.bed",sep="")))



getScores <- function(bed, length) {
  seqs <- collectSequences(hg19, bed, seq.length = length)
  mData <- prepareData(seqs)
  unstableScore(mData)
}


###
## First, runs stability over a set of sequence lengths and determines the optimal length.
AUCSScores <- double()
lengths <- c(seq(200, 1000, 100), seq(1500, 2500, 500))
expStab    <- c(rep(0, NROW(stable)), rep(1, NROW(unstable)))
for(length in lengths) {
  StabScores <- c(getScores(stable, length), getScores(unstable, length))
  StabScores[is.na(StabScores)] <- 0
  roc <- logreg.roc.calc(expStab, StabScores)
  AUCSScores <- c(AUCSScores, roc.auc(roc))
  #roc.plot(roc)
}

# Plots AUC v. legnth
data.frame(lengths, AUCSScores)

plot("AUCvLength.pdf")
 plot(lengths, AUCSScores, type="b")
dev.off()
