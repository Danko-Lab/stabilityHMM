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
## Second, given a fixed sequence length, chooses a cutoff threshold that optimizes classification accuracy.
expStab    <- c(rep(0, NROW(stable)), rep(1, NROW(unstable)))

length <- 500 #lengths[which.max(AUCSScores)]## BEST.
smax <- 0.0005
umin <- 0.05

SScores <- getScores(stable, length)
UScores <- getScores(unstable, length)

SScores[is.na(SScores)] <- 1
UScores[is.na(UScores)] <- 1

pdf("chooseThreshold.pdf")

# Plot histograms (Fig. 2B).
boxplot(SScores, UScores, names=c("Stable", "Unstable"))
vioplot(SScores, UScores, names=c("Stable", "Unstable"))
vioplot(log(SScores+0.01), log(UScores+0.01), names=c("Stable", "Unstable"), ylab="log(Unstable Score)")

histbreaks <- seq(0,1,0.05)
hist(SScores, breaks=histbreaks, main="Stable")
hist(UScores, breaks=histbreaks, main="Unstable")

plot(ecdf(SScores), col="dark blue", xlim=c(0,1), ylim=c(0,1), xlab="Unstable score")
par(new=TRUE)
plot(ecdf(UScores), col="red", xlim=c(0,1), ylim=c(0,1), xlab="")
abline(v=0.05)
abline(v=0.0005)

# ROC Plot
StabScores <- c(SScores, UScores)
StabScores[is.na(StabScores)] <- 0

roc <- logreg.roc.calc(expStab, StabScores)
roc.auc(roc)
roc.plot(roc)

## Fraction Correct ...
th <- c(seq(0.0,0.05,0.0005), 0.075, 0.1)
co <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE)+sum(UScores>x, na.rm=TRUE))/(NROW(SScores)+NROW(UScores))})
plot(th, co, type="b")

ssen <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE))/(NROW(SScores))})
sppv <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE))/(sum(UScores<x, na.rm=TRUE)+sum(SScores<x, na.rm=TRUE))})

usen <- sapply(th, function(x) {(sum(UScores>x, na.rm=TRUE))/(NROW(UScores))})
uppv <- sapply(th, function(x) {(sum(UScores>x, na.rm=TRUE))/(sum(SScores>x, na.rm=TRUE)+sum(UScores>x, na.rm=TRUE))})

plot(ssen, sppv, xlab="Sensitivity (Stable TU)", ylab="PPV", type="b", col="dark blue")
plot(usen, uppv, xlab="Sensitivity (Unstable TU)", ylab="PPV", type="b", col="red")

data.frame(th, ssen, sppv, usen, uppv, co)

## Accuracy dual threshold...
(sum(SScores < smax)+sum(UScores > umin))/(sum(SScores < smax)+sum(UScores > umin)+sum(SScores > umin)+sum(UScores < smax))
(sum(SScores > smax & SScores < umin)+sum(UScores > umin & UScores < smax))/(NROW(SScores)+NROW(UScores)) ## Fraction unclassified...

dev.off()

###
## Improve AUC by incorporating both a short and a long length??!?
## Strange -- it's actually worse!
score500  <- c(getScores(stable, 500), getScores(unstable, 500))
score1500 <- c(getScores(stable, 1500), getScores(unstable, 1500))
score500[is.na(score500)] <- 1
score1500[is.na(score1500)] <- 1

data_df <- data.frame(y= expStab, score_500= score500, score_1500= score1500)
stabmodel <- glm(y~., family=binomial, data=data_df)#[train,])

pred_stab <- predict(stabmodel, data_df)

roc.auc(logreg.roc.calc(expStab, pred_stab))
roc.auc(logreg.roc.calc(expStab, score500))
roc.auc(logreg.roc.calc(expStab, score1500))

## Strange that this does not give any improvement at all...

## FIN
