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
call <- paste("zcat refGene.hg19.bed.gz | awk 'BEGIN{OFS=","\t","} {print $1,$6==","+","?$2:$3-1,$6==","+","?$2+1:$3,$4,$5,$6}' | sort-bed - | bedtools closest -s -d -b stdin -a ", sep='"')
cell<-"k562"#"gm12878"
stable   <- rbind( read.table(pipe(paste(call, path,"tss_SS_",cell,"_plus.bed",sep=""))),
		read.table(pipe(paste(call, path,"tss_SS_",cell,"_minus.bed",sep=""))),
		read.table(pipe(paste(call, path,"tss_US_",cell,"_plus.bed",sep=""))),
		read.table(pipe(paste(call, path,"tss_SU_",cell,"_minus.bed",sep=""))))

unstable   <- rbind(read.table(pipe(paste(call, path,"tss_UU_",cell,"_plus.bed",sep=""))), 
                read.table(pipe(paste(call, path,"tss_UU_",cell,"_minus.bed",sep=""))),
                read.table(pipe(paste(call, path,"tss_SU_",cell,"_plus.bed",sep=""))),
		read.table(pipe(paste(call, path,"tss_US_",cell,"_minus.bed",sep=""))))

getScores <- function(bed, length) {
  seqs <- collectSequences(hg19, bed, seq.length = length)
  mData <- prepareData(seqs)
  unstableScore(mData)
}

###
## First, compute distal/ proximal
dist <- 1000
prox <- c(stable$V13 < dist, unstable$V13 < dist)

###
## Second, given a fixed sequence length, chooses a cutoff threshold that optimizes classification accuracy.
expStab    <- c(rep(0, NROW(stable)), rep(1, NROW(unstable)))

length <- 400 #lengths[which.max(AUCSScores)]## BEST.
smax <- 0.1#0.005
umin <- 0.001#0.05

SScores <- getScores(stable, length)
UScores <- getScores(unstable, length)

SScores[is.na(SScores)] <- 1
UScores[is.na(UScores)] <- 1

pdf("chooseThreshold.pdf")
vioplot(log(stable$V13+1, 10), log(unstable$V13+1, 10), names=c("Stable", "Unstable"))
abline(h=log(1000, 10))

# Plot histograms (Fig. 2B).
boxplot(SScores, UScores, names=c("Stable", "Unstable"))
vioplot(SScores, UScores, names=c("Stable", "Unstable"))
vioplot(log(SScores+0.01), log(UScores+0.01), names=c("Stable", "Unstable"))

histbreaks <- seq(0,1,0.05)
hist(SScores, breaks=histbreaks, main="Stable")
hist(UScores, breaks=histbreaks, main="Unstable")

plot(ecdf(SScores), col="dark blue", xlim=c(0,1), ylim=c(0,1), xlab="Unstable score")
par(new=TRUE)
plot(ecdf(UScores), col="red", xlim=c(0,1), ylim=c(0,1), xlab="")
abline(v=0.05)
abline(v=0.0005)

# ROC Plot
index <- sample(NROW(UScores), NROW(SScores))
StabScores <- c(SScores, UScores[index])
StabScores[is.na(StabScores)] <- 0

roc <- logreg.roc.calc(c(rep(0, NROW(stable)), rep(1, NROW(stable))), StabScores)
roc.auc(roc)
roc.plot(roc)

## Fraction Correct ...
th <- c(seq(0.0,0.05,0.001), 0.075, 0.1, 0.25, 0.5, 0.75, 1)
co <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE)+sum(UScores>x, na.rm=TRUE))/(NROW(SScores)+NROW(UScores))})
plot(th, co, type="b")

# Stable
ssen <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE))/(NROW(SScores))})
sppv <- sapply(th, function(x) {(sum(SScores<x, na.rm=TRUE))/(sum(UScores<x, na.rm=TRUE)+sum(SScores<x, na.rm=TRUE))})

ssen_d <- sapply(th, function(x) {(sum(SScores[stable$V13<dist]<x, na.rm=TRUE))/(NROW(SScores))})
sppv_d <- sapply(th, function(x) {(sum(SScores[stable$V13<dist]<x, na.rm=TRUE))/(sum(UScores[unstable$V13<dist]<x, na.rm=TRUE)+sum(SScores[stable$V13<dist]<x, na.rm=TRUE))})

# Unstable
usen <- sapply(th, function(x) {(sum(UScores>x, na.rm=TRUE))/(NROW(UScores))})
uppv <- sapply(th, function(x) {(sum(UScores>x, na.rm=TRUE))/(sum(SScores>x, na.rm=TRUE)+sum(UScores>x, na.rm=TRUE))})

usen_d <- sapply(th, function(x) {(sum(UScores[unstable$V13>dist]>x, na.rm=TRUE))/(NROW(UScores))})
uppv_d <- sapply(th, function(x) {(sum(UScores[unstable$V13>dist]>x, na.rm=TRUE))/(sum(SScores[stable$V13>dist]>x, na.rm=TRUE)+sum(UScores[unstable$V13>dist]>x, na.rm=TRUE))})

plot(ssen, sppv, xlab="Sensitivity (Stable TU)", ylab="PPV", type="b", col="light blue", lwd="dotted", xlim=c(0,1), ylim=c(0.3,1))
points(ssen_d, sppv_d, type="b", col="dark blue")#, xlab="Sensitivity (Stable TU)", ylab="PPV", type="b", col="dark blue")

plot(usen, uppv, type="b", col="pink", lwd="dotted", xlab="Sensitivity (Unstable TU)", ylab="PPV", xlim=c(0,1), ylim=c(0.3,1))
points(usen_d, uppv_d, type="b", col="dark red")#, xlab="Sensitivity (Unstable TU)", ylab="PPV", type="b", col="dark red")

data.frame(th, ssen, sppv, usen, uppv, co)
data.frame(th, ssen_d, sppv_d, usen_d, uppv_d, co)


###
## Improve AUC by incorporating both a short and a long length??!?
## Strange -- it's actually worse!
score500  <- c(getScores(stable, 500), getScores(unstable, 500))
score1500 <- c(getScores(stable, 1500), getScores(unstable, 1500))
score500[is.na(score500)] <- 0
score1500[is.na(score1500)] <- 0

fold_cv <- 0.8
train <- sample(NROW(score500), NROW(score500)*fold_cv)
test <- rep(TRUE, NROW(score500)); test[train] <- FALSE; test <- which(test)

data_df <- data.frame(y= c(rep(1, NROW(stable)), rep(0, NROW(unstable))), score_500= log(score500+0.01), score_1500= log(score1500+0.01), dist=log(1+c(stable$V13, unstable$V13)))
stabmodel <- glm(y~., family=binomial, data=data_df[train,])

pred_stab <- predict(stabmodel, data_df[test,])

roc.auc(logreg.roc.calc(data_df$y[test], pred_stab))
roc.auc(logreg.roc.calc(data_df$y[test], c(stable$V13, unstable$V13)[test]))
roc.auc(logreg.roc.calc(data_df$y[test], score500[test]))
roc.auc(logreg.roc.calc(data_df$y[test], score1500[test]))

roc.plot(logreg.roc.calc(data_df$y[test], pred_stab))

dev.off()

## FIN
