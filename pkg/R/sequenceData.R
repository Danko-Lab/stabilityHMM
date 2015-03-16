#

#' Collect DNA sequences for a set of potential transcripts
#'
#' BED file should be such that the 5' end corresponds to the TSS
#'
#' @param twobit.filename path to a 2bit file with sequence data (UCSC 2bit format)
#' @param bed data.frame with regions (UCSC BED6 format)
#' @param seq.length sequence length to use in model (default: 1.5kb)
#'
#' @return list of integer vectors, one per region indicated in the bed parameter (DNA bases encoded as A=1, C=2, G=3, T=4=.
#' @export
collectSequences <- function(twobit.filename, bed, seq.length = 1500) {
  twobit = twoBit::twobit.load(path.expand(twobit.filename))

  N = dim(bed)[1]
  result = vector(mode="list", length = N)

  is.minus = bed[,6] == '-'
  starts = bed[,2]
  ends = bed[,2] + seq.length
  starts[is.minus] = bed[is.minus, 3] - seq.length
  ends[is.minus] = bed[is.minus, 3]

  chroms = as.character(bed[,1])

  for (i in 1:N) {
    chrom = chroms[i]

    seq = twoBit::twobit.sequence(twobit, chrom, starts[i], ends[i])

    if (is.minus[i])
      seq = twoBit::twobit.reverse.complement(seq)

    result[[i]] <- twoBit::twobit.sequence.to.integer(seq)
  }

  return(result)
}

#' Convert sequence data into format used by HMM
#'
#' Convert sequence data obtained from collectSequences() into form used by HMM.
#'
#' @param sequence.lst list of sequences, as returned by collectSequences()
#' @export
prepareData <- function(sequence.lst) {
  res = lapply(sequence.lst, function(seq) c(seq + 1, 1))
  class(res) <- c("shmmData", "list")
  res
}

#' Characterize HMM path probabilities across a collection of sequences
#'
#' Alternative paths represent one of three conditions: (None) neither splice site nor
#' polyA sequence elements are found; (SS5 first) splicing sequence element was found first
#' with or without a polyA sequence element; (pA first) polyA sequence element found first
#' with or without a splice site sequence element.
#'
#' @param data sequence data prepared from applying prepareData() to the result of collectSequences()
#' @param hmm HMM model instance obtained from call to createStabilityHMM()
#'
#' @export
modelPathSummary <- function(data, hmm = createStabilityHMM()) {
  if (!all(class(data) == c("shmmData", "list"))) {
    stop("Unrecognized data class, please use the 'prepareData' and 'collectSequences' functions to create the input dataset.")
  }
  
  # estimate model parameters
  em.res = rqhmm::em.qhmm(hmm, data)

  # get re-normalized outgoing probabilities
  path.probs = rqhmm::get.transition.params.qhmm(hmm, 1)[2:6]
  path.probs = path.probs / sum(path.probs)

  # combine and label paths
  probs = c(
    path.probs[1],
    path.probs[2] + path.probs[3],
    path.probs[4] + path.probs[5])

  names(probs) <- c("None", "SS5 first", "pA first")
  attr(probs, "em.log") <- em.res
  return(probs)
}

#' Predict sequence stability
#'
#' Predict stability based on posterior probability of being in
#' the pA first paths
#'
#' @param data sequence data prepared from applying prepareData() to the result of collectSequences()
#' @param hmm HMM model instance obtained from call to createStabilityHMM()
#' @param n.threads integer number of threads used in parameter inference
#'
#' @return Numeric vector with per sequence probabilities of being an unstable transcript.
#' @export
unstableScore <- function(data, hmm = createStabilityHMM(), n.threads = 1) {
  sapply(data, function(seq) {
    post = rqhmm::posterior.from.state.qhmm(hmm, 1, seq, n_threads = n.threads)

    # outgoing transitions: B, End, SS5.A1, SS5, pA, pA.A2
    out.probs = rowSums(post)[2:6] # drop B
    out.probs = out.probs / sum(out.probs) # renormalize

    # unstable prob := prob of finding a pA first
    sum(out.probs[4:5])
  })
}
