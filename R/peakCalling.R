#'@title Peak calling with MeRIP-Seq.
#'
#'@return the peaks in GRangesList format.
#'@export
peakCalling <- function(bam_IP,
                        bam_input,
                        txdb = NULL,
                        genome = NULL,
                        bin_size = 25,
                        step_size = 25,
                        fragment_length = 100,
                        paired_end = FALSE,
                        strandness = c("unstrand", "1st_strand", "2nd_strand"),
                        gff = NULL,
                        p_cutoff = 1e-10,
                        poisson = TRUE,
                        plot_gc = FALSE,
                        parallel = 1){
  require(GenomicRanges)
  require(SummarizedExperiment)
  require(DESeq2)

  strandness <- match.arg(strandness)

  if(is.null(gff) & is.null(txdb)){
    stop("Please at least provide one among gff and TxDb for transcript annotation.")
  }

  if(is.null(txdb) & !is.null(gff)){
    txdb <- makeTxDbFromGFF(gff)
  }

  if(!(all(file.exists(bam_IP)) & all(file.exists(bam_input)))){
    stop("At least one bam file directories provided cannot be found.")
  }

  #Extract bins for count
  exByGene  <- exonsByiGenes(txdb)
  peakBins <- exonicBins(exByGene, bin_size, step_size)
  mcols(peakBins) <- NULL
  fragmentBins <- suppressWarnings(exonicFlank(peakBins, exByGene, fragment_length - floor(bin_size/2)))
  peakBins <- peakBins[names(fragmentBins)]

  #Count the bam files
  bam_dirs <- c(bam_IP, bam_input)
  asy <- featuresCounts(peakBins, bam_dirs, paired_end, strandness, parallel)
  rm(bam_dirs)

  #Create SummarizedExperiment
  se <- SummarizedExperiment(assays = list(counts = asy), rowRanges = peakBins)
  se$IP_input <- rep(c("IP", "input"), c(length(bam_IP), length(bam_input)))
  rm(asy, peakBins)

  #Identify Backgrounds
  se <- classifyBackground(se)

  #Estimate sample size factors
  se <- estimateColumnFactors(se)

  #Calculate GC contents
  gc <- NULL
  if(!is.null(genome)){
    gcFreq <- letterFrequency(DNAStringSet(Views(genome, unlist(fragmentBins))), letters = "GC")
    gcFreq_gr <- tapply(gcFreq, rep(1:length(fragmentBins), elementNROWS(fragmentBins)), sum)
    gc <- gcFreq_gr/sum(width(fragmentBins))
    rm(gcFreq, gcFreq_gr, fragmentBins)
  }

  #Fit GC content biases
  se <- estimateMatrixFactors(se, gc)

  #DESeq2 peak calling
  se$IP_input <- factor(se$IP_input) %>% relevel(.,"input") #Set reference level
  dds <- DESeqDataSet(se, ~ IP_input)
  normalizationFactors(dds) <- assays(se)[["sfm"]]

  if(poisson){
    dispersions(dds) <- 0
  }else{
    dds <- estimateDispersions(dds)
  }

  dds = nbinomWaldTest(dds)
  res <- results(dds, altHypothesis = "greater")
  pvalue <- res$pvalue
  rm(res, dds)

  ## Plot GC bias fits
  if(plot_gc) plotGCbias(se)

  #Return granges
  pvalue[is.na(pvalue)] <- 0
  gr_return <- rowRanges(se)
  mcols(gr_return)$IP_cov <- rowSums( assay(se)[,se$IP_input == "IP"] )
  mcols(gr_return)$input_cov <- rowSums( assay(se)[,se$IP_input == "input"] )
  mcols(gr_return)$pvalue <- -1*log(pvalue)
  mcols(gr_return)$peakCalling_result <- as.numeric(pvalue < p_cutoff)
  fol <- findOverlaps(gr_return, exByGene)
  mcols(gr_return)$Gene_indx <- NA
  mcols(gr_return)$Gene_indx[queryHits(fol)] <- names(exByGene)[subjectHits(fol)]
  mcols(gr_return)$bg <- NULL
  return(gr_return)
}
