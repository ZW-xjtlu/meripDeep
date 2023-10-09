#'@title Differential peak calling with MeRIP-Seq.
#'
#'@return the peaks in GRangesList format.
#'@export
diffAnalysis <- function(bam_IP,
                         bam_input,
                         bam_IP_treated,
                         bam_input_treated,
                         txdb = NULL,
                         genome = NULL,
                         bin_size = 25,
                         step_size = 25,
                         fragment_length = 100,
                         paired_end = FALSE,
                         strandness = c("unstrand", "1st_strand", "2nd_strand"),
                         gff = NULL,
                         p_cutoff = 1e-3,
                         lfc_threshold = 0,
                         alt_hypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
                         poisson = TRUE,
                         save_se = FALSE,
                         plot_gc = FALSE,
                         parallel = 1){
    require(GenomicRanges)
    require(GenomicFeatures)
    require(SummarizedExperiment)
    require(DESeq2)

    strandness <- match.arg(strandness)
    alt_hypothesis <- match.arg(alt_hypothesis)

    if(is.null(gff) & is.null(txdb)){
      stop("Please at least provide one among gff and TxDb for transcript annotation.")
    }

    if(is.null(txdb) & !is.null(gff)){
      txdb <- makeTxDbFromGFF(gff)
    }

    if(!(all(file.exists(bam_IP)) &
         all(file.exists(bam_input)) &
         all(file.exists(bam_IP_treated)) &
         all(file.exists(bam_input_treated)))){
      stop("At least one bam file directories provided cannot be found.")
    }

    #Extract bins for count
    exByGene  <- exonsByiGenes(txdb)
    peakBins <- exonicBins(exByGene, bin_size, step_size)
    mcols(peakBins) <- NULL
    fragmentBins <- suppressWarnings(exonicFlank(peakBins,
                                              exByGene,
                                              fragment_length - floor(bin_size/2)))
    peakBins <- peakBins[names(fragmentBins)]

    #Count the bam files
    bam_dirs <- c(bam_IP, bam_IP_treated, bam_input, bam_input_treated)
    asy <- featuresCounts(peakBins, bam_dirs, paired_end, strandness, parallel)
    rm(bam_dirs)

    #Create SummarizedExperiment
    se <- SummarizedExperiment(assays = list(counts = asy), rowRanges = peakBins)
    length_indx <- c(length(bam_IP), length(bam_IP_treated),
                     length(bam_input), length(bam_input_treated))
    se$IP_input <- rep(c("IP", "IP", "input", "input"), length_indx)
    se$Perturbation <- rep(c("C", "Treated", "C", "Treated"), length_indx)
    rm(asy, peakBins, length_indx)

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
    se$Perturbation <- factor(se$Perturbation) %>% relevel(.,"C") #Set reference level
    dds <- DESeqDataSet(se, ~ IP_input * Perturbation)
    normalizationFactors(dds) <- assays(se)[["sfm"]]

    if(poisson){
      dispersions(dds) <- 0
    }else{
      dds <- estimateDispersions(dds)
    }

    res <- nbinomWaldTest(dds) %>% results(., altHypothesis = alt_hypothesis,
                                              lfcThreshold = lfc_threshold)
    pvalue <- res$pvalue
    pvalue[rowData(se)$bg] <- 1 #Restrict differential modification on GMM foreground
    rm(res, dds)

    ## Plot GC bias fits
    if(plot_gc) plotGCbias(se)

    #Export intermediate data for differential analysis
    if(save_se) saveRDS(se, "se.rds")

    #Return reduced peaks
    if(length(p_cutoff) > 1){
      Peaklist <- lapply(p_cutoff, function(x) {
        peak <- reducePeaks(rowRanges(se)[which(pvalue < x)], exByGene)
        return(peak[sum(width(peak)) >= bin_size])
      })
      names(Peaklist) <- p_cutoff
      return(Peaklist)
    }else{
      peak <- reducePeaks(rowRanges(se)[which(pvalue < p_cutoff)], exByGene)
      return(peak[sum(width(peak)) >= bin_size])
  }
}
