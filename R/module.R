## A function to identify unmodified background using gaussian mixture model
classifyBackground <- function(se, gmm_cut = 5){
  require(mclust)
  IP_count <- assay(se[,se$IP_input == "IP"])
  input_count <- assay(se[,se$IP_input == "input"])

  input_sum <- rowSums(input_count)
  IP_sum <- rowSums(IP_count)
  indx_ratio <- (input_sum >= gmm_cut) & (IP_sum >= gmm_cut)
  logMratio <- log(IP_sum[indx_ratio] / input_sum[indx_ratio])
  rm(input_sum, IP_sum)

  gmm_fit <- Mclust(logMratio, G = 2)
  rm(logMratio)

  bg_class <- gmm_fit$parameters$mean
  bg_indx <- gmm_fit$classification == names(bg_class)[which.min(bg_class)]
  rm(gmm_fit, bg_class)
  rm(IP_count, input_count)

  rowData(se) <- DataFrame(bg = FALSE)
  rowData(se)$bg[which(indx_ratio)[bg_indx]] <- TRUE
  return(se)
}

## A function to estimate sequencing depth size factor from background
estimateColumnFactors <- function(se){
  require(magrittr)
  stopifnot(!is.null(rowData(se)$bg))
  se$sf <- assay(se)[rowData(se)$bg,] %>% apply(., 2, function(x) median(x[x>0]))
  return(se)
}

## A function to estimate feature specific size factors for GC content bias correction
estimateMatrixFactors <- function(se, gc = NULL, gc.knots = seq(from=.4, to=.6, length=3),
                                  gc.bk = c(0,1)){
  require(speedglm)
  require(splines)
  sfm <- matrix(nrow = nrow(se), ncol = ncol(se))
  fitm <- matrix(nrow = 200, ncol = ncol(se))
  fit_range <- data.frame(x=seq(0.2,0.8,length.out = 200))
  if(!is.null(gc)){
    for( i in seq_len(ncol(se)) ){
      count_i <- assay(se)[,i]
      model_Matrix_i <- data.frame(x=gc,y=count_i)
      fit_i <- speedglm(y ~ ns(x, knots=gc.knots, Boundary.knots=gc.bk),
                        family = poisson(link="log"),
                        data = model_Matrix_i,
                        surface = "direct")
      fitted_y <- predict.speedglm(fit_i, newdata = data.frame(x=gc))
      sfm[,i] <- exp(fitted_y)/mean(exp(fitted_y))  # 0 center offsets
      fitted_p <- predict.speedglm(fit_i, newdata = fit_range) #fitted values to generate plot
      fitm[,i] <- exp(fitted_p)/mean(exp(fitted_p))
    }
  }else{
    sfm <- matrix(1, nrow = nrow(se), ncol = ncol(se))
  }
  sfm <- t(t(sfm) * se$sf) #incorporate sequencing depth size factors
  assays(se, withDimnames=FALSE)$sfm <- sfm
  return(se)
}


## A function to count reads overlapped with features
featuresCounts <- function(features,
                           bam_dirs,
                           paired = FALSE,
                           strandness = c("unstrand",
                                          "1st_strand",
                                          "2nd_strand"),
                           parallel = 1,
                           yield_size = 5000000){
    require(GenomicAlignments)
    require(BiocParallel)

    ## Setup parallel number
    register(SerialParam())
    suppressWarnings( register(MulticoreParam(workers = parallel)) )
    register(SnowParam(workers = parallel))

    ## Setup bam file list
    bam_lst = BamFileList(file = bam_dirs,
                             asMates = !paired)
    yieldSize(bam_lst) = yield_size

    ## Count using summarizeOverlaps (equivalent method of Python htseqcount)
    if(strandness == "1st_strand"){
      preprocess_func <- reads_reverse
    }else{
      preprocess_func <- NULL
    }

    se <- summarizeOverlaps(
    features = features,
    reads = bam_lst,
    mode = "Union",
    inter.feature = FALSE,
    singleEnd = !paired,
    preprocess.reads = preprocess_func,
    ignore.strand = strandness == "unstrand",
    fragments = paired
    )

    return(assay(se))
}
