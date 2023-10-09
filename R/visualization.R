#'@title Plot transcript topology distribution for peaks
#'
#'@return save PDF files of tx topology
#'@export
plotTopology <- function(list_grl,
                         txdb,
                         savePrefix="topology"){
  require(ggplot2)
  topList <- lapply(list_grl,
                    function(x) topologyOnTranscripts(x,txdb))
  plot_df <- data.frame(topology = unlist(topList),
                        group = rep(names(topList), elementNROWS(topList)))
  ggplot(plot_df) + geom_density(aes(x = topology, color = group),
                                 data = plot_df) +
    theme_bw() + scale_color_brewer(palette = "Set2") +
    geom_vline(xintercept = c(0.33, 0.66), linetype = 2) +
    geom_text(aes(x=x,y=y,label=text),
              data = data.frame(
                x = c(0.165, 0.495, 0.825),
                y = c(-0.2, -0.2, -0.2),
                text = c("5'UTR","CDS", "3'UTR")
              )) + scale_x_continuous(breaks = c(0, 0.33, 0.66, 0.9))
  ggsave(paste0(savePrefix, ".pdf"), width = 6, height = 3)
}

#'@title Plot Poisson GLM fits for GC content bias detection
#'
#'@return save PDF files of GC bias curves
#'
plotGCbias <- function(se,
                       savePrefix = "gc_bias"){
  require(ggplot2)
  plot_df <- data.frame(glmFit = as.numeric(metadata(se)[["fitm"]]),
                        sample = rep(paste0("sample_", 1:ncol(se)), each = 200),
                        IP_input = rep(se$IP_input, each = 200),
                        gc = rep(seq(0.2,0.8,length.out = 200), ncol(se)))

  if(!is.null(se$Perturbation)) {
    plot_df$Perturbation <- rep(se$Perturbation, each = 200)
    ggplot(plot_df, aes(x=gc, y=glmFit, group = sample, colour = IP_input, linetype = Perturbation)) +
      geom_line() + scale_color_manual(values = c('#44AA99', '#332288')) +
      theme_bw() + labs(x = "GC content", y="Normalized Poisson GLM Fit")
  }else{
    ggplot(plot_df, aes(x=gc, y=glmFit, group = sample, colour = IP_input)) +
      geom_line() + scale_color_manual(values = c('#44AA99', '#332288')) +
      theme_bw() + labs(x = "GC content", y="Normalized Poisson GLM Fit")
  }

  ggsave(paste0(savePrefix,".pdf"), width = 5, height = 3)
}

