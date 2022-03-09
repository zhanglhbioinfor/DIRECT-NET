##' Plot links of DIRECT-NET on focused genes
#'
#' @param DIRECT_NET_Result DIRECT-NET result
#' @param gene_anno the genome annotation
#' @param marker the focused gene
#' @param cutoff the cutoff of importance scores, CREs are links with importance scores higher than cutoff
#' @param upstream the number of distance upstream to TSS
#' @param downstream the number of distance downstream to TSS
#' @importFrom cicero plot_connections
#' @return
#' @export
Plot_connections <- function(DIRECT_NET_Result, gene_anno, marker, cutoff = NULL, upstream = 250000, downstream = 250000){
  # take out the result of marker
  if (length(which(DIRECT_NET_Result$gene %in% marker)) > 0) {
    DIRECT_NET_Result <- DIRECT_NET_Result[which(DIRECT_NET_Result$gene %in% marker), ,drop = FALSE]
    conns <- DIRECT_NET_Result[,5:7]
    names(conns) <- c("Peak1","Peak2","coaccess")
    rownames(conns) <- NULL

    # normalize
    conns$coaccess <- log10(as.numeric(conns$coaccess)*100)
    if (is.null(cutoff)) {
      cutoff = max(0,max(conns$coaccess[which(DIRECT_NET_Result$type == "MC")]))
    }
    plot_connections(conns, DIRECT_NET_Result$Chr[1], as.numeric(DIRECT_NET_Result$Starts[1])-upstream, as.numeric(DIRECT_NET_Result$Starts[1])+downstream,
                     gene_model = gene_anno,
                     coaccess_cutoff = cutoff,
                     connection_width = .5,
                     connection_color = "#99000d",
                     peak_color = "#666666",
                     collapseTranscripts = "longest" )
  } else {
    message("Please try other markers!")
  }
}


