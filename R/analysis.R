#' generate links for Cytoscape
#' @param L_G_record A list of CRE-gene
#' @param L_TF_record A list of CRE-TF
#' @param P_L_G_record A list of CRE-gene about promoters
#' @param P_L_TF_record A list of CRE-TF about promoters
#' @param groups cell type information
#' @return TF-G dataframe
#' @export

generate_links_for_Cytoscape <- function(L_G_record, L_TF_record, P_L_G_record, P_L_TF_record, groups) {
  TF_G_C <- get_TF_G(L_G_record,L_TF_record, groups)
  P_TF_G_C <- get_TF_G(P_L_G_record,P_L_TF_record, groups)
  # take genes which has both distal TFs and proximal TFs
  intgenes <- intersect(TF_G_C$Gene, P_TF_G_C$Gene)
  TF_G_C <- TF_G_C[which(TF_G_C$Gene %in% intgenes), , drop = FALSE]
  P_TF_G_C <- P_TF_G_C[which(P_TF_G_C$Gene %in% intgenes), , drop = FALSE]
  TF_G_C_final <- rbind(TF_G_C,P_TF_G_C)
  types <- c(rep("distal",nrow(TF_G_C)),rep("proximal",nrow(P_TF_G_C)))
  TF_G_C_final$types <- types
  return(TF_G_C_final)
}

#' generate Node character for Cytoscape
#' @param network_links Links obtained by generate_links_for_Cytoscape
#' @param markers Dataframe with marker and cell type information
#' @return Node character dataframe
#' @export

generate_node_for_Cytoscape <- function(network_links,markers) {
  ### Node
  G <- unique(network_links$Gene)
  T <- unique(network_links$TF)
  Nodes <- c(T,G)
  label_type <- c(rep("TF",length(T)),rep("gene",length(G)))
  unik <- !duplicated(Nodes)
  Nodes <- Nodes[unik]
  label_type <- label_type[unik]

  colnames(markers) <- c("gene","group")
  names <- unique(markers$group)

  ### used to note cell type specific information
  DE_flag <- matrix(0,nrow = length(Nodes), ncol = length(names))
  for (i in 1:length(names)) {
    marker1<- markers[markers$group == names[i],]
    marker_gene <- as.character(marker1$gene)
    id <- which(Nodes %in% marker_gene)
    DE_flag[id,i] <- 1
  }
  group_type <- rep(NA,length(Nodes))
  for (i in 1:length(Nodes)) {
    id <- which(DE_flag[i,] == 1)
    name_i <- names[id[1]]
    if (length(id) == 1) {
      group_type[i] <- name_i
    } else if (length(id) == 2) {
      group_type[i] <- paste0(name_i,"_",names[id[2]])
    } else {
      for (j in 2:(length(id)-1)) {
        name_i <- paste0(name_i,"_",names[id[j]])
      }
      group_type[i] <- paste0(name_i,"_",names[id[length(id)]])
    }
  }
  Node_character <- cbind(label_type,group_type)
  colnames(Node_character) <- c("function_type","group")
  rownames(Node_character) <- Nodes
  return(Node_character)
}

#' Obtain TF and Gene relationships
#' @param L_G_record A list of CRE-gene
#' @param L_TF_record A list of CRE-TF
#' @param groups cluster information
#' @return TF-G dataframe
#'

get_TF_G <- function(L_G_record,L_TF_record,groups) {

  L_TF_G_record <- list()
  TF_G_C_record <- list()
  for (i in 1:length(L_TF_record)) {
    L_T <- L_TF_record[[i]]
    L_G <- L_G_record[[i]]
    genes <- rep(NA,nrow(L_T))
    for (j in 1:nrow(L_T)) {
      id <- which(L_G$loci %in% L_T$loci[j])
      if (length(id) > 0) {
        genes[j] <- L_G$gene[id]
      }
    }
    L_TF_G <- L_T[!is.na(genes), , drop = FALSE]
    L_TF_G$gene <- genes[!is.na(genes)]
    L_TF_G_record[[i]] <- L_TF_G

    TF_G <- L_TF_G[,2:3]
    rownames(TF_G) <- NULL
    #head(TF_G)
    unik <- !duplicated(TF_G)
    TF_G <- TF_G[unik,]
    TF_G_C_record[[i]] <- data.frame(TF = TF_G$TF, Gene = TF_G$gene, Cell_type = groups[i])
  }
  TF_G_C <- do.call(rbind,TF_G_C_record)
  unik <- !duplicated(TF_G_C)
  TF_G_C <- TF_G_C[unik,]
}


#' extract CREs-gene relations of markers
#' @param direct.net_result dataframe of the result of DIRECT-NET
#' @param markers two column dataframe data with marker genes and group information
#' @export

generate_CRE_Gene_links <- function(direct.net_result, markers) {
  # loci_gene corresponding relationship
  direct.net_result_CRE <- direct.net_result[which(direct.net_result$function_type == "HC"), , drop = FALSE]
  L_G_record <- list()
  P_L_G_record <- list() # promoter-gene
  uniqgroups <- unique(markers$group)
  for (i in 1:length(uniqgroups)) {
    marker1<- markers[markers$group == uniqgroups[i],]
    marker_gene <- as.character(marker1$gene)
    effctive_id <- which(direct.net_result_CRE$gene %in% marker_gene)
    L_G_i <- data.frame(loci = direct.net_result_CRE$Peak2[effctive_id], gene = direct.net_result_CRE$gene[effctive_id])
    L_G_record[[i]] <- L_G_i[!duplicated(L_G_i), , drop = FALSE]
    P_L_G_i <- data.frame(loci = direct.net_result_CRE$Peak1[effctive_id], gene = direct.net_result_CRE$gene[effctive_id])
    P_L_G_record[[i]] <- P_L_G_i[!duplicated(P_L_G_i), , drop = FALSE]
  }
  CRE_Gene <- list()
  CRE_Gene$distal <- L_G_record
  CRE_Gene$promoter <- P_L_G_record
  return(CRE_Gene)
}



#' extract CREs of markers
#' @param L_G_record a list of CRE-Gene relationships
#' @param P_L_G_record a list of Promoter-Gene relationships
#' @param da_peaks_list a list of DA of each group
#' @import chromVAR
#' @import motifmatchr
#' @import GenomicRanges
#' @export

generate_CRE <- function(L_G_record, P_L_G_record, da_peaks_list) {

  # extract overlapped peaks between DA and CREs of focused markers

  peaks_bed_record <- list() # peaks used to identify TFs bounded to
  P_peaks_bed_record <- list()

  if (is.null(da_peaks_list)) {
    for (i in 1:length(L_G_record)) {
      da_peaks_list[[i]] <- L_G_record[[i]]$loci
    }
  }
  L_G_record_new <- list()
  P_L_G_record_new <- list()
  for (i in 1:length(L_G_record)) {
    if (!is.null(da_peaks_list[[i]])) {
      if (!is.character(da_peaks_list[[i]])) {
        da <- rownames(da_peaks_list[[i]])
      }
      else {
        da <- da_peaks_list[[i]]
      }
      da <- gsub("-","_",da)
    }

    DA_HC_i <- intersect(da,L_G_record[[i]]$loci)
    if (length(DA_HC_i) > 0) {
      L_G_record_new[[i]] <- L_G_record[[i]][which(L_G_record[[i]]$loci %in% DA_HC_i),]
      overlap_gene <- L_G_record[[i]]$gene[which(L_G_record[[i]]$loci %in% DA_HC_i)]
      peaks1 <- strsplit(DA_HC_i,"_")
      peaks_bed <- matrix(0, nrow = length(DA_HC_i),ncol = 3)
      for (j in 1:length(DA_HC_i)) {
        for (k in 1:3) {
          peaks_bed[j,k] <- peaks1[[j]][k]
        }
      }
      colnames(peaks_bed) <- c("R.chrom","R.start","R.end")
      peaks_bed <- as.data.frame(peaks_bed)
      peaks_bed$R.start <- as.numeric(peaks_bed$R.start)
      peaks_bed$R.end <- as.numeric(peaks_bed$R.end)
      peaks_bed_record[[i]] <- peaks_bed
      ### promoters
      P_L_G_record_i <- P_L_G_record[[i]]
      P_L_G_record_i <- P_L_G_record_i[which(P_L_G_record_i$gene %in% overlap_gene),]
      P_L_G_record_new[[i]] <- P_L_G_record_i
      peaks2 <- strsplit(P_L_G_record_i$loci,"_")
      peaks_bed <- matrix(0, nrow = length(P_L_G_record_i$loci),ncol = 3)
      for (j in 1:length(P_L_G_record_i$loci)) {
        for (k in 1:3) {
          peaks_bed[j,k] <- peaks2[[j]][k]
        }
      }
      colnames(peaks_bed) <- c("R.chrom","R.start","R.end")
      peaks_bed <- as.data.frame(peaks_bed)
      peaks_bed$R.start <- as.numeric(peaks_bed$R.start)
      peaks_bed$R.end <- as.numeric(peaks_bed$R.end)
      P_peaks_bed_record[[i]] <- peaks_bed[!(duplicated(peaks_bed)), , drop = FALSE]
    }
  }
  CREs <- list()
  CREs$distal <- peaks_bed_record
  CREs$promoter <- P_peaks_bed_record
  CREs$L_G_record <- L_G_record_new
  CREs$P_L_G_record <- P_L_G_record_new
  return(CREs)
}

#' Identify TFs enriched in CREs of focus markers
#' @param peaks_bed_list A list of peaks bed file of each group
#' @param species Species used to detect TF
#' @param genome For example: BSgenome.Hsapiens.UCSC.hg19
#' @param markers two column dataframe data with marker genes and group information
#' @import chromVAR
#' @import motifmatchr
#' @import GenomicRanges
#' @export

generate_peak_TF_links <- function(peaks_bed_list, species, genome, markers) {
  motifs <- getJasparMotifs(species)
  L_TF_record <- list()
  for (i in 1:length(peaks_bed_list)) {
    if (!is.null(peaks_bed_list[[i]])) {
      peaks_bed <- peaks_bed_list[[i]]
      peaks_new <- GRanges(seqnames = peaks_bed$R.chrom,
                           ranges = IRanges(start = peaks_bed$R.start, end = peaks_bed$R.end))
      motif_ix <- matchMotifs(motifs, peaks_new, genome,out = "scores")
      S <- as.matrix(motif_ix@assays@data$motifScores)
      M <- as.matrix(motif_ix@assays@data$motifMatches)
      TF <- motif_ix@colData$name

      L_TF_list <- list()
      for (j in 1:nrow(M)) {
        if (sum(M[j,]) > 0) {
          p <- paste0(peaks_bed$R.chrom[j],"_",peaks_bed$R.start[j],"_",peaks_bed$R.end[j])
          # focus TFs in the markers
          TF_j = intersect(unique(markers$gene), TF[M[j,]])
          if (length(TF_j) > 0) {
            L_TF_list[[j]] <- data.frame(loci = p, TF = TF_j)
          }
        }
      }
      L_TF_record[[i]] <- do.call(rbind,L_TF_list)
    }
  }
  return(L_TF_record)
}


