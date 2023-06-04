#' Create aggregated data for DIRECT-NET
#'
#' Function to generate aggregated inputs for DIRECT-NET. \code{Aggregate_data}
#'
#' @param object Seurat object.
#' @param k_neigh Number of cells to be aggregated per group.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param reduction.name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, should accessibility values be normalized by size factor
#' @param verbose Logical, should warning and info messages be printed?
#' @param seed Random seed
#' @import Seurat
#' @return Aggregated Seurat object.
#' @export
#'
Aggregate_data <- function (object, k_neigh = 50, atacbinary = TRUE, max_overlap=0.8, reduction.name = NULL,
                            size_factor_normalize = TRUE, seed = 123, verbose = TRUE)
{
  if (!is.null(reduction.name)) {
      cell_coord <- object@reductions[[reduction.name]]
  } else {
      if ("RNA" %in% names(object@assays)) {
        cell_coord <- object@reductions$wnn.umap@cell.embeddings
      } else {
        cell_coord <- object@reductions$umap@cell.embeddings
      }
  }
  

  group <- as.character(Idents(object))
  uniqgroup <- unique(group)
  if ("RNA" %in% names(object@assays)) {
    rna_new <- matrix(0,nrow = nrow(object@assays$RNA@counts), ncol =1)
  }
  atac_new <- matrix(0,nrow = nrow(object@assays$ATAC@counts), ncol =1)
  cell_sample <- matrix(0,nrow = 1,ncol = k_neigh)
  for (i in 1:length(uniqgroup)) {
    if (verbose)  {
      message(paste0("Aggregating cluster ",uniqgroup[i]))
    }
    subobject <- subset(object, idents = uniqgroup[i])
    sub_index <- which(group %in% uniqgroup[i])
    cell_coord_i <- cell_coord[sub_index,]
    sub_aggregated_data <- generate_aggregated_data(subobject, cell_coord_i, k_neigh, atacbinary, max_overlap, seed, verbose)

    sub_cell_sample <- sub_aggregated_data$cell_sample
    if ("RNA" %in% names(object@assays)) {
      rna_new <- cbind(rna_new,sub_aggregated_data$rna)
    }
    atac_new <- cbind(atac_new,sub_aggregated_data$atac)
    if (ncol(sub_cell_sample) < k_neigh) {
      sub_cell_sample_new <- as.matrix(sub_cell_sample)
      sub_cell_sample_new <- cbind(sub_cell_sample_new,matrix(0,nrow = 1,ncol = k_neigh - ncol(sub_cell_sample_new)))
    } else {
      sub_cell_sample_new <- apply(sub_cell_sample, 2, function(x) {
        sub_index[x]#for each column return original index
      })
      sub_cell_sample_new <- as.data.frame(sub_cell_sample_new)
      sub_cell_sample_new <- as.matrix(sub_cell_sample_new)
    }
    cell_sample <- rbind(cell_sample,sub_cell_sample_new)
  }
  if ("RNA" %in% names(object@assays)) {
    rna_new <- rna_new[,-1]
  }
  atac_new <- atac_new[,-1]
  cell_sample <- cell_sample[-1,]

  ######### normalization

  if (size_factor_normalize) {
    if ("RNA" %in% names(object@assays)) {
      rna_new <- t(t(log(rna_new+1))/estimateSizeFactorsForMatrix(rna_new))
    }
    atac_new <- t(t(log(atac_new+1))/estimateSizeFactorsForMatrix(atac_new))
  }
  new_data <- list()
  if ("RNA" %in% names(object@assays)) {
    new_data$rna <- rna_new
  }
  new_data$atac <- atac_new
  new_data$cell_sample <- cell_sample
  return (new_data)
}


#' Create aggregated data for a certain cluster
#'
#' Function to generate aggregated inputs of a cetrain cluster. \code{generate_aggregated_data}
#' takes as input sparse data. This function will aggregate binary accessibility scores (or gene expression)
#' per cell cluster, if they do not overlap any existing group with more than 50% cells.
#'
#' @param object Seurat object.
#' @param cell_coord similarity matrix or dimiension reductions.
#' @param k_neigh Number of cells to aggregate per group.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param seed Random seed
#' @param verbose Logical, should warning and info messages be printed?
#' @importFrom FNN knn.index
#' @import Matrix
#'
#' @return Aggregated data.
#' @export
#'
generate_aggregated_data <- function (object, cell_coord, k_neigh = 50, atacbinary = TRUE, max_overlap=0.8, seed = 123, verbose = TRUE)
{
  if (nrow(cell_coord) > k_neigh) {
    # Create a k-nearest neighbors map
    nn_map <- as.data.frame(FNN::knn.index(cell_coord,
                                           k = (k_neigh - 1)))
    row.names(nn_map) <- row.names(cell_coord)
    nn_map$agg_cell <- 1:nrow(nn_map)
    good_choices <- 1:nrow(nn_map)

    if (verbose)
      message("Sample cells randomly.")

    #Sample cells randomly
    set.seed(seed)
    choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]

    it <- 0
    ##Slow (contain calculating of overlapping between cell groups)
    while (length(good_choices) > 0 & it < nrow(cell_coord)/((1-max_overlap)*k_neigh)) {
      it <- it + 1
      choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]

      #calculate overlapping between cell groups
      combs <- data.frame(1:(nrow(cell_sample) - 1), nrow(cell_sample))
      shared <- apply(combs, 1, function(x) {    #Slow
        (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x,
                                                                      ]))))
      })
      if (max(shared) < max_overlap * k_neigh) {
        chosen <- new_chosen
      }
    }

    #aggregating both scRNA-seq and scATAC-seq counts of cells within one group
    if ("RNA" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$RNA@counts)
      rna_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(rna_old)) %in% cell_sample[x,,drop=FALSE])
      rna_mask <- Matrix::Matrix(rna_mask)
      rna_new <- rna_old %*% rna_mask
      rna_new <- as.matrix(rna_new)
    }

    atac_old <- object@assays$ATAC@counts
    # binarize
    if (atacbinary) {
      atac_old <- atac_old > 0
    }

    atac_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(atac_old)) %in% cell_sample[x,,drop=FALSE])
    atac_mask <- Matrix::Matrix(atac_mask)
    atac_new <- atac_old %*% atac_mask
    atac_new <- as.matrix(atac_new)

  } else {
    if ("RNA" %in% names(object@assays)) {
      rna_old <- as.matrix(object@assays$RNA@counts)
      rna_new <- rowSums(rna_old)
      rna_new <- as.matrix(rna_new) }

    atac_old <- object@assays$ATAC@counts
    # binarize
    if (atacbinary) {
      atac_old <- atac_old > 0
    }
    atac_new <-rowSums(atac_old)
    atac_new <- as.matrix(atac_new)
    cell_sample <- as.data.frame(t(matrix(seq(from = 1, to = nrow(cell_coord)))))
  }
  new_data <- list()
  if ("RNA" %in% names(object@assays)) {
    new_data$rna <- rna_new
  }
  new_data$atac <- atac_new
  new_data$cell_sample <- cell_sample
  return (new_data)
}


#' Function to calculate the size factor for the single-cell data
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#' @importFrom stats median
#' @import slam
#' @import Matrix
#' @export
estimateSizeFactorsForMatrix <- function(counts, locfunc = median, round_exprs=TRUE,  method="mean-geometric-mean-total")
{
  #library("slam")
  if (isSparseMatrix(counts)){
    estimateSizeFactorsForSparseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs, method=method)
  }else{
    estimateSizeFactorsForDenseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs,  method=method)
  }

}

#' Convert a slam matrix to a sparseMatrix
#' @param simpleTripletMatrix A slam matrix
#' @import slam
#' @export
asSparseMatrix = function (simpleTripletMatrix) {
  retVal = sparseMatrix(i=simpleTripletMatrix[["i"]],
                        j=simpleTripletMatrix[["j"]],
                        x=simpleTripletMatrix[["v"]],
                        dims=c(simpleTripletMatrix[["nrow"]],
                               simpleTripletMatrix[["ncol"]]))
  if (!is.null(simpleTripletMatrix[["dimnames"]]))
    dimnames(retVal) = simpleTripletMatrix[["dimnames"]]
  return(retVal)
}

#' Convert a sparseMatrix from Matrix package to a slam matrix
#' @param sp_mat The matrix for the aggregated single cell data
#' @import slam
#' @export
asSlamMatrix = function (sp_mat) {
  sp <- Matrix::summary(sp_mat)
  simple_triplet_matrix(sp[,"i"], sp[,"j"], sp[,"x"], ncol=ncol(sp_mat), nrow=nrow(sp_mat), dimnames=dimnames(sp_mat))
}

#' Convert a sparseMatrix from Matrix package to a slam matrix
#' @param x The sparseMatrix data
#' @import Matrix
#' @export
isSparseMatrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}


#' Estimate size factors for each column, given a sparseMatrix from the Matrix package
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
#' @import slam
#' @importFrom stats median
#' @export
estimateSizeFactorsForSparseMatrix <- function(counts,
                                               locfunc = median,
                                               round_exprs=TRUE,
                                               method="mean-geometric-mean-total"){
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  CM <- asSlamMatrix(CM)

  if (method == "weighted-median"){

    log_medians <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowapply_simple_triplet_matrix(CM, function(x) { mean(log(CM)) })

    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    stop("Error: method 'median' not yet supported for sparse matrices")
  }else if(method == 'mode'){
    stop("Error: method 'mode' not yet supported for sparse matrices")
  }else if(method == 'geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

#' Estimate size factors dense matrix
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#' @importFrom stats median
#' @export
estimateSizeFactorsForDenseMatrix <- function(counts, locfunc = median, round_exprs=TRUE, method="mean-geometric-mean-total"){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}


#' Find the most commonly occuring relative expression value in each cell
#'
#' Converting relative expression values to mRNA copies per cell requires
#' knowing the most commonly occuring relative expression value in each cell
#' This value typically corresponds to an RPC value of 1. This function
#' finds the most commonly occuring (log-transformed) relative expression value
#' for each column in the provided expression matrix.
#'
#' @param relative_expr_matrix a matrix of relative expression values for
#' values with each row and column representing genes/isoforms and cells,
#' respectively. Row and column names should be included.
#' Expression values should not be log-transformed.
#' @param relative_expr_thresh Relative expression values below this threshold
#' are considered zero.
#' @return an vector of most abundant relative_expr value corresponding to the
#' RPC 1.
#' @details This function estimates the most abundant relative expression value
#' (t^*) using a gaussian kernel density function. It can also optionally
#' output the t^* based on a two gaussian mixture model
#' based on the smsn.mixture from mixsmsn package
#' @export
#' @examples
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' t_estimate = estimate_t(HSMM_fpkm_matrix)
#'}

estimate_t <- function(relative_expr_matrix, relative_expr_thresh = 0.1) {
  #apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr) 10^mean(dmode(log10(relative_expr[relative_expr > relative_expr_thresh]))))) #avoid multiple output
}


#' use gaussian kernel to calculate the mode of transcript counts
#' @param x log tranformed relative expression
#' @param  breaks control parameter
#' @importFrom stats density
dmode <- function(x, breaks="Sturges") {
  if (length(x) < 2) return (0);
  den <- stats::density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )
}


#' run DIRECT-NET on Seurat object
#'
#' @param object Seurat object.
#' @param peakcalling  call peak
#' @param macs2.path  path to macs2
#' @param fragments  fragments file
#' @param k_neigh Number of cells to aggregate per group.
#' @param atacbinary Logical, should accessibility values be binarized
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param reduction.name The reduction name of extracting the cell coordinates used for aggregating.
#' @param size_factor_normalize Logical, whether need to do size normalization
#' @param genome.info the TSS information of genome, e.g. hg19, hg38
#' @param focus_markers the focused genes
#' @param params the list of parameters used in Xgboost
#' @param nthread  the number of threads can be manually specified in Xgboost trainning stage, default is 2
#' @param early_stop Logical, whether use early stop rule on validation data to reduce overfitting
#' @param HC_cutoff the threshold of high functional CREs
#' @param LC_cutoff the threshold of low functional CREs
#' @param rescued Logical, whether to rescue highly correlated CREs
#' @param seed Random seed
#' @param verbose Logical, should warning and info messages be printed
#' @import Signac
#' @import Seurat
#' @import Matrix
#' @import xgboost
#' @import ggplot2
#' @importFrom  stats quantile
#' @importFrom cicero find_overlapping_coordinates
#' @return a Seurat object with new links assay.
#' @export
Run_DIRECT_NET <- function(object,peakcalling = FALSE, macs2.path = NULL, fragments = NULL, k_neigh = 50, atacbinary = TRUE, max_overlap=0.8, reduction.name = NULL, size_factor_normalize = FALSE, genome.info, focus_markers, params = NULL, nthread = 2, early_stop = FALSE, HC_cutoff = NULL, LC_cutoff = NULL, rescued = FALSE,seed = 123, verbose = TRUE) {
    ########################################################### step 0. Peak calling
    
    if(peakcalling) {
        if (verbose)  {
          message("Calling Peak")
        }
        object$cluster <- Idents(object)
        if(is.null(macs2.path)) {
          message("Please give the path to macs2!")
        }
        peaks <- CallPeaks(
          object = object,
          group.by = "cluster",
          macs2.path = macs2.path
        )
        if(is.null(fragments)) {
          message("Please input fragments!")
        }
        new_atac_data <- FeatureMatrix(
          fragments = fragments, # fragments of original, we need to do aggregation
          features = peaks
        )
        object@assays$ATAC@counts <- new_atac_data
        if (verbose)  {
          message("Peak calling finished")
        }
    }
    
  ########################################################### step 1. Aggregation
  if (verbose)  {
    message("Generating aggregated data")
  }
  if ("aggregated.data" %in% names(Misc(object))) {
    agg.data <- Misc(object, slot = 'aggregated.data')
  } else {
    agg.data <- Aggregate_data(object,k_neigh = k_neigh, atacbinary = atacbinary, max_overlap=max_overlap, reduction.name = NULL,size_factor_normalize = size_factor_normalize)
    Misc(object, slot = 'aggregated.data') <- agg.data
  }


  ########################################################### step2.  Run_DIRECT_NET
  #library(cicero)
  #library(Matrix)
  #library(data.table)
  ## import gbm package
  #library(xgboost)      # a faster implementation of gbm
  #library(ggplot2)      # model visualization
  options(stringsAsFactors = FALSE)

  ########## Implement gbm for each gene  ##########
  # parameter list
  if (is.null(params)) {
    params <- list(
      eta = 0.3,
      max_depth = 6,
      min_child_weight = 1,
      subsample = 1,
      colsample_bytree = 1,
      lambda = 1
    )
  }

  if ("rna" %in% names(agg.data)) {
    data_rna <- as.matrix(agg.data$rna)
    rna <- rownames(data_rna)
    rna <- lapply(rna, function(x) strsplit(x,"[.]")[[1]][1])
    rna <- unlist(rna)
    #rna <-toupper(rna)
    rownames(data_rna) <- rna
    unik <- !duplicated(rna)# filter out different transcript
    data_rna <- data_rna[unik,]
  }

  data_atac <- as.matrix(agg.data$atac)
  rownames(data_atac) <- gsub("-", "_", rownames(data_atac))
  peaks <- rownames(data_atac)

  ########## Obtain candidate regions of focus markers ##########
  genes <- lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
  genes <- lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
  genes <- unlist(genes)
  genome.info$genes <- genes
  unik <- !duplicated(genes)
  genome.info <- genome.info[unik,]

  focus_markers <- lapply(focus_markers, function(x) strsplit(x,"[.]")[[1]][1])
  focus_markers <- unique(unlist(focus_markers))
  focus_markers <- genome.info$genes[which(genome.info$genes %in% focus_markers)]

  genome.info.used <- genome.info[which(genome.info$genes %in% focus_markers),]
  Chr <- genome.info.used$Chrom
  Starts <- genome.info.used$Starts
  Ends <- genome.info.used$Ends

  DIRECT_NET_Result <- list()
  TXs <- list()
  TYs <- list()
  for (i in 1:length(focus_markers))
  {
    if (verbose)  {
      message(paste0("Inferring links for ",focus_markers[i]))
    }

    p1 <- paste(Chr[i],":",Starts[i]-500,"-",Starts[i],sep = "")
    p2 <- paste(Chr[i],":",Starts[i]-250000,"-",Starts[i]+250000,sep = "")
    promoters <- find_overlapping_coordinates(peaks, p1)
    enhancers <- find_overlapping_coordinates(peaks, p2)
    enhancers <- setdiff(enhancers,promoters)

    ########## build data matrix of each gene used in gbm model ##########
    if ("rna" %in% names(agg.data)) {
      idx <- which(rownames(data_rna) == focus_markers[i])
    } else {
      idx <- 1
    }

    if ((length(promoters) > 0 && length(enhancers) > 1) && length(idx) != 0){
      id1 <- match(promoters,peaks)
      id1 <- id1[!is.na(id1)]
      id2 <- match(enhancers,peaks)
      id2 <- id2[!is.na(id2)]
      id2_new <- setdiff(id2,id1)
      X <- data_atac[id2_new,]
      TXs[[i]] <- X
      Y <- data_atac[id1,]
      if (length(id1) > 1){
        Y <- colSums(Y)
      }
      Y <- t(as.matrix(Y))
      rownames(Y) <- peaks[id1[1]]
      TYs[[i]] <- Y
      if ("rna" %in% names(agg.data)) {
        Z <- data_rna[idx,]
        Z <- t(as.matrix(Z))
        rownames(Z) <- focus_markers[i]
      } else {
        Z <- Y
      }
      flag <- 1
    } else {
      flag <- 0
      message(paste0("There are less than two peaks detected within 500 kb for ",focus_markers[i]))
    }

    if (flag == 1) {
      X <- as.matrix(X)
      if (ncol(X) == 1) {
        X <- t(X)
      }
      Y <- as.matrix(Y)
      Z <- as.matrix(Z)
      rownames(Z) <- rownames(Y)
      ### whether use early stop rule
      if (early_stop) {
        if (length(size_factor)/5 < 100) {
          message("The number of cells is too small to split!")
        }
        set.seed(seed)
        cv_idx <- sample(1:5, size=length(size_factor), replace=T)
        test_idx <- which(cv_idx == 1)
        validation_idx <- which(cv_idx == 2)
        x_train <- as.matrix(X[,-c(test_idx,validation_idx)])
        x_test <- as.matrix(X[,test_idx])
        x_validation <- as.matrix(X[,validation_idx])
        y_train <- Y[-c(test_idx,validation_idx)]
        y_test <- Y[test_idx]
        y_validation <- Y[validation_idx]
        
        dtrain = xgb.DMatrix(data=t(x_train), label = as.numeric(y_train))
        dtest = xgb.DMatrix(data=t(x_test), label= as.numeric(y_test))
        dvalidation = xgb.DMatrix(data=t(x_validation), label= as.numeric(y_validation))
        watchlist1 = list(train=dtrain, test=dvalidation)
        watchlist2 = list(train=dtrain, test=dtest)
        xgb_v <- xgb.train(params = params,
                           data = dtrain,
                           watchlist = watchlist1,
                           nrounds = 100,
                           nthread = nthread,
                           objective = "reg:linear",
                           verbose = 0)
        cv1 <- xgb_v$evaluation_log
        rmse_d <- cv1$test_rmse-cv1$train_rmse
        rmse_dd <- abs(rmse_d[2:100]-rmse_d[1:99])/rmse_d[1]
        stop_index <- which(rmse_dd == min(rmse_dd))
        # train final model
        xgb.fit.final <- xgboost(
          params = params,
          data = t(X),
          label = as.numeric(Z),
          nrounds = stop_index[1],
          nthread = nthread,
          objective = "reg:squarederror",
          verbose = 0
        ) 
        } else {
        # train final model
        xgb.fit.final <- xgboost(
          params = params,
          data = t(X),
          label = as.numeric(Z),
          nrounds = 100,
          nthread = nthread,
          objective = "reg:squarederror",
          verbose = 0
        )
      }
     
      # create importance matrix
      tryCatch({
        importance_matrix <- xgb.importance(model = xgb.fit.final)
        Imp_peak <- importance_matrix$Feature
        Imp_peak <- as.vector(Imp_peak)
        Imp_value <- importance_matrix$Gain
        Imp_peak_h <- Imp_peak
        Imp_value_h <- Imp_value
        conns_h <- list()
        conns_h$Peak1 <- as.character(rownames(Y))
        conns_h$Peak2 <-  as.character(Imp_peak_h)
        conns_h$Importance <- Imp_value_h
        conns_h <- as.data.frame(conns_h)
        colnames(conns_h) <- c("Peak1","Peak2","Importance")

      }, error = function(e){
      })
      if (length(ncol(conns_h))){
        DIRECT_NET_Result[[i]] <- conns_h
     }
    }
  }
  conns <- do.call(rbind,DIRECT_NET_Result)
  ######################### variable selection ###########
  if (is.null(HC_cutoff)) {
    HC_cutoff = max(stats::quantile(conns$Importance,0.50),0.001)
  }
  if (is.null(LC_cutoff)) {
    LC_cutoff = min(0.001,stats::quantile(conns$Importance,0.25))
  }
  for (i in 1:length(DIRECT_NET_Result)) {
    if (!is.null(DIRECT_NET_Result[[i]])) {
      conns_h <- DIRECT_NET_Result[[i]]
      Imp_value <- conns_h$Importance
      index1 <- which(Imp_value > HC_cutoff) # HC
      index2 <- intersect(which(Imp_value > LC_cutoff), which(Imp_value <= HC_cutoff)) # MC
      index3 <- which(Imp_value <= LC_cutoff) # LC

      #### do data frame:gene, starts, end, peak1,peak2,importance,function_type
      function_type <- rep(NA, length(Imp_value))
      function_type[index1] <- "HC"
      function_type[index2] <- "MC"
      function_type[index3] <- "LC"
      # rescue highly correlated CREs
      if (rescued) {
          if (i <= length(TXs)) {
              X <- TXs[[i]]
              Y <- TYs[[i]]
              CPi <- abs(cor(t(X)))
              for (p in 1:nrow(CPi)) {
                CPi[p,p] <- 0
              }
              # focus on HC rows
              hic_index <- which(rownames(X) %in% conns_h$Peak2[index1])
              other_index <- which(rownames(X) %in% conns_h$Peak2[-index1])
              CPi_sub <- CPi[hic_index, other_index, drop = FALSE]
              flag_matrix <- matrix(0,nrow = nrow(CPi_sub), ncol = ncol(CPi_sub))
              flag_matrix[which(CPi_sub > 0.25)] <- 1
              correlated_index <- which(colSums(flag_matrix) > 0)
              if (!is.null(correlated_index)) {
                function_type[conns_h$Peak2 %in% rownames(X)[other_index[correlated_index]]] <- "HC"
              }
          }
      }
      DIRECT_NET_Result[[i]] <- cbind(data.frame(gene = focus_markers[i], Chr = Chr[i], Starts = Starts[i], Ends = Ends[i]),cbind(conns_h,function_type = function_type))
    }
  }
  DIRECT_NET_Result_all <- do.call(rbind,DIRECT_NET_Result)
  DIRECT_NET_Result_all$Starts <- as.numeric(DIRECT_NET_Result_all$Starts)
  DIRECT_NET_Result_all$Ends <- as.numeric(DIRECT_NET_Result_all$Ends)
  DIRECT_NET_Result_all$Importance <- as.numeric(DIRECT_NET_Result_all$Importance)
  # save result
  Misc(object, slot = 'direct.net') <- DIRECT_NET_Result_all
  return(object)
}

