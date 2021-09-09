#' Using clustering results based on feature selection to perform model selection.
#'
#' @param Y A gene expression matrix
#' @param tops A numeric vector containing a list of numbers corresponding to top genes; e.g., tops = c(500, 1000, 2000).
#' @param cluster The initial cluster labels NA values are allowed. This can directly from the \code{Consensus} function.
#' @return mse and the SC3 clustering result.
#' @examples
#' data(Yan)
#' k = length(unique(trueclass))
#' Y = process_Y(Y, thre = 2) # preprocess the data
#' set.seed(123)
#' rixs = sample(nrow(Y), 500)
#' cixs = sample(ncol(Y), 40)
#' Y = Y[rixs, cixs]
#' con_res = Consensus(Y, k=k)
#' # not run
#' # mod_res = Select_Model_short_SC3(Y, cluster = con_res$cluster, top = c(100, 200))
#' @importFrom stats median
#' @importFrom matrixStats colVars
#' @export
Select_Model_short_SC3 = function(Y, cluster, tops = c(500, 1000, 2000)){
    # make sure the labels and the cells are one-to-one match
    stopifnot(ncol(Y) == length(cluster))
    F_res = cal_F2(Y, classes = cluster)
    Fs = F_res$F_scores
    F_orders = order(Fs, decreasing = TRUE)
    Y = as.matrix(Y)
    k = length(unique(cluster[!is.na(cluster)]))
    if(all(Y %%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Ynorm = log2(sweep(Y, 2, L, FUN="/") + 1)
    }else{
        Ynorm = Y
    }
    message("Start the validation processing. ")
    ntop = length(tops)
    pb = txtProgressBar(min = 0, max = ntop, style = 3)
    mse = SC3_res = NULL
    for (i in seq_len(ntop)){
        setTxtProgressBar(pb, i)
        top = tops[i]
        ix = F_orders[seq_len(top)]
        # Need to check the column variance; otherwise cannot calculate the cell-cell distance, thus add one index.
        col_sds = colVars(Y[ix,])
        col_ix = which(col_sds == 0)
        if (length(col_ix) > 0){
            for(col_id in col_ix){
                id_add = which.max(Y[, col_id])
                ix = c(ix, id_add)
            }
        }
        markers = rownames(Y)[ix]
        ## predict the labels and fit the regression
        SC3_rslt = SC3_Clust(Y, k = k, input_markers = markers)
        SC3_res[[toString(top)]] = SC3_rslt

        # use linear regression for validation
        message("linear regression validation \n")
        mse = c(mse, cal_MSE(Ynorm, SC3_rslt$cluster))
    }
    close(pb)
    names(mse) = tops
    res = list(mse=mse,  SC3_res = SC3_res)
    return(res)
}



#' Using clustering results (from TSCAN) based on feature selection to perform model selection.
#'
#' @param Y A gene expression matrix
#' @param tops A numeric vector containing a list of numbers corresponding to top genes; e.g., tops = c(500, 1000, 2000).
#' @param cluster The initial cluster labels NA values are allowed. This can directly from the \code{Consensus} function.
#' @param minexpr_percent  The threshold used for processing data in TSCAN. Using it by default.
#' @param cvcutoff  The threshold used for processing data in TSCAN. Using it by default.
#' @return mse and the TSCAN clustering result.
#' @examples
#' data(Yan)
#' k = length(unique(trueclass))
#' Y = process_Y(Y, thre = 2) # preprocess the data
#' set.seed(123)
#' rixs = sample(nrow(Y), 500)
#' cixs = sample(ncol(Y), 40)
#' Y = Y[rixs, cixs]
#' con_res = Consensus(Y, k=k)
#' # not run
#' # mod_res = Select_Model_short_TSCAN(Y, cluster = con_res$cluster, top = c(100, 200))
#' @importFrom stats median
#' @export
Select_Model_short_TSCAN = function(Y, cluster, minexpr_percent = 0.5, cvcutoff = 1, tops = c(500, 1000, 2000)){
    # make sure the labels and the cells are one-to-one match
    stopifnot(ncol(Y) == length(cluster))
    F_res = cal_F2(Y, classes = cluster)
    message("Start calculating the gene-level significance. ")
    Fs = F_res$F_scores
    F_orders = order(Fs, decreasing = TRUE)
    Y = as.matrix(Y)
    k = length(unique(cluster[!is.na(cluster)]))
    if(all(Y %%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Ynorm = log2(sweep(Y, 2, L, FUN="/") + 1)
    }else{
        Ynorm = Y
    }
    message("Start the validation processing. ")
    ntop = length(tops)
    pb = txtProgressBar(min = 0, max = ntop, style = 3)
    mse = TSCAN_res = NULL
    tscan_original = TSCAN_Clust(Y, k, minexpr_percent = minexpr_percent, cvcutoff = cvcutoff)
    TSCAN_res[["origin"]] = tscan_original
    for (i in seq_len(ntop)){
        setTxtProgressBar(pb, i)
        top = tops[i]
        ix = F_orders[seq_len(top)]
        # # Need to check the column variance; otherwise cannot calculate the cell-cell distance, thus add one index.
        # col_sds = colVars(Y[ix,])
        # col_ix = which(col_sds == 0)
        # if (length(col_ix) > 0){
        #     for(col_id in col_ix){
        #         id_add = which.max(Y[, col_id])
        #         ix = c(ix, id_add)
        #     }
        # }
        markers = rownames(Y)[ix]
        ## predict the labels and fit the regression
        TSCAN_rslt = TSCAN_Clust(Ynorm, k = k, input_markers = markers)
        TSCAN_res[[toString(top)]] = TSCAN_rslt

        # use linear regression for validation
        message("linear regression validation \n")
        mse = c(mse, cal_MSE(Ynorm, TSCAN_rslt$cluster))
    }
    close(pb)
    names(mse) = tops
    res = list(mse=mse,  TSCAN_res = TSCAN_res)
    return(res)
}

