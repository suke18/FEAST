#########################################################################################
###### ------------------------- Process Expression data ------------------------- ######
#########################################################################################


#' Standard way to preprocess the gene expression data in SC3 method
#'
#' @param Y A gene expression data (Raw count matrix)
#' @return A processed gene expression matrix. It is already log transformed and normalized by CPM
log_TPM_SC3 = function(Y){
    Y = as.matrix(Y)
    #' Followed by the step from SC3 to do the gene filtering.
    #' > 2 counts and < 6%; expressed in more than 94% of the cells
    row_means = rowMeans(Y)
    row_expr_rate = rowMeans(Y!=0)
    ncells = ncol(Y)
    id1 = intersect(which(row_means > 2), which(row_expr_rate < 0.06))
    id2 = intersect(which(row_means > 0), which(row_expr_rate > 0.94))
    id3 = which(colSds(Y) < 1e-5)
    rem_id = Reduce(union, list(id1, id2, id3))
    if (! is.null(rem_id)){Y = Y[-rem_id, ]}
    L = colSums(Y)/1e6
    Y = log(sweep(Y, 2, L, FUN="/") + 1)
    return(Y)
}



#' Standard way to preprocess the count matrix. It is the QC step for the genes.
#'
#' @param Y A gene expression data (Raw count matrix)
#' @param thre The threshold of minimum number of cells expressing a certain gene (default =2)
#' @return A processed gene expression matrix. It is \emph{not log transformed}
#' @examples
#' data(Yan)
#' YY = process_Y(Y, thre=2)
#' @export
process_Y = function(Y, thre = 2){
    Y = as.matrix(Y)
    row_exprs_rate = rowMeans(Y != 0)
    row_sds = rowVars(Y)
    ncell = ncol(Y)
    rem_id1 = which(row_exprs_rate <= thre/ncell)
    rem_id2 = which(row_sds == 0)
    rem_id = union(rem_id1, rem_id2)
    if (length(rem_id) > 0) {Y = Y[-rem_id, ]}
    return(Y)
}





#' Standard way to preprocess the count matrix. It is the QC step for the genes.
#'
#' @param Ynorm A normalized gene expression matrix. If not, we will normalize it for you.
#' @param cluster The clustering outcomes. Specifically, they are cluster labels.
#' @param return_mses True or False indicating whether returning the MSE.
#' @return The MSE of the clustering centers with the predicted Y.
#' @examples
#' data(Yan)
#' Ynorm = Norm_Y(Y)
#' cluster = trueclass
#' MSE_res = cal_MSE(Ynorm, cluster)
#' @export
cal_MSE = function(Ynorm, cluster, return_mses = FALSE){
    Xregressor = model.matrix(~as.factor(cluster)-1)
    beta = solve(t(Xregressor)%*%Xregressor) %*% (t(Xregressor)%*%t(Ynorm))
    Yfit = Xregressor %*% beta
    res = t(Ynorm) - Yfit
    mse = mean(res^2)
    if (return_mses){
        mses = res^2
        res = list(mse = mse, mses = mses)
    }else{
        res = mse
    }
    return(res)
}


#' Normalize the count expression matrix by the size factor and take the log transformation.
#'
#' @param Y a count expression matrix
#' @return a normalized matrix
#' @examples
#' data(Yan)
#' Ynorm = Norm_Y(Y)
#' @export
Norm_Y = function(Y){
    L = colSums(Y)/median(colSums(Y))
    Ynorm = log(sweep(Y, 2, L, FUN = "/") + 1)
    return(Ynorm)
}



