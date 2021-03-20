######################################################################
##### -------- Consensus clustering a backup clustering -------- #####
######################################################################


#' FEAST main function
#'
#' @param Y A expression matrix. Raw count matrix or normalized matrix.
#' @param k The number of input clusters (best guess).
#' @param num_pcs The number of top pcs that will be investigated through the consensus clustering.
#' @param dim_reduce dimension reduction methods chosen from pca or svd.
#' @param split boolean. If T, using subsampling to calculate the gene-level significance.
#' @param batch_size when split is true, need to claim the batch size for spliting the cells.
#' @param BPPARAM parameters for BiocParallel. e.g. bpparam().
#' @return the rankings of the gene-significance.
#' @examples
#' data(Yan)
#' k = length(unique(trueclass))
#' ixs = FEAST(Y, k=k)
#' @export
FEAST = function (Y, k = 2, num_pcs = 10, dim_reduce = c("pca", "svd"), split = FALSE, batch_size =1000, BPPARAM=bpparam()){
    if (all(Y%%1 == 0)) {
        L = colSums(Y)/median(colSums(Y))
        Ynorm = log(sweep(Y, 2, L, FUN = "/") + 1)
    }else{
        Ynorm = Y
    }

    # dimention reduction part
    row_ms = rowMeans(Ynorm, na.rm = TRUE)
    row_sds = rowSds(Y, na.rm = TRUE)
    cv_scores = row_sds / row_ms
    gene_ranks = order(cv_scores, decreasing = TRUE, na.last = TRUE)
    # this top number of features for pca can be adjusted. Alternatively using top 1000 genes by rowMeans.
    top = round(nrow(Y)*0.33)
    gene_ixs = gene_ranks[1:top]
    YYnorm = Ynorm[gene_ixs, ]; ncells = ncol(Ynorm)
    dim_reduce = match.arg(dim_reduce)

    # starting dimention reduction.
    pc_res = switch(dim_reduce,
                    pca = prcomp(t(YYnorm))$x,
                    svd = svd(t(YYnorm))$u)

    # setup for parallel computing.
    # if (bpworkers(BPPARAM) ==1){
    #     num_cores = 2
    #     BPPARAM = SnowParam(workers = num_cores, type = "FORK")
    # }

    message("start consensus clustering ...")
    # consensus clustering for less cells (<5000)
    BPPARAM$progressbar = TRUE
    if (ncol(Ynorm) < 5000 & split == FALSE){
        mat_res = matrix(0, ncol = ncells, nrow = ncells)
        # write function for bplapply.
        bp_fun = function(i){
            tmp_pca_mat = pc_res[, 1:i]
            if (i == 1) {
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "V", verbose = FALSE))
            }
            else {
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
            }
            if (is.null(res)){
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
            }
            clusterid = apply(res$z, 1, which.max)
            return(clusterid)
        }
        pc_cluster = bplapply(1:num_pcs, bp_fun, BPPARAM = BPPARAM)
        pc_mat = lapply(pc_cluster, vector2matrix)
        con_mat = Reduce("+", pc_mat)

        # final clustering
        res = suppressWarnings(Mclust(con_mat, G = k, modelNames = "VII",  verbose = FALSE))
        if (is.null(res)) {
            res = suppressWarnings(Mclust(con_mat, G = k, verbose = FALSE))
        }
        cluster = apply(res$z, 1, function(x) {
            id = which(x > 0.95)
            if (length(id) == 0) {
                return(NA)
            }
            else {
                return(id)
            }
        })
        F_scores = cal_F2(Ynorm, classes = cluster)$F_scores
    }
    else{
        split_k = round(ncells / batch_size)
        chunk_ixs = suppressWarnings(split(sample(ncol(Y)), 1:split_k))
        # write function for bplapply.
        bp_fun = function(i){
            cell_ixs = chunk_ixs[[i]]
            tmp_pca_mat = pc_res[cell_ixs, 1:3]
            res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
            if (is.null(res)){
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
            }
            clusterid = apply(res$z, 1, which.max)
            return(clusterid)
        }
        chunk_cluster = bplapply(1:split_k, bp_fun, BPPARAM=BPPARAM)
        F_res_all = lapply(1:split_k, function(j){
            cell_ixs = chunk_ixs[[j]]
            tmp_mat = Ynorm[, cell_ixs]
            tmp_cluster = chunk_cluster[[j]]
            F_scores = cal_F2(tmp_mat, classes = tmp_cluster)$F_scores
            return(F_scores)
        })
        F_mat = Reduce(cbind, F_res_all)
        F_scores = rowMeans(F_mat, na.rm = TRUE)
    }
    ixs = order(F_scores, decreasing = TRUE)
    return(ixs)
}




#' FEAST main function (fast version)
#'
#' @param Y A expression matrix. Raw count matrix or normalized matrix.
#' @param k The number of input clusters (best guess).
#' @param num_pcs The number of top pcs that will be investigated through the consensus clustering.
#' @param split boolean. If T, using subsampling to calculate the gene-level significance.
#' @param batch_size when split is true, need to claim the batch size for spliting the cells.
#' @param BPPARAM parameters for BiocParallel. e.g. bpparam().
#' @return the rankings of the gene-significance.
#' @examples
#' data(Yan)
#' k = length(unique(trueclass))
#' res = FEAST_fast(Y, k=k)
#' @export
FEAST_fast = function (Y, k = 2, num_pcs = 10, split = FALSE, batch_size =1000, BPPARAM=bpparam()){
    if (all(Y%%1 == 0)) {
        L = colSums(Y)/median(colSums(Y))
        Ynorm = log(sweep(Y, 2, L, FUN = "/") + 1)
    }else{
        Ynorm = Y
    }

    # dimention reduction part
    row_ms = rowMeans(Ynorm, na.rm = TRUE)
    gene_ranks = order(row_ms, decreasing = TRUE, na.last = TRUE)
    # this top number of features for pca can be adjusted. Using 1000 from rowMeans for fast calculation
    top = 1000
    gene_ixs = gene_ranks[1:top]
    YYnorm = Ynorm[gene_ixs, ]; ncells = ncol(Ynorm)

    # using the fast version from irlba
    Ynorm_scale = scale(YYnorm, scale=TRUE, center=TRUE)
    pca_res = prcomp_irlba(t(Ynorm_scale), n=10)
    pc_res = pca_res$x

    # setup for parallel computing. register for bplapply
    # if (bpworkers(BPPARAM) ==1){
    #     num_cores = 2
    #     BPPARAM = SnowParam(workers = num_cores, type = "FORK")
    # }

    # pb = txtProgressBar(min=1, max=num_pcs-1,style=3)
    # f = function(){
    #     count = 0
    #     function(...) {
    #         count <<- count + length(list(...)) - 1
    #         setTxtProgressBar(pb,count)
    #         Sys.sleep(0.0001)
    #         flush.console()
    #         c(...)
    #     }
    # }

    message("start consensus clustering ...")
    # consensus clustering for less cells (<5000)
    BPPARAM$progressbar = TRUE
    if (ncells < 5000 & split == FALSE){
        mat_res = matrix(0, ncol = ncells, nrow = ncells)
        bp_fun = function(i){
            tmp_pca_mat = pc_res[, 1:i]
            if (i == 1) {
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "V", verbose = FALSE))
            }
            else {
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
            }
            if (is.null(res)){
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
            }
            clusterid = apply(res$z, 1, which.max)
            return(list(clusterid))
        }
        pc_cluster = bplapply(1:num_pcs, bp_fun, BPPARAM=BPPARAM)
        pc_mat = lapply(pc_cluster, function(x){
            vector2matrix(x[[1]])
        })
        con_mat = Reduce("+", pc_mat)
        message("!!! done !!!")
        # sometimes hclust can be unbalanced, pam can be slow. Thus, using kmeans of is an option.Might be influenced by the outliers.
        km = kmeans(con_mat, k)
        cluster = km$cluster
        # res = suppressWarnings(Mclust(con_mat, G = k, modelNames = "VII",  verbose = F))
        # if (is.null(res)) {
        #     res = suppressWarnings(Mclust(con_mat, G = k, verbose = F))
        # }
        # cluster = apply(res$z, 1, function(x) {
        #     id = which(x > 0.95)
        #     if (length(id) == 0) {
        #         return(NA)
        #     }
        #     else {
        #         return(id)
        #     }
        # })
        message("calculate F scores")
        F_scores = cal_F2(Ynorm, classes = cluster)$F_scores
    }
    else{
        split_k = round(ncells / batch_size)
        chunk_ixs = suppressWarnings(split(sample(ncol(Y)), 1:split_k))
        bp_fun = function(i){
            cell_ixs = chunk_ixs[[i]]
            tmp_pca_mat = pc_res[cell_ixs, 1:3]
            res = suppressWarnings(Mclust(tmp_pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
            if (is.null(res)){
                res = suppressWarnings(Mclust(tmp_pca_mat, G = k, verbose = FALSE))
            }
            clusterid = apply(res$z, 1, which.max)
            return(list(clusterid))
        }
        chunk_cluster = bplapply(1:split_k, bp_fun, BPPARAM = BPPARAM)
        F_res_all = lapply(1:split_k, function(j){
            cell_ixs = chunk_ixs[[j]]
            tmp_mat = Ynorm[, cell_ixs]
            tmp_cluster = chunk_cluster[[j]]
            F_scores = cal_F2(tmp_mat, classes = tmp_cluster)$F_scores
            return(F_scores)
        })
        F_mat = Reduce(cbind, F_res_all)
        F_scores = rowMeans(F_mat, na.rm = TRUE)
    }
    ixs = order(F_scores, decreasing = TRUE, na.last = TRUE)
    return(ixs)
}
