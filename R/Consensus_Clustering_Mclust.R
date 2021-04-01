#' function for convert a vector to a binary matrix
#' @param vec a vector.
#' @return a n by n binary matrix indicating the adjacency.
vector2matrix = function(vec){
    mat = matrix(0, nrow = length(vec), ncol = length(vec))
    diag(mat) = diag(mat) + 1
    classes = unique(vec)
    for (class in classes){
        tmp_ix = which(vec == class)
        # find all pair index of a class
        pair_ix = t(combn(tmp_ix, 2))
        pair_ix = rbind(pair_ix, pair_ix[,c(2,1)])
        mat[pair_ix] = 1
    }
    return(mat)
}


#' Consensus Clustering
#'
#' @param Y A expression matrix. It is recommended to use the raw count matrix.
#' @param num_pcs The number of top pcs that will be investigated on through consensus clustering.
#' @param k The number of input clusters (best guess).
#' @param top_pctg Top percentage of features for dimension reduction
#' @param thred For the final GMM clustering, the probability of a cell belonging to a certain cluster.
#' @return the clustering labels and the featured genes.
#' @examples
#' data(Yan)
#' set.seed(123)
#' rixs = sample(nrow(Y), 500)
#' cixs = sample(ncol(Y), 40)
#' Y = Y[rixs, cixs]
#' con = Consensus(Y, k=5)
#' @export
Consensus = function(Y, num_pcs = 10, top_pctg = 0.33, k =2, thred = 0.9){
    if (all(Y %%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Y = log(sweep(Y, 2, L, FUN="/") + 1)
    }
    # select some genes (by top 50% cv) and do pca
    rm_ix = which(rowVars(Y) == 0)
    if (length(rm_ix) > 0) Y = Y[-rm_ix, ]
    row_ms = rowMeans(Y, na.rm = TRUE)
    row_sds = rowSds(Y, na.rm = TRUE)
    cv_scores = row_sds / row_ms
    gene_ranks = order(cv_scores, decreasing = TRUE, na.last = TRUE)
    top = round(nrow(Y) * top_pctg)
    ixs = gene_ranks[seq_len(top)]
    Y = Y[ixs, ]
    pca_out = prcomp(t(Y))
    # consensus clustering
    mat_res = matrix(0, ncol = ncol(Y), nrow = ncol(Y))
    message("start consensus clustering ...")
    pb = txtProgressBar(min = 0, max = num_pcs, style = 3)
    for (i in seq_len(num_pcs)){
        setTxtProgressBar(pb, i)
        pca_mat = pca_out$x[, seq_len(i)]
        if (i ==1){
            res = suppressWarnings(Mclust(pca_mat, G = k, modelNames = "V", verbose = FALSE))
        }else{
            res = suppressWarnings(Mclust(pca_mat, G = k, modelNames = "VVV", verbose = FALSE))
        }
        if (is.null(res)) res = suppressWarnings(Mclust(pca_mat, G = k, verbose = FALSE))
        clusterid = apply(res$z, 1, which.max)
        mat = vector2matrix(clusterid)
        mat_res = mat_res + mat
    }
    # final step of clustering
    res = suppressWarnings(Mclust(mat_res, G = k, modelNames = "VII", verbose = FALSE))
    if (is.null(res)){
        res = suppressWarnings(Mclust(mat_res, G = k, verbose = FALSE))
    }
    close(pb)
    cluster = apply(res$z, 1, function(x){
        id = which(x > thred)
        if (length(id) == 0){
            return(NA)
        }else{
            return(id)
        }
    })
    return(list(mat_res = mat_res, cluster = cluster))
}

