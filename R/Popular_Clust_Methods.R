#########################################################################################
######### ---------- popular clustering methods: SC3, Seurat, TSCAN. ---------- #########
#########################################################################################
# These methods only involve one-step of feature selection.
# Among these, SC3 can specify the number of clusters, but the other two d
# I didn't consider CIDR and Monocle, because their performance is not good.



#' SC3 Clustering
#'
#' @param Y A expression matrix. It is recommended to use the raw count matrix.
#' @param k The number of clusters. If it is not provided, k is estimated by the default method in SC3.
#' @param input_markers A character vector including the featured genes. If they are not presented, SC3 will take care of this.
#' @return the clustering labels and the featured genes.
#' @importFrom matrixStats colVars
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import SC3
#' @export
SC3_Clust = function(Y, k = NULL, input_markers = NULL){
    # Y is the count matrix
    sce = SingleCellExperiment(
        assays = list(
            counts = as.matrix(Y),
            logcounts = log2(as.matrix(Y) + 1))
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
    if (! is.null(input_markers)){
        rowD = rowData(sce)
        sce = sc3_prepare(sce)
        sc3_gene_filter = rep(FALSE, nrow(sce))
        ix = match(input_markers, rownames(sce))
        col_sds = colVars(Y[ix,])
        col_ix = which(col_sds == 0)
        if (length(col_ix) > 0){
            for(col_id in col_ix){
                id_add = which.max(Y[, col_id])
                ix = c(ix, id_add)
            }
        }
        markers = rownames(Y)[ix]
        sc3_gene_filter[ix] = TRUE
        rowD$sc3_gene_filter = sc3_gene_filter
        markers = rownames(rowD)[rowD$sc3_gene_filter == TRUE]
        rowData(sce) = rowD
    }else{
        sce = sc3_prepare(sce)
        rowD = rowData(sce)
        markers = rownames(rowD)[rowD$sc3_gene_filter == TRUE]
    }
    if(is.null(k)){
        sce = sc3_estimate_k(sce)
        cat("estimated k clusters: ", metadata(sce)$sc3$k_estimation, "\n")
        k = metadata(sce)$sc3$k_estimation
    }
    sce = sc3_calc_dists(sce)
    sce = sc3_calc_transfs(sce)
    sce = sc3_kmeans(sce, ks = k)
    sce = sc3_calc_consens(sce)
    colTb = data.frame(colData(sce), stringsAsFactors = FALSE)
    cluster = colTb[,1]; names(cluster) = rownames(colTb)
    return(list(cluster = cluster, markers = markers))
}



#' TSCAN Clustering
#'
#' @param Y A expression matrix. It is recommended to use the raw count matrix.
#' @param k The number of clusters. If it is not provided, k is estimated by the default method in SC3.
#' @param input_markers A character vector including the featured genes. If they are not presented, SC3 will take care of this.
#' @param minexpr_percent minimum expression threshold (default = 0.5).
#' @param cvcutoff the cv cutoff to filter the genes (default = 1).
#' @return the clustering labels and the featured genes.
#' @examples
#'  data(Yan)
#'  k = length(unique(trueclass))
#'  # TSCAN_res = TSCAN_Clust(Y, k=k)
#' @import TSCAN
#' @importFrom stats lm
#' @export
TSCAN_Clust = function(Y, k, minexpr_percent = 0.5, cvcutoff = 1, input_markers = NULL) {
    # Here k is a vector for BIC calculation in Mclust function.
    if (all(Y %%1 == 0)){
        dat = log2(Y + 1)
    }else{
        dat = Y
    }
    # Remove genes with sd == 0. Didn't follow the original criterion
    if (is.null(input_markers)){
        dat = preprocess(dat, takelog = FALSE, minexpr_percent = minexpr_percent, cvcutoff = cvcutoff)
    }else{
        dat = dat[input_markers, ]
    }
    # from the TSCAN source code.
    markers = input_markers
    sdev <- prcomp(t(dat), scale = TRUE)$sdev[seq_len(20)]
    x <- seq_len(20)
    optpoint <- which.min(vapply(2:10, function(i) {
        x2 <- pmax(0, x - i)
        sum(lm(sdev ~ x + x2)$residuals^2)
    }, numeric(1)))
    pcadim = optpoint + 1
    tmpdata <- t(apply(dat, 1, scale))
    colnames(tmpdata) <- colnames(dat)
    tmppc <- prcomp(t(tmpdata), scale = TRUE)
    pcareduceres <- t(tmpdata) %*% tmppc$rotation[, seq_len(pcadim)]
    res <- suppressWarnings(Mclust(pcareduceres, G = k,
                                   modelNames = "VVV"))
    # This part is not robust
    if (is.null(res)){
        res = suppressWarnings(Mclust(pcareduceres, G = k, verbose = FALSE))
    }
    clusterid <- apply(res$z, 1, which.max)
    list(cluster = clusterid, markers = markers)
}



#' #' Seurat Clustering
#' #'
#' #' @param Y A expression matrix. It is recommended to use the raw count matrix.
#' #' @param knn of nearest neighbors (knn).
#' #' @param input_markers A character vector including the featured genes. If they are not presented, Seurat will take care of this.
#' #' @param resolution clustering resolution (0, 1) (default = 0.5).
#' #' @return the clustering labels and the featured genes.
#' Seurat_Clust = function(Y, knn = 10, resolution = 0.5, input_markers = NULL){
#'     #' input_markers is a character vector
#'     seuset = CreateSeuratObject(counts = Y)
#'     seuset = NormalizeData(object = seuset, normalization.method = "LogNormalize")
#'     if (! is.null(input_markers)){
#'         seuset@assays$RNA@var.features  = input_markers
#'         markers = input_markers
#'     }else{
#'         seuset = FindVariableFeatures(object = seuset, selection.method = "vst")
#'         markers = seuset@assays$RNA@var.features
#'     }
#'     seuset = ScaleData(object = seuset)
#'     seuset = RunPCA(object = seuset, features = VariableFeatures(object = seuset))
#'     if (! is.null(knn)){
#'         seuset = FindNeighbors(object = seuset, dims = seq_len(knn))
#'     }else{
#'         seuset = FindNeighbors(object = seuset)
#'     }
#'     if (! is.null(resolution)){
#'         seuset = FindClusters(object = seuset, resolution = resolution)
#'     }else{
#'         seuset = FindClusters(object = seuset)
#'     }
#'     cluster = Idents(seuset)
#'     tb = table(cluster)
#'     percentage = tb/sum(tb)
#'     return(list(cluster = cluster, percentage = percentage, markers = markers))
#' }


