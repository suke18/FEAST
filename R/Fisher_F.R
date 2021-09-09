########################################################################################
######### ---------- Gene-level significance investigation (scores) ---------- #########
########################################################################################
# We implement two ways for computing the gene-level scores.
# The rankings of the genes represent the significance at gene-level.
# The two methods are actually monotone. The rankings from fisher scores and F scores are essentially the same.
# The implementation is fast which will take seconds to obtain the results.


#' Calculate the gene-level fisher score.
#'
#' @param Y A gene expression matrix
#' @param classes The initial cluster labels NA values are allowed. This can directly from the \code{Consensus} function.
#' @return The score vector
cal_Fisher2 = function(Y, classes){
    #' This is from the paper https://arxiv.org/pdf/1202.3725.pdf
    #' Vector based calculation
    if(all(Y %%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Y = log2(sweep(Y, 2, L, FUN="/") + 1)
    }
    classes = as.factor(classes)
    Y = as.matrix(Y)
    if (is.null(classes) || length(classes)!=ncol(Y)){
        stop("please input right clustering labels ! ")
    }else{
        tb = table(classes)
        nks = as.numeric(tb)
    }
    row_ms = rowMeans(Y, na.rm = TRUE)
    row_mclasses = row_Varclasses = matrix(0, ncol = length(tb), nrow = nrow(Y))

    for (i in seq_len(length(tb))){
        tmp_class = names(tb)[i]
        cix = which(classes==tmp_class)
        tmp_Y = Y[, cix]
        if (length(cix)  > 1){
            row_mclasses[, i] = rowMeans(tmp_Y, na.rm = TRUE)
            row_Varclasses[, i] = rowVars(tmp_Y, na.rm = TRUE)
        }
    }
    scores = nks %*% t((row_mclasses - row_ms)^2) / (nks %*% t(row_Varclasses))
    return(scores)
}




#' Calculate the gene-level F score and corresponding significance level.
#'
#' @param Y A gene expression matrix
#' @param classes The initial cluster labels NA values are allowed. This can directly from the \code{Consensus} function.
#' @return The score vector
#' @examples
#' data(Yan)
#' cal_F2(Y, classes = trueclass)
#' @importFrom stats median
#' @export
cal_F2 = function(Y, classes){
    if (all(Y%%1 == 0)){
        L = colSums(Y) / median(colSums(Y))
        Ynorm = log2(sweep(Y, 2, L, FUN="/") + 1)
    }else{Ynorm = Y}
    cid = which(is.na(classes))
    if (length(cid) > 0){
        Ynorm = Ynorm[, -cid]
        classes = classes[-cid]
    }
    classes = as.factor(classes)
    unique_classes = unique(classes)
    k = length(unique(classes))
    row_class_mean = matrix(0, ncol = k, nrow = nrow(Y))

    row_mean = rowMeans(Ynorm)
    k = length(unique(classes))
    pb = txtProgressBar( min = 0, max = k, style = 3)
    for (i in seq_len(k)){
        setTxtProgressBar(pb, i)
        ix = which(classes == unique_classes[i])
        if (length(ix) > 1){
            tmp_mean =  rowMeans(Ynorm[, ix])
            row_class_mean[, i] = tmp_mean
        }else{
            row_class_mean[, i] = Ynorm[, ix]
        }
    }; close(pb)
    colnames(row_class_mean) = unique_classes
    ### make sure the classes are matched; otherwise, causing error ###
    table_class = table(classes)[unique_classes]
    BBS = table_class  %*% t((row_class_mean - row_mean)^2)
    BBS = BBS[1,]
    TSS = rowSums((Ynorm - row_mean) ^2)
    ESS = TSS - BBS
    df1 = k-1; df2 = ncol(Ynorm)-k
    F_scores = (BBS/df1) / (ESS/ df2)
    ps = pf(F_scores, df1, df2, lower.tail = FALSE)
    return(list(F_scores = as.numeric(F_scores), ps = ps))
}

