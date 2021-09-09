#' Align the cell types from the prediction with the truth.
#'
#' @param tt0 a N*N table.
#' @return the matched (re-ordered) table
#' @examples
#' vec1 = rep(1:4, each=100)
#' vec2 = sample(vec1)
#' tb = table(vec1, vec2)
#' #tb_arg = align_CellType(tb)
#' @importFrom methods is
#' @export
align_CellType = function(tt0) {
    N_rows = nrow(tt0)
    N_cols = ncol(tt0)
    tt = tt0
    if(N_rows<N_cols) {
        tt = t(tt0)
        N_rows = nrow(tt)
        N_cols = ncol(tt)
    }

    K = ncol(tt)
    res = array(NA, dim=dim(tt))

    ## first find the largest entry
    ix.row = rep(NA, nrow(tt))
    ix.col = rep(NA, ncol(tt))
    for(i in seq_len(K)) { ## loop through the columns. Note that now we have N_rows<N_cols
        iii = which(tt == max(tt), arr.ind = TRUE)
        if(is.matrix(iii)) ## multiple max numbers
            iii = iii[1,]
        ix.row[i] = iii[1]
        ix.col[i] = iii[2]
        tt[iii[1],] = -1
        tt[,iii[2]] = -1
    }
    if(N_rows>N_cols) { # fill in extra rows
        ix.row[(K+1):N_rows] = setdiff(seq_len(N_rows), ix.row)
    }
    tt = tt0
    if(nrow(tt0)<ncol(tt0))
        tt = t(tt0)
    res = tt[ix.row, ix.col]
    if(nrow(tt0)<ncol(tt0))
        res = t(res)
    return(res)

}



#' Calculate the purity between two vectors.
#'
#' @param x a vector.
#' @param y a vector. x and y are with the same length.
#' @return the purity score
Purity = function(x, y) {
    x = as.vector(x)
    y = as.vector(y)

    tbl = table(x, y)
    tbl.major = apply(tbl, 1, max)
    n = sum(tbl)

    purity = sum(tbl.major)/ n
    return(purity)

}



#' Calculate 3 metrics and these methods are exported in C codes.
#' flag = 1 --- Rand index, flag = 2 --- Fowlkes and Mallows's index, flag = 3 --- Jaccard index
#'
#' @param cl1 a vector
#' @param cl2 a vector
#' @param randMethod a string chosen from "Rand", "FM", or "Jaccard"
#' @return a numeric vector including three values
cal_metrics = function(cl1, cl2, randMethod = c("Rand", "FM", "Jaccard"))
{
    if(!is.vector(cl1)){
        stop("cl1 is not a vector!\n");
    }
    if(!is.vector(cl2)){
        stop("cl2 is not a vector!\n");
    }
    if(length(cl1) != length(cl2)){
        stop("two vectors have different lengths!\n");
    }

    len <- length(randMethod)
    if(len == 0){
        stop("The argument 'randMethod' is empty!\n") }

    # unique values of elements in 'cl1'
    cl1u <- unique(cl1)
    # number of clusters in partition 1
    m1 <- length(cl1u)

    # unique values of elements in 'cl2'
    cl2u <- unique(cl2)
    # number of clusters in partition 2
    m2 <- length(cl2u)

    n <- length(cl1)
    randVec <- rep(0, len)
    names(randVec) <- randMethod
    for(i in seq_len(len))
    {
        randMethod[i] <- match.arg(arg = randMethod[i], choices = c("Rand", "FM", "Jaccard"))

        flag <- match(randMethod[i], c("Rand", "FM", "Jaccard"))

        c.res <- .C("cal_3_metrics",
                    as.integer(cl1),
                    as.integer(cl1u),
                    as.integer(cl2),
                    as.integer(cl2u),
                    as.integer(m1),
                    as.integer(m2),
                    as.integer(n),
                    as.integer(flag),
                    r = as.double(0))
        randVec[i] <- c.res$r
    }
    return(randVec)
}


aricode_NMI = function (c1, c2, variant =
                            c("max", "min", "sqrt", "sum", "joint")){
    variant <- match.arg(variant)
    entropy = function (c1, c2) {
        res <- sortPairs(c1, c2)
        N <- length(c1)
        H.UV <- -sum(res$nij * log(res$nij))/N + log(N)
        H.U <- -sum(res$ni. * log(res$ni.))/N + log(N)
        H.V <- -sum(res$n.j * log(res$n.j))/N + log(N)
        res <- list(UV = H.UV, U = H.U, V = H.V, sortPairs = res)
        res
    }
    H <- entropy(c1, c2)
    MI <- -H$UV + H$U + H$V
    D <- switch(variant, max = max(H$U, H$V), sqrt = sqrt(H$U *H$V),
                min = min(H$U, H$V), sum = 0.5 * (H$U + H$V), joint = H$UV)
    res <- MI/D
    res
}




#' Calculate the a series of the evaluation statistics.
#'
#' @param vec1 a vector.
#' @param vec2 a vector. x and y are with the same length.
#' @return a vector of evaluation metrics
#' @examples
#' vec2 = vec1 = rep(1:4, each = 100)
#' vec2[1:10] = 4
#' acc = eval_Cluster(vec1, vec2)
#' @export
eval_Cluster = function(vec1, vec2){
    vec1 = as.numeric(as.factor(vec1))
    vec2 = as.numeric(as.factor(vec2))
    ari = mclust::adjustedRandIndex(vec1, vec2)
    metrics3 = cal_metrics(vec1, vec2)
    # ri = metrics3[1]
    fmi = metrics3[2]
    jaccard = metrics3[3]
    # nmi = aricode_NMI(vec1, vec2)
    purity = Purity(vec1, vec2)
    res = as.numeric(c(ari, purity, jaccard, fmi))
    names(res) = c("ARI", "Purity", "Jaccard", "FM")
    return(res)
}








