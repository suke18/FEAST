## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
suppressMessages(library(FEAST))

## ----quick, eval = TRUE-------------------------------------------------------
data(Yan)
k = length(unique(trueclass))
Y = process_Y(Y, thre = 2)
ixs = FEAST(Y, k=k)
# look at the features
Ynorm = Norm_Y(Y)
par(mfrow = c(3,3))
for (i in 1:9){
  tmp_ix = ixs[i]
  tmp_gene = rownames(Ynorm)[tmp_ix]
  boxplot(as.numeric(Ynorm[tmp_ix, ])~trueclass, main = tmp_gene, xlab="", ylab="", las=2)
}

## ----load_data, eval = FALSE--------------------------------------------------
#  load("Pathto/Deng.RData")
#  Deng # load the Deng dataset, which includes 6 cell types (268 cells).
#  trueclass = colData(Deng)$cellTypes
#  k = length(unique(trueclass))
#  Y = assays(Deng)$counts

## ----consensus, eval = TRUE---------------------------------------------------
#Y = process_Y(Y, thre = 2) # preprocess the data if needed
con_res = Consensus(Y, k=k)

## ----gene-level, eval = TRUE--------------------------------------------------
F_res = cal_F2(Y, con_res$cluster)
ixs = order(F_res$F_scores, decreasing = TRUE) # order the features

## ----validation, eval = FALSE-------------------------------------------------
#  ## clustering step
#  tops = c(500, 1000, 2000)
#  cluster_res = NULL
#  for (top in tops){
#      tmp_ixs = ixs[1:top]
#      tmp_markers = rownames(Y)[tmp_ixs]
#      tmp_res = SC3_Clust(Y, k = k, input_markers = tmp_markers)
#      cluster_res[[toString(top)]] = tmp_res
#  }
#  ## validation step
#  Ynorm = Norm_Y(Y)
#  mse_res = NULL
#  for (top in names(cluster_res)){
#      tmp_res = cluster_res[[top]]
#      tmp_cluster = tmp_res$cluster
#      tmp_mse = cal_MSE(Ynorm = Ynorm, cluster = tmp_cluster)
#      mse_res = c(mse_res, tmp_mse)
#  }
#  names(mse_res) = names(cluster_res)

## ----demo, eval = FALSE-------------------------------------------------------
#  SC3_original = SC3_Clust(Y, k=k)
#  id = which.min(mse_res)
#  eval_Cluster(SC3_original$cluster, trueclass)
#  eval_Cluster(cluster_res[[id]]$cluster, trueclass)

## ----demopca, eval = FALSE----------------------------------------------------
#  SC3_original = SC3_Clust(Y, k=k)
#  id = which.min(mse_res)
#  eval_Cluster(SC3_original$cluster, trueclass)
#  eval_Cluster(cluster_res[[id]]$cluster, trueclass)

## ----multiple-group, eval = FALSE---------------------------------------------
#  Y = assays(Deng)$counts
#  trueclass = Deng$cellTypes
#  k = length(unique(trueclass))
#  Y = process_Y(Y, thre = 2) # preprocess the data
#  con_res = Consensus(Y, k=k)
#  mod_res = Select_Model_short_SC3(Y, cluster = con_res$cluster, top = c(200, 500, 1000, 2000))
#  # to visualize the result, one needs to load ggpubr library.
#  library(ggpubr)
#  Visual_Rslt(model_cv_res = mod_res, trueclass = trueclass)

## -----------------------------------------------------------------------------
sessionInfo()

