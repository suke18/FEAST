#' Using clustering results based on feature selection to perform model selection.
#'
#' @param model_cv_res model selection result from \code{Select_Model_short_SC3}.
#' @param trueclass The real class labels
#' @return a list of mse dataframe, clustering accuracy dataframe, and ggplot object.
#' @examples
#' \dontrun{
#' data(Yan)
#' k = length(unique(trueclass))
#' Y = process_Y(Y, thre = 2) # preprocess the data
#' con_res = Consensus(Y, k=k)
#' mod_res = Select_Model_short_SC3(Y, cluster = con_res$cluster, top = c(200, 500, 1000, 2000))
#' library(ggpubr)
#' Visual_Rslt(model_cv_res = mod_res, trueclass = trueclass)
#' }
#' @export
Visual_Rslt = function(model_cv_res, trueclass){
    mse = model_cv_res$mse
    sc3_res = model_cv_res$SC3_res
    if (missing(trueclass) || length(sc3_res[[1]]$cluster) != length(trueclass)){
        stop("Please input the right true classes !")
    }
    accur = sapply(sc3_res, function(x) eval_Cluster(x$cluster, trueclass))
    accur = accur[c(1,2,3), ]
    df1 = data.frame(lens = names(mse), mse = as.numeric(mse))
    g1 = ggbarplot(df1, x = "lens", y = "mse", ylab = "MSE", fill = "lightblue", color = "grey80",
              ylim = c(min(df1$mse), max(df1$mse)),
              title = "Validation by MSE", font.main = c(18, "bold"),
              xlab = "# top features")+
        theme(plot.title = element_text(size=18, face="bold"),
              legend.position = "top", legend.box.just = "left",
              legend.title = element_text(size = 13, face = "bold"),
              axis.title.x=element_text(size=16, face='bold'),
              axis.text.x=element_text(size=15,angle=270),
              axis.title.y=element_text(size=16,face='bold'),
              axis.text.y=element_text(size=15),
              legend.text=element_text(size=13))
    measures = rep(rownames(accur), ncol(accur))
    lens = rep(colnames(accur), each = nrow(accur))
    df2 = data.frame(ms = as.numeric(accur), measures = measures, lens = lens, stringsAsFactors = FALSE)
    df2$measures = factor(df2$measures, levels = unique(df2$measures))
    df2$lens = factor(df2$lens, levels = unique(df2$lens))
    g2 = ggline(df2, "lens", "ms",  legend= "bottom",
           linetype = "measures", size = 1, shape = "measures", palette = get_palette("npg", 3),
           xlab = "# top features", ylab = "measurement", legend.title = '', color = "measures", title = "FEAST by SC3") +
        scale_y_continuous(breaks =seq(0, 1, 0.1), labels =seq(0, 1, 0.1), limits =c(0, 1)) + theme(legend.position = c(0.45, 0.11), legend.direction = "horizontal") +
        theme(plot.title = element_text(size=18, face="bold"),
              legend.position = "top", legend.box.just = "left",
              legend.title = element_text(size = 13, face = "bold"),
              axis.title.x=element_text(size=16, face='bold'),
              axis.text.x=element_text(size=15,angle=270),
              axis.title.y=element_text(size=16,face='bold'),
              axis.text.y=element_text(size=15),
              legend.text=element_text(size=13)) + guides(col = guide_legend(nrow = 1))

    gg = ggarrange(plotlist=list(g2, g1),ncol=2, common.legend=FALSE)
    return(list(mse = df1, measure = df2, ggobj = gg))
}
