# PcaStandards.R

# version from CosineScore_v1.4.1.R


#' PCA plot on the most variable portion of the standard expression dataset
#'
#' This function will generate a principal component analysis (PCA) plot of
#' the IQR-filtered expression values that were used to generate the cosine
#' similarity scores.
#' @param label vector to be used for the point colours
#' @param label.name name of the label
#' @param exps an expression matrix of the IQR-filtered values as
#' obtained by the function CosineSimScore().
#' @param text.label a vector of characters to label each point.
#' @param col.palette a vector of colours to be used. There are 41
#' default colours.
#' @return The function will plot two panels, a PCA plot on the left and a
#' legend on the right. This is to accommodate that fact that the cell types
#' names are NOT abbreviated and the legend might not fit in the plot area.
#' @keywords pca
#' @seealso \code{\link[CellScore]{CosineSimScore}} for details on cosine
#'   similarity calculation.
#' @export
#' @importFrom stats prcomp sd
#' @importFrom graphics layout legend plot points text
#' @examples
#' ## Load the expression set for the standard cell types
#' library(Biobase)
#' library(hgu133plus2CellScore) # eset.std
#'
#' ## Locate the external data files in the CellScore package
#' rdata.path <- system.file("extdata", "eset48.RData", package = "CellScore")
#' tsvdata.path <- system.file("extdata", "cell_change_test.tsv",
#'                             package = "CellScore")
#'
#' if (file.exists(rdata.path) && file.exists(tsvdata.path)) {
#'
#'    ## Load the expression set with normalized expressions of 48 test samples
#'    load(rdata.path)
#'
#'    ## Import the cell change info for the loaded test samples
#'    cell.change <- read.delim(file= tsvdata.path, sep="\t",
#'                              header=TRUE, stringsAsFactors=FALSE)
#'
#'    ## Combine the standards and the test data
#'    eset <- combine(eset.std, eset48)
#'
#'    ## Generate cosine similarity for the combined data
#'    ## NOTE: May take 1-2 minutes on the full eset object
#'    ## so we subset it for 4 cell types
#'    pdata <- pData(eset)
#'    sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER")
#'    eset.sub <- eset[, sel.samples]
#'    cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)
#'
#'    PcaStandards(cs$pdataSub$experiment_id, "Experiment ID", cs$esetSub.IQR)
#' }

PcaStandards <- function(label, label.name, exps,  text.label=NULL,
                         col.palette=c("blue", "magenta", "green",
                                       "red", "goldenrod",
                                       "mediumslateblue", "olivedrab",
                                       "navyblue", "plum",
                                       "tomato", "thistle",
                                       "limegreen", "burlywood4",
                                       "cornflowerblue", "deeppink",
                                       "chartreuse", "forestgreen",
                                       "darkslateblue", "blueviolet",
                                       "gray50", "darkorange",
                                       "black", "lightsalmon4",
                                       "mediumseagreen", "palegreen4",
                                       "palevioletred4","peachpuff4",
                                       "plum4", "mediumspringgreen",
                                       "darkred", "khaki4","lawngreen",
                                       "lightseagreen","orange",
                                       "orchid3", "sienna4","snow4",
                                       "turquoise3","wheat3",
                                       "goldenrod2","darkorange3")) {

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    if (!is.vector(label) || !is.matrix(exps) || length(label) != ncol(exps))
        stop(paste("The argument 'label' (a vector) should have",
                   "as many as elements as columns (samples)",
                   "in the argumnet 'exps' (an expression matrix),",
                   "exiting function", fun.main))


    ############################################################################
    ## PART I. Data preparation, PCA and colour mapping
    ############################################################################
    ## Remove any rows that have zero variance
    sel <- apply(exps, 1, function(x) sd(x) != 0)
    plot.me <- exps[sel,]

    ## Because the input matrix has values from 0-1
    ## these values are YuGene normalized, DONT USE SCALING
    pca <- prcomp(t(plot.me), scale=FALSE)

    ## Extract coordinates from pca object
    pca.comp <- pca$x[, c(1, 2)]

    ## Get proportion of variance explained by PC1 and PC2
    pca.sum <- summary(pca)
    pc1 <- pca.sum$importance[2,1] *100
    pc2 <- pca.sum$importance[2,2] *100

    ## Colour mapping
    map.table <- .colourMapping(label, col.palette)

    ############################################################################
    ## PART II. Plot
    ############################################################################
    old.par <- par(no.readonly=TRUE)

    ## Set up two-plot layout
    layout(matrix(c(1,2), nrow=1, ncol=2))

    ## A. PCA plot
    plot(pca.comp, type="n",
         xlab=paste0("PC1 (", round(pc1, 1), "%)"),
         ylab=paste0("PC2 (", round(pc2, 1), "%)"),
         main ="",
         sub=paste(nrow(exps), "probes;", ncol(exps),"samples"))

    ## Plot the samples as points & mark them with text labels if provided
    ## Note for NM: removed the for-loop and the revrese order of the vector
    points(pca.comp[, 1], pca.comp[, 2], col=map.table$col, pch=16, cex=2)
    if (length(text.label) == nrow(pca.comp)){
        text(pca.comp[, 1], pca.comp[, 2], labels=text.label, col="black",
             cex=0.6)
    }

    ## B. Legend
    plot(pca.comp, type="n", xaxt="n", yaxt="n", xlab="", ylab="",
         main="", bty="n", cex.main=0.8)

    ## Set number of columns based on number of groups in legend
    legend.columns <- round(length(levels(as.factor(label))) / 30)
    map.table <- unique(map.table)
    map.table <- map.table[order(map.table$group), ]
    legend("topleft", cex=0.7,
           legend=map.table$group,
           fill=map.table$col,
           ncol=legend.columns, title=label.name)

    ## Reset graphical parameters
    par(old.par)
}
