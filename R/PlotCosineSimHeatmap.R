# from version CosineScore_v1.4.1.R

#' Plot heatmap of the cosine similarity score
#'
#' This function plots a triangular heatmap of the cosine similarity scores.
#' @param data a data.frame of cosine similarity scores, as generated
#' by the function CosineSimScore().
#' @param desc a single character, with description for the file name.
#'   Suggested are "general.groups", "subgroups", and "samples".
#' @param width the width of the output pdf, in inches.
#' @param height the height of the output pdf, in inches.
#' @param x the x-position of the heatmap legend. It may be necessary to change
#' the value to position the legend in a suitable place on the plot.
#' @param y the y-position of the heatmap legend. It may be necessary to change
#' the value to position the legend in a suitable place on the plot.
#' @return This function will print a pdf of the cosine similarity scores
#' in the current working directory.
#' @keywords cosine similarity score, heatmap
#' @seealso \code{\link[CellScore]{CosineSimScore}} for details on cosine
#'   similarity calculation.
#' @export
#' @importFrom grDevices pdf dev.off
#' @importFrom squash hkey bluered distogram
#' @importFrom stats dist hclust as.dendrogram order.dendrogram reorder
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
#'    ## NOTE: May take 1-2 minutes on the full eset object,
#'    ## so we subset it for 4 cell types
#'    pdata <- pData(eset)
#'    sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER", 
#'                  "ASC", "NPC", "MSC")
#'    eset.sub <- eset[, sel.samples]
#'    cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)
#'
#'    ## Generate pdf of cosine similarity heatmap in the working directory
#'    PlotCosineSimHeatmap(cs$cosine.general.groups, "general groups",
#'                         width=7, height=7, x=-3.5, y=1)
#' }


PlotCosineSimHeatmap <- function(data, desc="xx", width=20, height=20,
                                 x=-30, y=3) {

    ###########################################################################
    ## PART 0. Check function arguments
    ###########################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotSymetricMatrix0to1(data, "data", fun.main)

    ############################################################################
    ## PART I. Data preparation
    ############################################################################
    ## As the matrix is symetric, only cluster rows and extend to columns
    ## 1) USE average linkage for clustering and
    ## 2) reorder accroding to the row means (as done in heatmap.2)
    hcr <- hclust(dist(data), method="average")
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, rowMeans(data, na.rm=TRUE))
    rowInd <- order.dendrogram(ddr)
    if (nrow(data) != length(rowInd))
        stop(paste("Dendrogram ordering gave index of wrong length," ,
                   "exiting function", fun.main))
    plot.me <- data[rowInd, rowInd]

    ############################################################################
    ## PART II. Plot
    ############################################################################
    pdf(file=paste0("CosineSimilarityHeatmap_", gsub(" ", "_" , desc),".pdf"),
        width=width,
        height=height)

    ## Display a color-coded triangular distance matrix
    map <- distogram(plot.me, n=9, key=FALSE,
                     colFn=bluered,
                     main=paste("cosine similiarity for", desc,
                                "\n(average linkage-hierarchical clustering)"),
                     )
    ## Add colour key                )
    hkey(map, title = "cosine similarity", x=x, y=y, side=1)

    dev.off()
}
