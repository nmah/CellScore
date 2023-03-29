# ScatterplotCellScoreComponents.R

# from CSReport_v1.4.3.R

#' Scatterplot of the the donor-like and target-like scores
#'
#' This function will plot the components of the CellScore, namely the donor-
#' like and the target-like scores. The function will only plot the scores for
#' the test samples (annotated by the cellscore$column sub_cell_type1).
#' Standards are not included.
#' @param cellscore a data.frame of CellScore values as calculated by
#'   CellScore()
#' @param cell.change a data.frame with 3 columns: start cell type, test cell
#' type, target cell type
#' @param index.plot a logical variable, with TRUE meaning sample index should
#'   be plotted for easy identification of spots. Default is FALSE.
#'   This is useful if you want to see where the samples are located on the
#'   plot.
#' @return This function outputs the plot on the active graphical device
#'   and returns invisibly NULL.
#' @keywords cellscore scatterplot
#' @seealso \code{\link[CellScore]{CellScore}} for details on CellScore.
#' @importFrom graphics par plot text legend
#' @export
#' @examples
#' ## Load the expression set for the standard cell types
#' library(Biobase)
#' library(hgu133plus2CellScore) # eset.std
#'
#' ## Locate the external data files in the CellScore package
#' rdata.path <- system.file("extdata", "eset48.RData", package = "CellScore")
#' tsvdata.path <- system.file("extdata", "cell_change_test.tsv",
#'                              package = "CellScore")
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
#'    sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER", 
#'                  "ASC", "NPC", "MSC", "iPS", "piPS")
#'    eset.sub <- eset[, sel.samples]
#'    cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)
#'
#'    ## Generate the on/off scores for the combined data
#'    individ.OnOff <- OnOff(eset.sub, cell.change, out.put="individual")
#'
#'    ## Generate the CellScore values for all samples
#'    cellscore <- CellScore(data=eset.sub, transitions=cell.change, scores.onoff=individ.OnOff$scores,
#'                           scores.cosine=cs$cosine.samples)
#'
#'    ## Make the scaterplot of CellScore components
#'    ScatterplotCellScoreComponents(cellscore, cell.change, FALSE)
#' }

ScatterplotCellScoreComponents <- function(cellscore, cell.change,
                                           index.plot=FALSE) {
    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)

    ############################################################################
    ## PART I. Extract and format the data for plotting
    ############################################################################
    ## Get the test CellScore-s from valid transitions defined by cell.change
    plot.data <- extractTransitions(cellscore, cell.change)
    ## Get colours and plot symbols
    map.table <- .group2ColourSymbolMapping(plot.data$sub_cell_type1)

    ############################################################################
    ## PART II. Plot
    ############################################################################
    old.par <- par(no.readonly=TRUE)

    par(mfcol=c(1,2))  # put plot on one side and legend on the other
    par(mar=c(4, 6, 4, 2) + 0.1) # c(bottom, left, top, right)

    ## Set axis range so that the x- and y-axes have the same scale
    the.min <- min(c(0.7, unlist(plot.data$donor.like, plot.data$target.like)) )
    xlim <- ylim <- c(the.min, 2)
    ## Set the title
    main.title <- sprintf("CellScore Components for %d test samples",
                          nrow(map.table))
    plot(plot.data$donor.like, plot.data$target.like, type="p",
         pch=map.table$pch, col=map.table$col,
         xlab="Donor-like score", ylab="Target-like score",
         cex=1,
         cex.lab=1.5,
         cex.axis=1.4,
         cex.main=1.5,
         main=main.title,
         xlim=xlim,
         ylim=ylim
    )
    ## Add labels, for tracking samples
    if (index.plot == TRUE) {
        text(plot.data$donor.like, plot.data$target.like,
             labels=plot.data$index, cex=0.7)
    }

    ## Add legend
    plot(plot.data$donor.like, plot.data$target.like, type="n", xaxt="n",
         yaxt="n", xlab="", ylab="", bty="n", main="", cex.main=0.8)
    map.table <- unique(map.table)
    legend("left", legend=map.table$group, col=map.table$col,
           pch=map.table$pch, cex=1.1, pt.cex=2,
           title="Derived Cell Type",
           ncol=ceiling(nrow(map.table)/23))

    ## Reset graphical parameters
    par(old.par)

    invisible()
}

## group2ColourSymbolMapping
##
## Local function that creates one data.frame with three columns: group, color
## and symbol. This will be used to map the groups (unique values in the input
##  vector) to the colour/symbol needed for the plot.

.group2ColourSymbolMapping <- function(data.in){

    col.df <- .colourMapping(data.in)
    size.groups <- length(unique(col.df$group))
    size.palette <- length(unique(col.df$col))

    ## Symbols will be recycled if there are more groups than colours in palette
    ## Change pch symbol after running out of colours
    temp.symbols <- rep(c(16, 4, 6, 3, 0), each=size.palette)
    size.symbols <- length(unique(temp.symbols))
    if (size.groups > size.palette*size.symbols ) {
        temp.symbols <- rep(temp.symbols,
                            ceiling(size.groups/(size.palette*size.symbols)))
    }

    ## Collect colour and symbol mappings in one table
    data.frame(group=col.df$group,
               col=col.df$col,
               pch=temp.symbols[as.numeric(factor(col.df$group))],
               stringsAsFactors=FALSE)
}

