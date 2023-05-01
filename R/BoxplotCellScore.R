# BoxplotCellScore.R

# from CSReport_v1.4.3.R

#' Boxplot of the CellScore values across test samples
#'
#' This function will plot a boxplot of the CellScore values for each selected
#' transition (defined in the cell.change data frame). The function will only
#' plot the scores for the test samples of valid subtypes (as annotated by
#' cellscore$sub_cell_type1). Scores for the standards are not included. Note
#' that if a subtype is specified by two different transitions, the coresponding
#' scores will be plotted in both transitions.
#' @param cellscore a data.frame of CellScore values as calculated
#'   by CellScore().
#' @param cell.change a data frame containing three columns, one for the
#'   start (donor) test and target cell type. Each row of the data
#'   frame describes one transition from the start to a target cell type.
#' @seealso \code{\link[CellScore]{CellScore}} for details on CellScore
#'   calculation.
#' @return Invisibly, it returns list of the CellScore values by groups
#'   (in the same order as on the plot)
#' @keywords cellscore boxplot
#' @importFrom graphics axis par boxplot mtext
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
#'    sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER")
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
#'    ## Make the boxplot of CellScore values
#'    BoxplotCellScore(cellscore, cell.change)
#' }


BoxplotCellScore <- function(cellscore, cell.change) {

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)

    ############################################################################
    ## PART I. Extract and format the data for plotting
    ############################################################################
    ## Get the test CellScore-s from valid transitions defined by cell.change
    plot.data <- extractTransitions(cellscore, cell.change)

    ## remove any cell transitions that specify general_cell_types as the test
    ## cell types otherwise empty box plots will be plotted [KT: this should be
    ## moved into .extractTransitions ?!]
    sel <- !(cell.change$test %in% cellscore$general_cell_type)
    if (sum(sel) == 0){
        stop(paste("Extracting data resulted with no valid transitions",
                   "exiting function", fun.main))
    }
    cc <- cell.change[sel,]

    cc.order <- cc[order(cc$target, cc$start, decreasing=TRUE),]
    cc.order$cxkey <- apply(cc.order, 1, function(x) (paste(x, collapse="_")))
    cc.order$transition <- apply(cc.order[, c("start", "target")], 1,
                           function(x) (paste(x, collapse="->")))
    ## Map the test cell types to colours
    col.table <- .colourMapping(cc.order$test)
    cc.order$col <- col.table$col

    ## filter data by test cell type
    plot.list <- lapply(cc.order$cxkey,
                        function(x){
                            plot.data[plot.data$cxkey.subcelltype %in% x,
                                      "CellScore"]
                        })
    n.list <- length(plot.list)

    ############################################################################
    ## PART II. Plot
    ############################################################################
    old.par <- par(no.readonly = TRUE)

    par(mar=c(4, 8, 4,8) + 0.1) #c(bottom, left, top, right)
    cex.lab <- 1.5
    ## Boxplot
    boxplot(plot.list, las=2,
            main="",
            cex.lab=cex.lab,
            cex.main=1.5,
            names=cc.order$transition,
            col=cc.order$col,
            plot=TRUE,
            horizontal=TRUE,
            #     boxwex=0.4,
            ylim=c(-1.2,1.2),
            xlab="CellScore")

    ## Add axis labels
    axis(side=4, at=1:n.list, labels=cc.order$test, las=2)

    ## Text annotations: transitions on left side, sub_cell_type1 on right side
    mtext(side=4, expression(bold("Derived Celltype")), at=n.list+2,
          line=1, las=2, cex=1)
    mtext(side=2, expression(bold("Transition")), at=n.list+2,
          line=1, las=2, cex=1)

    ## Reset graphical parameters
    par(old.par)

    invisible(plot.list)
}
