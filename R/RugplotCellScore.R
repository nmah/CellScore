# RugplotCellScore.R

# from CSReport_v1.4.3.R

#' RugplotCellScore
#'
#' This function will plot a rugplot of all CellScore values for each transition
#' selected in the cell.change data frame. The function will only plot the
#' scores for the test samples (annotated by the
#' cellscore$column sub_cell_type1). Standards are not included.
#' Samples are coloured by a secondary property, which must be a single column
#' in the cellscore data frame.
#' @param cellscore a data.frame of CellScore values as calculated by
#'   CellScore().
#' @param cell.change a data frame containing three columns, one for the
#'   start (donor) test and target cell type. Each row of the data.
#'   frame describes one transition from the start to a target cell type.
#' @param colour.by the name of the column in the cellscore argument
#'   that contains the secondary property.
#' @return This function outputs the plot on the active graphical device
#'   and returns invisibly NULL.
#' @keywords cellscore boxplot
#' @seealso \code{\link[CellScore]{CellScore}} for details on CellScore.
#' @export
#' @importFrom graphics axis layout legend par plot points rect
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
#'    ## Rugplot of CellScore, colour samples by transition induction method
#'    RugplotCellScore(cellscore, cell.change,
#'                     "transition_induction_method")
#'  }

RugplotCellScore <- function(cellscore, cell.change, colour.by=NULL) {

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)
    if (is.null(colour.by) || !(colour.by %in% colnames(cellscore))){
        colour.by = "transition_induction_method"
        warning(paste(fun.main, ":",
                      "The argument 'colour.by' was not specified or is not a
                      valid column of the argument 'cellscore'.",
                      "It is now set to 'transition_induction_method'."))
    }

    ############################################################################
    ## PART I. Extract and format the data for plotting
    ############################################################################
    ## Get the test CellScore-s from valid transitions defined by cell.change
    plot.data <- extractTransitions(cellscore, cell.change)
    ## Get the sec. property and replace NA with "none"
    sec.property.vector <- plot.data[, colour.by]
    sec.property.vector[is.na(sec.property.vector)] <- "none"
    ## Map the sec. property to colours
    map.table <- .colourMapping(sec.property.vector)

    ## Secondary property name for legend, replace "_" with space, and
    ## capitalize first letter of each word
    sec.property.name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",
                              gsub("_", " ", colour.by, fixed=TRUE),
                              perl=TRUE)


    plot.data <- data.frame(plot.data,
                            plot_group=paste0(plot.data$experiment_id, "_",
                                              plot.data$sub_cell_type1),
                            map.table, stringsAsFactors=FALSE)

    ## Sort the table by target, donor, study
    plot.data <- plot.data[order(plot.data$target, plot.data$donor_tissue,
                                 plot.data$experiment_id), ]

    group.df <- unique(plot.data[, c("plot_group", "experiment_id",
                                     "sub_cell_type1")])
    rownames(group.df) <- group.df$plot_group

    ############################################################################
    ## PART II. Plot
    ############################################################################
    .doPlot(plot.data, group.df, sec.property.name)

    invisible()
}

.doPlot <- function(plot.data, group.df, legend.title){
    n.groups <- nrow(group.df)
    old.par <- par(no.readonly=TRUE)

    ## Setup plot layout
    layout(matrix(c(1,1,2), nrow=1, ncol=3))
    par(mar=c(5, 7, 4, 5) + 0.1)

    ## A. Rugplot
    plot(0, 0,
         xlim=c(min(plot.data$CellScore), max(plot.data$CellScore)),
         ylim=c(0, n.groups),
         type="n",
         ylab="",
         axes=FALSE,
         frame.plot=TRUE,
         #  bty="o",
         #  ann=F,
         yaxp=c(1, n.groups, n.groups),
         #  yaxt="n",
         xlab="CellScore")
    ## Set the plot region to grey
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col="grey90")

    ## Set size of symbol "|" and axis labels
    cex.points <- ifelse(n.groups > 50, 0.8, 1)
    ## Add axes
    axis(side=1, labels=TRUE)
    axis(side=4, at=c(1:n.groups), las=2, cex.axis=cex.points,
         labels=group.df$experiment_id)
    axis(side=2, at=c(1:n.groups), las=2, cex.axis=cex.points,
         labels=group.df$sub_cell_type1)

    ## Add the points by groups
    lapply(seq_len(n.groups),
           function(i) {
               sel <- plot.data$plot_group %in%  group.df[i, "plot_group"]
               points(plot.data[sel, "CellScore"], rep(i, sum(sel)),
                      pch="|", col= plot.data[sel, "col"], cex=cex.points)
           })

    ## B. Legend
    plot(1,1, type="n", xaxt="n", yaxt="n", xlab="", ylab="",
         main="", bty="n", cex.main=0.8)
    col.table = unique(plot.data[,c("group","col")])
    legend("topleft", fill=col.table$col, legend=col.table$group,
           title=legend.title)

    ## Reset graphical parameters
    par(old.par)
}
