# rugplotDonorTargetTest.R

# from CSReport_v1.4.3.R

#' rugplotDonorTargetTest
#'
#' This function is called by CellScoreReport to make a rugplot showing the
#' CellScore of all test samples, in relation to the standards. Donor and target
#' individual CellScore values are plotted in one horizontal lane, then test
#' CellScore values are are in another horizontal lane. Z-score cutoffs based
#' on the target standards are shown as dashed vertical lines.
#' @param test.data a data.frame of CellScore values as calculated by
#'   CellScore(), for only plot group of test samples.
#' @param cellscore a data.frame of CellScore values as calculated by
#'   CellScore().
#' @return This function outputs the plot on the active graphical device
#'   and returns invisibly NULL.
#' @keywords cellscore target rugplot test sample
#' @export
#' @seealso \code{\link[CellScore]{CellScore}} for details on CellScore.
#' @importFrom graphics abline axis legend par plot points rect rug text
#' @importFrom grDevices densCols
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats density sd
#' @examples
#' \dontrun{
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
#'    ## Generate the on/off scores for the combined data
#'    individ.OnOff <- OnOff(eset, cell.change, out.put="individual")
#'
#'    ## Generate cosine similarity for the combined data
#'    ## NOTE: May take 1-2 minutes on the full eset object
#'    cs <- CosineSimScore(eset, cell.change, iqr.cutoff=0.05)
#'
#'    ## Generate the CellScore values for all samples
#'    cellscore <- CellScore(data=eset, transitions=cell.change, scores.onoff=individ.OnOff$scores,
#'                           scores.cosine=cs$cosine.samples)
#'    ## Get the CellScore fvalues rom valid transitions defined by cell.change
#'    ## table
#'    plot.data <- extractTransitions(cellscore, cell.change)
#'
#'    ## Define a plot group variable
#'    plot.data$plot_group <- paste(plot.data$experiment_id,
#'                                 plot.data$cxkey.subcelltype, sep="_")
#'    ## Sort the scores 1) by target 2) by donor 3) by study
#'    plot.data.ordered <- plot.data[order(plot.data$target,
#'                                        plot.data$donor_tissue,
#'                                        plot.data$experiment_id), ]
#'    ## How many plot_groups are there?
#'    table(plot.data$plot_group)
#'
#'    ## pick one plot_group to plot
#'    group <- unique(plot.data$plot_group)[4]
#'
#'   ## Select scores for only one plot group
#'   test.data <- plot.data.ordered[plot.data.ordered$plot_group %in% group, ]
#'
#'   ## Plot
#'   rugplotDonorTargetTest(test.data, cellscore)
#'
#' }
#' }

rugplotDonorTargetTest <- function(test.data, cellscore) {

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(test.data, 'test.data', fun.main)
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)

    ############################################################################
    ## PART I. Data preparation
    ############################################################################
    ## Get the standard data for this transition and calculate Z-scores
    celltype <- list(donor=test.data$start[1],
                     target=test.data$target[1],
                     test=test.data$sub_cell_type1[1])
    data.list <- zscore.list <- list(target=NULL, donor=NULL, test=NULL)
    for (group in names(data.list)) {
        if (group == "test"){
            data.list[[group]] <- test.data$CellScore
        }else{
            sel.trans <- cellscore$start == celltype$donor &
                cellscore$target == celltype$target
            sel.group <- cellscore$general_cell_type == celltype[[group]]
            data.list[[group]] <- cellscore$CellScore[sel.trans & sel.group]
        }
        zscore.list[[group]] <- sapply(data.list[[group]],
                                       function(x) .zscore(data.list$target, x))
    }

    ############################################################################
    ## PART II. Plot
    ############################################################################
    ## A. Density plot of CellScore
    .doDensityPlot(data.list, celltype)

    ## B. Another rug plot representation: just plot lines [BETTER]
    title <- paste0(test.data$experiment_id[1], ": ", "Transition from ",
                    celltype$donor," -> ", celltype$target)
    .doRugPlot(data.list, zscore.list, celltype, title)

    invisible()
}

## doDensityPlot
##
## Local function that plots density curves of CellScore values (data.list)
## by groups defined in celltype.

.doDensityPlot <- function(data.list, celltype){
    dens.list <- lapply(data.list, density)
    plot(dens.list$target, type="n",
         main=paste0("CellScores: ","Transition from ",
                     celltype$donor," -> ", celltype$target),
         ylim=c(0, max(dens.list$target$y)),
         xlim=c(min(dens.list$donor$x), max(dens.list$target$x)))

    ## Set colours and add coloured points
    col.list <- .getMainColours("all", FALSE)
    lapply(names(dens.list),
           function(group){
               ## Plot density if there is more than one value
               if (length(data.list[[group]]) > 1) {
                   points(dens.list[[group]], col=col.list[[group]], type="l")
               }
           })

    ## Add rug of test scores
    rug(data.list$test, col=col.list$test)

    ## Add legend
    leg.vector <- sapply(names(col.list),
                         function(x){
                             paste0(celltype[[x]], " (",
                                    length(data.list[[x]]),")")
                         })
    legend("top", fill=unlist(col.list), legend=leg.vector)
}

## dorugPlot
##
## Local function that displays rugplot of CellScore values (data.list)
## by groups defined in celltype. It actually plots the Z-score of CellScore
## values (zscore.list). Title and should be provided.

.doRugPlot <- function(data.list, zscore.list, celltype, title){

    xlim <- c(min(unlist(data.list[c("target","test")])),
              max(unlist(data.list[c("target","test")])))
    ylim <- c(0, 3.5)
    par(mar=c(5, 5, 4, 5) + 0.1)

    ## Set up plot area
    plot(zscore.list$test, zscore.list$test,
         xlim=xlim,
         ylim=ylim,
         type="n",
         ylab="",
         axes=FALSE,
         frame.plot=TRUE,
         main="CellScores of Target and Test Cell Type",
         yaxp=c(1, rep(length(data.list), 2)),
         xlab="CellScore"
    )
    axis(side=1, labels=TRUE)
    axis(side=2, at=c(1:2), labels=c(celltype$test, "Standards"), las=2,
         cex.axis=0.9)

    ## Set the plot region to grey
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = "grey90")

    ## Standards in one line, use density colouring
    ## Set colours and add coloured points
    colramp.list <- .getMainColours("all", TRUE)
    lapply(names(data.list),
           function(group){
               col.group <- colramp.list[[group]]
               y.value <- 1
               if (group != "test"){
                   col.group <- densCols(data.list[[group]], colramp=col.group)
                   y.value <- 2
               }
               points(data.list[[group]],
                      rep(y.value, length(data.list[[group]])),
                      col=col.group,
                      pch="|",
                      cex=1.8)
           })

    ## Labels range of z-scores for the standard
    target.mean <- mean(data.list$target)
    target.std <- sd(data.list$target)
    abline(v=target.mean, lty=2, lwd=2, col="darkgreen")

    ## How many std to plot?
    zscore.std.plot <- abs(min(unlist(zscore.list[c("target","test")])))
    for (j in seq_len(zscore.std.plot)){

        list.sd = list(pos=target.mean + j*target.std,
                       neg=target.mean - j*target.std)

        for (element in list.sd){
            abline(v=element, lty=2, col="darkgreen")
            # plot points as bg for text
            points(element, rep(0,length(element)), pch=16, cex=1, col="white")
            # plot text
            text(element, rep(0, length(element)), labels=j, cex=0.7,
                 col="royalblue")
        }
    }

    ## Add legend
    col.list <- .getMainColours("all", FALSE)
    leg.vector <- sapply(names(col.list),
                         function(x){
                             paste0(celltype[[x]], " (",
                                    length(data.list[[x]]),")")
                         })
    legend("topleft", fill=unlist(col.list[c("target", "test")]),
           border=FALSE, bg="grey90", pt.cex=1.5, cex=1.1,
           legend=leg.vector[c("target", "test")],
           title=title)
}
