# scatterplotDonorTargetTest.R

# from CSReport_v1.4.3.R

#' scatterplotDonorTargetTest
#'
#' This function is called by CellScoreReport to make a scatterplot of test and
#' standard samples (donor and target).
#' @param test.data a data.frame of CellScore values as calculated by
#'   CellScore(), for a group of test samples.
#' @param cellscore a data.frame of CellScore values as calculated by
#'   CellScore().
#' @param index.plot a logical variable, with TRUE meaning sample index should
#'   be plotted for easy identification of spots. Default is FALSE.
#' @return This function outputs the plot on the active graphical device
#'   and returns invisibly NULL.
#' @keywords scatterplot donor target
#' @export
#' @importFrom graphics par plot text points legend
#' @importFrom grDevices densCols
#' @importFrom RColorBrewer brewer.pal
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
#'
#'    ## How many plot_groups are there?
#'    table(plot.data$plot_group)
#'
#'    ## pick one plot_group to plot
#'    group <- unique(plot.data$plot_group)[4]
#'
#'    ## Select scores for only one plot group
#'    test.data <- plot.data.ordered[plot.data.ordered$plot_group %in% group, ]
#'
#'    ## save current graphical parameters
#'	 old.par <- par(no.readonly=TRUE)
#'
#'    ## Plot: this will plot a 2-paneled plot
#'	 par(mfrow=c(1,2))
#'    scatterplotDonorTargetTest(test.data, cellscore, FALSE)
#'  
#'    ## Reset graphical parameters
#'    par(old.par)
#'    
#' }
#' }

scatterplotDonorTargetTest <- function(test.data, cellscore, index.plot=FALSE) {

    ###########################################################################
    ## PART 0. Check function arguments
    ###########################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(test.data, 'test.data', fun.main)
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)

    ############################################################################
    ## PART I. Data preparation
    ############################################################################
    celltype <- list(donor=test.data$start[1],
                     target=test.data$target[1],
                     test=test.data$sub_cell_type1[1])

    ## Collect all data tables in one list
    ## Subset columns using column names: "donor.like","target.like","index"
    sel.cols <- match(c("donor.like","target.like","index"),
                      colnames(cellscore))
    data.list <- list(target=NULL, donor=NULL, test=NULL)
    for (group in names(data.list)) {
        if (group == "test"){
            data.list[[group]] <- test.data[, sel.cols]
        }else{
            sel.trans <- cellscore$start == celltype$donor &
                cellscore$target == celltype$target
            sel.group <- cellscore$general_cell_type == celltype[[group]]
            data.list[[group]] <- cellscore[ sel.trans & sel.group, sel.cols]
        }
    }

    ## Colour mapping
    col.table <- .colourMapping(test.data$sub_cell_type1)

    ############################################################################
    ## PART II. Plot
    ############################################################################
    ## A. Scatter plot
    ## Set boundaries so that the x- and y-axes have the same scale
    ## Limits should be based on the standards
    xylim <- c(min(c(0.7, unlist(data.list[c("donor", "target")]))), 2)

    par(mar=c(4, 6, 4, 2) + 0.1)
    plot(xylim, xylim, type="n",
         xlab="Donor-like score", ylab="Target-like score",
         cex.lab=1.5, cex.axis=1.3, cex.main=1.5,
         main="CellScore Components", xlim=xylim, ylim=xylim )

    lapply(names(data.list),
           function(group){
               ## 1. Set the colour and plotting symbol
               if (group == "test"){
                   col.group <- col.table$col
                   pch.group <- 4
               }else{
                   col.group <-
                       try(densCols(data.list[[group]],
                                    colramp=.getMainColours(group, TRUE)),
                           silent=TRUE)
                   if (class(col.group) == "try-error") {
                       col.group <- .getMainColours(group, FALSE)
                   }
                   pch.group <- 20
               }

               ## 2. Show samples by group
               points(data.list[[group]][, -3], # remove the index
                      col=col.group, pch=pch.group, cex=1.5)

               ## 3. Add text annotaions for tracking samples
               if (index.plot) {
                   text(data.list[[group]][, -3],
                        labels=data.list[[group]][, 3], cex=0.7)
               }

           })

    ## B. Legend
    ## Plot it in a 2nd field since it could be very large
    plot(test.data$donor.like, test.data$target.like, type="n", xaxt="n",
         yaxt="n", xlab="", ylab="", bty="n", main="", cex.main=0.8)

    col.table = unique(col.table)
    leg.vector <- sapply(names(celltype),
                         function(x){
                             if (x != "test"){
                                 type <- celltype[[x]]
                             }else{
                                 type <- col.table$group
                             }
                             paste0(type, " (", nrow(data.list[[x]]),")")
                         })
    legend("topleft",
           fill=c(.getMainColours("donor", FALSE),
                  .getMainColours("target", FALSE),
                  NA),
           border=FALSE,
           legend=unlist(leg.vector),
           pch=c(NA, NA, 4), col=c(NA, NA, col.table$col),
           pt.cex=1.5, cex=1.1,
           title=paste0(test.data$experiment_id[1], ": ", "Transition from ",
                        celltype$donor," -> ", celltype$target))

    invisible()
}
