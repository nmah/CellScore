# OnOff.R

# from version OnOff_v1.3.1.R

#' On/off score
#'
#' This function calculates the on/off score for cell transitions. The score
#' takes into account the cell type spcefific and most variable portion of the
#' detected transcriptome. It can be calculated for a sample or group of
#' samples representing specific (standard or engineered) cell type.
#' @param eset an ExpressionSet containing data matrices of normalized
#'   expression data, present/absent calls, a gene annotation data frame and a
#'   phenotype data frame.
#' @param cell.change a data frame containing three columns, one for the start
#'   (donor) test and target cell type. Each row of the data frame describes one
#'   transition from the start to a target cell type.
#' @param out.put a character flag with two possible values, "marker.list" and
#'   "individual". The former means the on/off scores will be aggregated accross
#'   cell groups and also the marker genes for each cell transition
#'   (in cell.change) will be calculated, while the latter will generate
#'   the on/off scores for all individual samples.
#' @param min.diff.cutoff a real number that represents the minimum difference
#'   between the fraction of present calls in donor vs target (in the
#'   standards), in order to define the markers for a given cell transition.
#'   Default is 0.8.
#' @param test.cutoff a real number in (0, 1] that is the minimum fraction of
#'   present calls in a test sample/group to decide if a gene is present in a
#'   test sample/group. Default is stringently set at 0.95.
#' @return This function returns a list of two objects, as follows:
#' %\describe{
#'   \item{scores}{a data.frame of on/off scores for each cell group
#'    given in cell.change(out.put="marker.list") or for each individual sample
#'    (out.put="idividual")}
#'   \item{markers}{a list of marker genes for the selected cell transitions
#'   in cell.change (out.put="marker.list") or NULL (out.put="individual")}
#'   %}
#' @keywords onoff score markers
#' @export
#' @seealso \code{\link[hgu133plus2CellScore]{hgu133plus2CellScore}} for details on the
#'   specific ExpressionSet object that shoud be provided as an input.
#' @importClassesFrom Biobase ExpressionSet
#' @importMethodsFrom Biobase fData pData
#' @importFrom Biobase assayDataElement
#' @importFrom utils setTxtProgressBar txtProgressBar
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
#'    ## Generate a marker list
#'    group.OnOff <- OnOff(eset, cell.change, out.put="marker.list")
#'
#'    ## Calculate the on/off scores for individual samples
#'    individ.OnOff <- OnOff(eset, cell.change, out.put="individual")
#' }

OnOff <- function(eset, cell.change,
                  out.put = c("marker.list", "individual"),
                  min.diff.cutoff=0.8,
                  test.cutoff=0.95){

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    summarized_experiment <- .stopIfCantCoerceToSummarizedExperiment(eset, 'eset', fun.main)
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)
    .stopIfNotNumeric0to1(min.diff.cutoff, 'min.diff.cutoff', fun.main)
    .stopIfNotNumeric0to1(test.cutoff, 'test.cutoff', fun.main)

    ############################################################################
    ## PART I. Filter samples according to phenoData
    ##         Filter input cell transitions, acording to phenoData
    ##         o score comparisons should only contain abbreviations in the
    ##         start-test-target columns
    ##         o first check in the pdata table that all start, test and stop
    ##         groups exist in the table.
    ##         o if any groups are missing, remove these rows from the
    ##         score.comparisons table, so that ALL of the score
    ##         comparisons in the table have samples in all groups
    ############################################################################
    annot <- rowData(summarized_experiment)
    pdata <- .filterPheno(colData(summarized_experiment), fun.main, 'na')
    calls <- assay(summarized_experiment[, rownames(pdata)], "calls")
    rownames(calls) <- annot$probe_id
    score.comparisons <- .filterTransitions(cell.change, pdata, fun.main,
                                            "valid.names")

    ############################################################################
    ## PART II. Calculate on/off scores
    ############################################################################
    switch(out.put,
           ## A. Calculate individual OnOff Scores
           individual = .calculateIndOnOff(score.comparisons, calls, pdata,
                                           min.diff.cutoff, test.cutoff),
           ## B. Calculate OnOff group scores and marker lists for transitions
           marker.list = .calculateGroupOnOff(score.comparisons, calls, pdata,
                                              annot, min.diff.cutoff, test.cutoff))
}

## calculateIndOnOff
##
## Local function that calculates the on/off scores for each individual sample
.calculateIndOnOff <- function(score.comparisons, calls, pdata,
                               diff.cutoff, test.cutoff){

    parent.fun <- deparse(sys.calls()[[sys.nframe() - 1]])

    ## Select only unique (start, traget) pairs
    selected <- !duplicated(score.comparisons[, c("start", "target")])
    score.comparisons <- score.comparisons[selected, ]
    rownames(score.comparisons) <- 1:nrow(score.comparisons)
    message(paste( "To process", nrow(score.comparisons), "cell transitions"))

    ## For each transition/comparison get the test samples
    ## NOTE: NOW includes all available samples as test samples, so there are
    ## scores calculated for standard samples, used later to define the cutoffs
    test.samples <- rownames(pdata)
    test.size <- length(test.samples)

    out <- apply(score.comparisons[, c("start", "target")], 1,
                 function(instance){
                     ## Print message about the progress
                     message(sprintf("Now processing transition %s -> %s",
                                     instance["start"], instance["target"]))
                     ## Open a progress bar
                     pb <- txtProgressBar(min=0, max=test.size, style=3)

                     ## a0. Get the calls for the start and target cell type
                     calls.st <- .subsetCallsByGroup(instance, calls, pdata,
                                                     parent.fun)
                     ## b0. Get fractions of present calls for both cell types
                     frac.st <- lapply(calls.st, .calculateFractions)

                     temp <- lapply(seq_along(test.samples),
                                    function(i){
                                        ## Show progress bar
                                        setTxtProgressBar(pb, i)
                                        ## Current test sample
                                        test <- test.samples[i]
                                        ## a1. Get the calls for the test sample
                                        calls.test <- calls[, test, drop=FALSE]
                                        ## b1. Get fractions of present calls
                                        ##     for the test sample
                                        frac.test <-
                                            .calculateFractions(calls.test)
                                        fractions <- c(frac.st,
                                                       list(test=frac.test))
                                        ## c. Mark the cell type-specific
                                        ##     probestes
                                        probesets <-
                                            .getOnOffMarkersByGroup(fractions,
                                                                    diff.cutoff,
                                                                    test.cutoff)
                                        ## d. Calculate and format the on/off
                                        ##    scores
                                        scores <- .formatOnOffScores(probesets)

                                        data.frame(c(instance, test=test,
                                                     scores),
                                                   stringsAsFactors=FALSE)
                                    })
                     ## Close progres bar
                     close(pb)

                     data.frame(pdata[test.samples,], do.call("rbind", temp),
                                stringsAsFactors=FALSE)
                 })

    scores <- do.call("rbind", out)
    list(scores=scores, markers=NULL)
}

## calculateGroupOnOff
##
## Local function that calculates the on/off scores for each individual sample
.calculateGroupOnOff <- function(score.comparisons, calls, pdata, annot,
                                 diff.cutoff, test.cutoff){
    ##  Check that each standard groups have at least 3 samples
    ##  in each group, otherwise reject the comparison
    ##  NOTE: these groups need not be labeled as "standard", eg. iPS is
    ##  used as a start point but are also test cells
    parent.fun <- deparse(sys.calls()[[sys.nframe() - 1]])
    pdata.std <- pdata[pdata$category == "standard", ]
    score.comparisons <- .filterTransitions(score.comparisons, pdata.std,
                                            parent.fun, "sample.counts")
    if (length(unique(pdata$platform_id)) > 1){
        warning("Multiple array platforms exist in the phenotype data.")
    }

    ##  Loop through all selected comparisons to obtain the
    ##  corresponinding on/off marker probesets
    out <- apply(score.comparisons[, c("start", "target", "test")], 1,
                 function(instance){
                     ## a. Get the calls for each cell type
                     calls.list <- .subsetCallsByGroup(instance, calls, pdata,
                                                       parent.fun)
                     ## b. Get fractions of present calls for each cell type
                     fractions <- lapply(calls.list, .calculateFractions)
                     ## c. Mark the cell type-specific probestes
                     probesets <- .getOnOffMarkersByGroup(fractions, diff.cutoff,
                                                          test.cutoff)
                     ## d. Format the marker probesets into a data.frame
                     ##    and extend it with gene annotaions
                     markers <- .formatOnOffMarkers(instance, probesets, annot)
                     ## e. Calculate and format the on/off scores
                     scores <- .formatOnOffScores(probesets)
                     scores <- data.frame(c(instance, scores),
                                          stringsAsFactors=FALSE)

                     list(markers=markers, scores=scores)
                 })

    ## Combine marker and score lists
    onoff.markers <- do.call('rbind', lapply(out, '[[', 'markers'))
    sub.onoff.markers <- onoff.markers[, c("comparison", "group", "probe_id")]
    onoff.markers <- onoff.markers[!duplicated(sub.onoff.markers),]
    onoff.scores <- do.call('rbind', lapply(out, '[[', 'scores'))

    ## Final output
    list(scores=onoff.scores, markers=onoff.markers)
}

## subsetCallsByGroup
##
## Local function that subsets the calls matrix (given by the argumnet 'calls')
## for the probesets expressed in each cell group (type) invovled in a
## given cell transition (defined by the argument 'instance').
## Returns list of filtered call matrices.
.subsetCallsByGroup <- function(instance, calls, pdata, calling.fun){
    temp.list <- names(instance)
    names(temp.list) <- temp.list
    lapply(temp.list,
           function(group){
               pdata.sel <- .filterPheno(pdata, calling.fun,
                                         flag=ifelse(group == "test",
                                                     "anygroup", "group"),
                                         flag.value=instance[group])

               calls[, rownames(pdata.sel)]
           })
}

## getOnOffMarkersByGroup
##
## Local function that determines the marker probesets for a
## given cell transition, as the ones that pass certain cutoffs based on the
## expression summary (given as the argument 'fractions') over all cell samples.
## Returns list of filtered call matrices.
.getOnOffMarkersByGroup <- function(fractions, diff.cutoff, test.cutoff){

    sel <- vector(4, mode="list")
    names(sel) <- c("diff", "start", "target", "test")

    ## The distance between the start and target fractions
    sel[["diff"]] <- fractions[["start"]] - fractions[["target"]]

    ## Select the marker genes based on the given fraction
    ## difference cutoff
    sel[["start"]] <- sel[["diff"]] > diff.cutoff
    sel[["target"]] <- sel[["diff"]] < -diff.cutoff

    ## Select which genes are present in test matrix based on the
    ## given (arbitrary) cut.off. Default is 0.95, higher than this is
    ## very stringent. Alternative can be the highest fractions of
    ## the selected markers, like 3rd quantile
    sel[["test"]] <- fractions[["test"]] > test.cutoff
    sel
}

## formatOnOffMarkers
##
## Local function that formats the probeset markers into a tabular format,
## that included gene annotations (symbol, name and entrez id)
## Returns a data.frame
.formatOnOffMarkers <- function(instance, probesets, annotations){
    markers <- lapply(c("start", "target"),
                      function(group){
                          stopifnot(identical(names(probesets[[group]]),
                                              annotations$probe_id))
                          info <- paste(instance[c("start", "target")],
                                        collapse = "->")
                          count.group <- sum(probesets[[group]])
                          data.frame(comparison=rep(info, count.group),
                                     group=rep(instance[group], count.group),
                                     annotations[probesets[[group]],,drop=FALSE],
                                     stringsAsFactors=FALSE)
                      })
    do.call("rbind", markers)
}

## formatOnOffScores
##
##
## Local function that formats the on/off scores and its score components
## (loss and gain of marakers) into a tabular format
## Returns a data.frame
.formatOnOffScores <- function(probesets){
    start.count <- sum(probesets[["start"]])
    target.count <- sum(probesets[["target"]])
    test.start.count <- sum(probesets[["test"]] & probesets[["start"]])
    test.target.count <- sum(probesets[["test"]] & probesets[["target"]])
    start.loss <- 1 - test.start.count/start.count
    target.gain <- test.target.count/target.count

    data.frame(markers.start=start.count,
               markers.target=target.count,
               start.mkrs.in.test=test.start.count,
               target.mkrs.in.test=test.target.count,
               loss.start.mkrs=start.loss,
               gain.target.mkrs=target.gain,
               OnOffScore=start.loss + target.gain,
               stringsAsFactors=FALSE)
}
