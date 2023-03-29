# extractTransitions.R

# from CSReport_v1.4.3.R

#' Extract scores for given cell transitions
#'
#' This function extracts the values of the CellScore for all the test samples
#' of a given set of (valid) cell transition. While it can be used as
#' standalone, it serves as an internal function for several other CellScore
#' functions.
#' @param cellscore a data.frame of CellScore values as calculated by
#'   the function CellScore().
#' @param cell.change a data frame containing three columns, one for the
#'   start (donor) test and target cell type. Each row of the data.
#'   frame describes one transition from the start to a target cell type.
#' @return This function returns a data frame with the same columns as the input
#' data frame cellscore, extended with additional column that is used as a
#' single identifier of each valid cell transition. Technically, the output is
#'  subselection of the input data frame.
#' @keywords cellscore
#' @export
#' @seealso \code{\link[CellScore]{CellScore}} for details on CellScore
#'   calcualtion.
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
#'    ## Extract the scores for the transitions given in cell.change
#'    cellscore.cc <- extractTransitions(cellscore, cell.change)
#'
#'    ## View the sub_cell_type1 in the extracted object, it should be the same
#'    ## as the test cell types named in cell.change
#'    table(cellscore.cc$sub_cell_type1)
#' }

extractTransitions <- function(cellscore, cell.change) {

    ############################################################################
    ## PART 0. Check function arguments
    ############################################################################
    fun.main <- deparse(match.call()[[1]])
    .stopIfNotDataFrame(cellscore, 'cellscore', fun.main)
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)

    ############################################################################
    ## PART I. Get complex key for the transition table & for the score table
    ############################################################################
    key.cellchange <- .concatenateByRows(cell.change)
    key.subcelltype <- .concatenateByRows(data.frame(cellscore$start,
                                                     cellscore$sub_cell_type1,
                                                     cellscore$target))

    ############################################################################
    ## PART II. Extract only the test cell lines for valid cell transitions
    ############################################################################
    ## [LOOKS LIKE OVERKILL] TODO: some of these checkings should be done as a
    ## preprocesing step on the eset object and only after the object passes
    ## all the checks scores can be calculated and extracted (to omit the check
    ## in each function for NAs or invalid combinations)
    sel1 <- (!is.na(cellscore$start == cellscore$donor_tissue)) &
            cellscore$category == "test"
    sel2 <- !is.na(cellscore$cosine.donor)
    sel3 <- key.subcelltype %in% key.cellchange
    sel <- sel1 & sel2 & sel3

    ############################################################################
    ## PART III. Get subselection of the cellscore table & append the key
    ############################################################################
    data.frame(cellscore[sel, ],
               cxkey.subcelltype=key.subcelltype[sel],
               stringsAsFactors=FALSE)

}


