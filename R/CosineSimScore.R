# CosineSimScore.R
#
# from version CosineScore_v1.4.1.R

#' Cosine similarity score
#'
#' This function calculates the cosine similarity for cell transitions.
#' @param eset an ExpressionSet containing data matrices of normalized
#'   expression data, present/absent calls, a gene annotation data frame and a
#'   phenotype data frame.
#' @param cell.change a data frame containing three columns, one for the
#'   start (donor) test and target cell type. Each row of the data
#'   frame describes one transition from the start to a target cell type.
#' @param iqr.cutoff set the threshold for top most variable genes which should
#'   be included for the cosine similarity calculation. Default is the top 10%
#'   genes, expressed as a fraction. All samples that are annotated as standards
#'   will be used for the iqr calculation.
#' @return This function returns a list of five objects, as follows:
#' \describe{
#'   \item{pdataSub}{the phenotype data frame describing the standard samples}
#'   \item{esetSub.IQR}{the expression value matrix, as filtered by IQR
#'   threshold}
#'   \item{cosine.general.groups}{a numeric matrix of cosine similarity
#'   between the centroids of all groups defined by eset@general_cell_types}
#'   \item{cosine.subgroups}{a numeric matrix of cosine similarity
#'   between the centroids of all gsubroups defined by eset@sub_cell_types1}
#'   \item{cosine.samples}{a numeric matrix of cosine similarity between
#'    general groups, subgroups and individual samples.}
#'}
#' @keywords cosine similarity
#' @seealso \code{\link[hgu133plus2CellScore]{hgu133plus2CellScore}} for details on the
#'   specific ExpressionSet object that shoud be provided as an input.
#' @export
#' @importClassesFrom Biobase ExpressionSet
#' @importMethodsFrom Biobase exprs pData
#' @importFrom lsa cosine
#' @importFrom stats IQR quantile median
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
#' }

CosineSimScore <- function(eset, cell.change, iqr.cutoff=0.1) {

    ###########################################################################
    ## PART 0. Check function arguments
    ###########################################################################
    fun.main <- deparse(match.call()[[1]])
    summarized_experiment <- .stopIfCantCoerceToSummarizedExperiment(eset, "eset", fun.main)
    .stopIfNotDataFrame(cell.change, 'cell.change', fun.main)
    .stopIfNotNumeric0to1(iqr.cutoff, 'min.diff.cutoff', fun.main)

    ## Preallocate output list
    result <- vector(mode="list", length=5)
    names(result) <- list("pdataSub", "esetSub.IQR", "cosine.general.groups",
                          "cosine.subgroups", "cosine.samples")

    ###########################################################################
    ## PART I. Filter samples according to phenoData
    ###########################################################################
    ## o phenoData table contains which samples should be used in the analysis
    ##    o the samples which should be used in the analysis will have an
    ##     assigned category eset@phenoData$category, as "standard" or "test"
    ##    o non-assigned samples with NA values will be ignored
    ## o NOTE
    ##   assigning NA values to samples is an easy way to eliminate samples
    ##   from the analysis, without having to remove them from all input tables
    ##   (eg removing from eset, pdata, calls)
    pdata <- colData(summarized_experiment)
    ## filter out samples with missing category and/or general cell type
    pdata.sel <- .filterPheno(pdata, fun.main, "na")

    ############################################################################
    ## PART II. Filter genes in dataset
    ############################################################################
    ## A. For the cosine scoring, keep only the standard reference samples as
    ## defined by the phenotable.
    selStd <- pdata.sel$category %in% "standard"
    if (sum(selStd) == 0)
        stop(paste("No standard reference samples found, exiting function",
             fun.main))
    result$pdataSub <- pdataStd <- pdata.sel[selStd, ]
    ynormStd <- assay(summarized_experiment[, selStd], 'exprs')

    ## B. Filter out not variable probes,
    ##    o variance in terms of IQR of group-wise median expressions
    ##    o want to get rid of probes that are below a given IQR threshold

    ## Calculate medians per group and put into table
    groupFac <- as.factor(pdataStd$general_cell_type)
    groupMed <- .calculateCentroids(ynormStd, "median", 1, groupFac)

    ## Calculate IQR of the group medians and filter by the given IQR threshold
    groupMedIQR <- apply(groupMed, 1, IQR)
    selGenes <- groupMedIQR >= quantile(groupMedIQR, probs=1 - iqr.cutoff)
    if (sum(selGenes) == 0) {
        ## this should almost never happen
        stop(paste("No gene passed the IQR-based filtering, exiting function",
             fun.main))
    }

    ## Final gene-filtered dataset of standards BASED ON TOP 10% OF IQR
    result$esetSub.IQR <- ynormStd[selGenes, ]

    ## We keep for score calculation the gene-filtered dataset with standards
    ## and test samples
    ynormIQR <- assay(summarized_experiment[selGenes, rownames(pdata.sel)], 'exprs')

    ############################################################################
    ## PART III. Calculating metrics
    ############################################################################
    ##  1. Euclidean distance between categories
    ##  2. Cosine similarity between categories
    ## -DONT USE-3. cosine "score" between categories
    ##      o 0 < cosine score <1 where range is set to 1-min(cosine.similarity)
    ## What is this good for?
    ##   a. Visualization and exploratory data analysis
    ##   b. Generates values to be used for cell scoring

    ## A. Calcualte Centroids across (general and sub-) cell types
    centroidExpList <-
        lapply(names(result[3:5]),
               function(x){
                   groupFac <-
                       switch(x,
                              "cosine.subgroups" =
                                  as.factor(pdataStd$sub_cell_type1),
                              "cosine.general.groups" =
                                  as.factor(pdataStd$general_cell_type),
                              "cosine.samples" = NULL
                       )
                   if (is.null(groupFac)){
                       ## All samples - no centroid calculation
                       ynormIQR
                   } else{
                       .calculateCentroids(result$esetSub.IQR, "mean", 1,
                                          groupFac)
                   }
               })
    names(centroidExpList) <- names(result[3:5])

    ## B. Calculate and save the Cosine similarities
    ##    Calculate only once and the subset by groups
    result[["cosine.samples"]] <- .getAllCosine(centroidExpList)
    result[3:4] <-
        lapply(names(result[3:4]),
               function(name){
                   .subsetSymetricMatrix(result[["cosine.samples"]],
                                         colnames(centroidExpList[[name]]))
               })
    result
}

## getAllCosine
##
## Local function that combines the given list of expression matrices and
## calcualted cosine similarity between all samples

.getAllCosine <- function(data.list){
    ## Combine all list elements in one big matrix
    data <- do.call("cbind", data.list)
    stopifnot(!is.null(data))

    ## Calculate the similarity; lsa::cosine()
    cosineSim <- cosine(data)
    ## Remove duplicated cell types (should not happen!)
    selected <- !duplicated(colnames(cosineSim))

    cosineSim[selected, selected]
}
