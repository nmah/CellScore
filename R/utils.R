# utils.R

#' Helper functions called from the exported package functions

################################################################################
## FUNCTION: .zscore
################################################################################
#' Calculates z-score. We use z-score as a measure of how close the standards
#' cluster to each other
#' @param x a vector of values from sample population
#' @param y value of test sample
#' @return z-score
#' @keywords z-score
#' @importFrom stats sd
.zscore <- function(x, y) {
    # sd.pop <- sd(x)*sqrt((length(x)-1)/(length(x))) # population sd
    sd.sample <- sd(x) # sample sd
    mean.samples <- mean(x)
    #(y-mean.samples)/sd.pop
    (y - mean.samples)/sd.sample
}

################################################################################
## FUNCTION: .calculateCentroids
################################################################################
#' Calculates group centroids, using either median or mean aggregation across
#' rows or coulmns of a given numeric matrix.
#' @param es a numeric matrix, e.g. an expression matrix,
#' @param aggregation a character that codes for the aggregation function
#' name ("mean" or "median")
#' @param dimension an integer that identifies the dimension
#'   along which aggregation is applied (1: rows or 2 : columns)
#' @param group a factor that keeps the group labels
#' @keywords aggregation centroids
.calculateCentroids <- function(es, aggregation = c('mean', 'median'),
                                dimension = c(1, 2), group){
    ## Check the group argument has appropriate size
    switch(as.character(dimension),
           "1" = stopifnot(length(group) == ncol(es)),
           "2" = stopifnot(length(group) == nrow(es)))
    ## Aggreagtion
    temp <- apply(es, dimension,
                  function(x){
                      tapply(x, group,
                             function(y){
                                 do.call(aggregation, list(y))
                             })
                  })
    temp <- t(temp)

    ## Check the order of elements along the unmodified dimension
    switch(as.character(dimension),
           "1" = stopifnot(identical(rownames(temp), rownames(es))),
           "2" = stopifnot(identical(colnames(temp), colnames(es))))
    temp
}

################################################################################
## FUNCTION: .concatenateByRows
################################################################################
#' Collapses the the values of each row of a given data frame in one string.
#' @param df a data.frame
#' @return  a character vector
#' @keywords concatenate
.concatenateByRows <- function(df){
    apply(df, 1, function(x) paste(x, collapse="_"))
}

################################################################################
## FUNCTION: .subsetSymetricMatrix
################################################################################
#' Collapses the the values of each row of a given data frame in one string.
#' @param data a matrix
#' @param which.samples a character vector with subset of (row) column names
#' @return  a submatrix
#' @keywords submatrix
.subsetSymetricMatrix <- function(data, which.samples){
    stopifnot(isSymmetric(data) &&
                  identical(rownames(data), colnames(data)) &&
                  !is.null(which.samples) &&
                  sum(selected <- rownames(data) %in% which.samples))

    data[selected, selected]
}


################################################################################
## FUNCTION: .filterPheno
################################################################################
#' Filters samples from phenotype data with missing info
#' @param pheno a data.frame with phenotype descriptions
#' @param calling.fun a character - name of the parent function calling this
#'   function
#' @param flag a character - type of filter, one of the following: 'na',
#'   'group', subgroup', anygroup' (the last option checks for any match in
#'   group or subgroup)
#' @param flag.value a character - cell name needed for the 'specifc' filter
#' @return filtered data.frame
#' @keywords filter
.filterPheno <- function(pheno, calling.fun,
                         flag = c("na", "group", "subgroup", "anygroup"),
                         flag.value=NULL){
    ## If filtering by (sub)groups, the flag.value must be specified
    stopifnot(flag == "na" || !is.null(flag.value))

    sel <- switch(flag,
                  na = !is.na(pheno$category) & !is.na(pheno$general_cell_type),
                  group = pheno$general_cell_type %in% flag.value,
                  subgroup = pheno$sub_cell_type1 %in% flag.value,
                  anygroup = pheno$general_cell_type %in% flag.value |
                      pheno$sub_cell_type1 %in% flag.value)

    if (sum(sel) == 0){
        print(sel)
        stop(paste("No samples left after phenotype filtering, exiting function",
                   calling.fun))
    }

    pheno[sel, ]
}

################################################################################
## FUNCTION: .filterTransitions
################################################################################
#' Selects cell transitions that have donor and target samples defined in pheno
#' or that at have at least 3 samples each
#' @param cell.change a data.frame with selected transitions
#' @param pheno a data.frame with phenotype descriptions
#' @param calling.fun a character - name of the parent function calling this
#'   function
#' @param flag a character - type of filter, either 'sample.counts' or
#'   'valid.names'
#' @return filtered data.frame
#' @keywords filter
.filterTransitions <- function(cell.change, pheno, calling.fun,
                               flag=c("valid.names", "sample.counts")){
    sel <- switch(flag,
                  valid.names = {
                      sel.start <- cell.change$start %in% pheno$general_cell_type
                      sel.target <- cell.change$target %in% pheno$general_cell_type
                      sel.start & sel.target
                  },
                  sample.counts = {
                      sample.counts <- table(pheno$general_cell_type)
                      sample.counts <- names(sample.counts[sample.counts >= 3])

                      sel.start <- cell.change$start %in% sample.counts
                      sel.target <- cell.change$target %in% sample.counts
                      sel.start & sel.target
                  })

    if (sum(sel) == 0){
        stop(paste("No cell transitions left after filtering, exiting function",
                   calling.fun))
    }

    cell.change[sel,]

}

################################################################################
## FUNCTION: .calculateFractions
################################################################################
#' Applies subseting by columns and aggregation by rows
#' on the input data.frame/matrix of binary values (such as the calls matrix)
#' @param df a data.frame
#' @return  a numeric vector
#' @keywords fractions
.calculateFractions <- function(df, list.names=NULL){
    ## Filter df by columns using the provided set of column names
    ## Keep the data.frame format even when single column left
    if (!is.null(list.names)){
        df <- df[, colnames(df) %in% list.names, drop=FALSE]
    }
    ## Calulculate fractions by rows, which for binary matrix is the row mean
    switch(as.character(ncol(df) > 1),
           "TRUE" = rowMeans(df, na.rm=TRUE),
           "FALSE" = df[, 1])
}

################################################################################
## FUNCTION: .getMainColours
################################################################################
#' To get the colours for each group in a given tarnsition, used for plotting
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
.getMainColours <- function(group=c("all", "donor", "target", "test"),
                            colramp=FALSE){
    col.list <- list(donor=brewer.pal(9,"Reds"),
                     target=brewer.pal(9,"Blues"),
                     test="limegreen")

    col.list[1:2] <- lapply(col.list[1:2],
                            function(x){
                                if (colramp){
                                    colorRampPalette(x[-(1:3)])
                                }else{
                                    x[6]
                                }
                            })

    if (group == "all") return(col.list)
    col.list[[group]]
}

################################################################################
## FUNCTION: .colourMapping
################################################################################
#' Returns a data frame with colour mapping for each unique value in the input
#' charcter vector data.in. Colours can be provided as second argument. This
#' colour vector will be recycled if it is shorter than the set of unique values
#' in data.in.
#' @param data.in a character vector
#' @param col.palette a character vector of colours
#' @return a data.frame
#' @keywords colours
.colourMapping <- function(data.in, col.palette=NULL){
    size.groups <- length(unique(data.in))

    if (is.null(col.palette)){
        col.palette <- c("saddlebrown","limegreen","orange", "blue","magenta",
                         "red", "cyan","yellow","lightslateblue","midnightblue",
                         "deeppink4","purple", "cornflowerblue","black",
                         "paleturquoise","steelblue", "chartreuse4", "khaki",
                         "green","salmon","pink")
    }
    size.palette <- length(col.palette)

    ## Colours will be recycled if there are more than 21 groups
    if (size.groups > size.palette) {
        col.palette <- rep(col.palette, ceiling(size.groups/size.palette))
    }

    ## Collect colour mappings in one table
    df <- data.frame(group=data.in,
                     col=col.palette[as.numeric(factor(data.in))],
                     stringsAsFactors=FALSE)

}

################################################################################
## FUNCTION: .stopIfNotExpressionSet
################################################################################
#' @importFrom methods is
.stopIfNotExpressionSet <- function(x, x.name, fun.name){
    if (!is(x,"ExpressionSet")) {
        stop(paste("In the function", fun.name, "the", x.name,
                   "argument shoud be an ExpressionSet."))
    }
}

################################################################################
## FUNCTION: .stopIfNotExpressionSet
################################################################################
#' @importFrom methods is
.stopIfCantCoerceToSummarizedExperiment <- function(x, x.name, fun.name){
    if (is(x, "SummarizedExperiment")) {
           return(x)
    }
    if (is(x,"ExpressionSet")) {
      # Attempt to coerce to a SummarizedExperiment
      return(as(x, "RangedSummarizedExperiment"))
    }
    stop(paste("In the function", fun.name, "the", x.name,
               "argument shoud be a SummarizedExperiment or an ExpressionSet."))
}

################################################################################
## FUNCTION: .stopIfNotDataFrame
################################################################################
.stopIfNotDataFrame <- function(x, x.name, fun.name){
    if (!is.data.frame(x)) {
        stop(paste("In the function", fun.name, "the", x.name,
                   "argument shoud be a data.frame."))
    }
}

################################################################################
## FUNCTION: .stopIfNotNumeric0to1
################################################################################
.stopIfNotNumeric0to1 <- function(x, x.name, fun.name) {
    if (!is.numeric(x) || x >= 1 || x <= 0) {
        stop(paste("In function", fun.name, "the", x.name,
                   "argument shoud be a real number in (0,1]."))
    }
}

################################################################################
## FUNCTION: .stopIfNotSymetricMatrix0to1
################################################################################
.stopIfNotSymetricMatrix0to1 <- function(x, x.name, fun.name){
    if (!isSymmetric(as.matrix(x)) ||
        !identical(rownames(x), colnames(x)) ||
        all(x < 0) || all(x > 1)) {
        stop(paste("In function", fun.name, "the", x.name,
                   "argument shoud be a symetric numeric matrix with values",
                   "in [0,1]."))
    }

}

################################################################################
## FUNCTION: .stopIfNotBinaryMatrix
################################################################################
.stopIfNotBinaryMatrix <- function(x, x.name, fun.name){
    if (!all(x == matrix(as.numeric(as.logical(x)), nrow=nrow(x)))){
        stop(paste("In function", fun.name, "the", x.name,
                   "argument shoud be a binary numeric matrix."))
    }
}
