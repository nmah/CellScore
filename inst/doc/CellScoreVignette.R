## ----label='version', include=FALSE, cache=FALSE-------------------------
#suppressPackageStartupMessages(library(Biobase))
library(Biobase)
pkg.ver <- package.version("CellScore")

## ----label='Setup', include=FALSE, cache=FALSE----------------------

## Save the current working directory
dir.main <- getwd()
## Set the name of the directory in which figures will be saved (if any)
dir.figures <- 'figures'

## global chunk options
library(knitr)
opts_chunk$set(
    concordance=FALSE,
    cashe=2,
    ## cache is only valid with a specific version of R and session info
    ## cache will be kept for at most a month (re-compute the next month)
    cache.extra=list(R.version,
                     sessionInfo(),
                     format(Sys.Date(), '%Y-%m')
                     ),
    autodep=TRUE,
    fig.path=paste0(dir.figures,"/"),
    tidy=FALSE,
    size="small",
    message=FALSE,
    warning=FALSE
)
options(width=70, promp="R> ", continue="+  ", useFancyQuotes=FALSE)

## ----eval=FALSE-----------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite()

## ----eval=FALSE-----------------------------------------------------
#  biocLite(c("hgu133plus2CellScore", "CellScore"))

## ----eval=TRUE, echo=FALSE, cache=FALSE-----------------------------
options(BIOCINSTALLER_ONLINE_DCF=FALSE)
pc_info <- ifelse(grepl(getwd(),"kate"),
  "Intel Core i7-7500U CPU @ 2.70GHz",
  "Intel Core i7-5600U CPU @ 2.60GHz")

## ----eval=TRUE, echo=1:3, cache=FALSE, include=TRUE-----------------
library(Biobase)
library(CellScore)
library(hgu133plus2CellScore) # loads eset.std

## ----eval=TRUE------------------------------------------------------
## Locate the external data files in the CellScore package
rdata.path <- system.file("extdata", "eset48.RData",
                          package="CellScore")
tsvdata.path <- system.file("extdata", "cell_change_test.tsv",
                            package="CellScore")

if (file.exists(rdata.path) && file.exists(tsvdata.path)) {
    ## Load the normalized expressions of 48 test samples
    load(rdata.path)

    ## Import the cell change info for the loaded test samples
    cell.change <- read.delim(file=tsvdata.path, sep="\t",
                              header=TRUE, stringsAsFactors=FALSE)
    print("Content of cell.change")
    print(cell.change)

    ## Combine the standards and the test data
    eset <- combine(eset.std, eset48)
    print("dim(eset) returns")
    print(dim(eset))
}

## ----label='grouponoff', eval=TRUE, echo=FALSE----------------------
group.OnOff <- OnOff(eset, cell.change, out.put="marker.list")

## ----eval=TRUE------------------------------------------------------
group.OnOff <- OnOff(eset, cell.change, out.put="marker.list")
summary(group.OnOff)

## ----eval=TRUE------------------------------------------------------
head(group.OnOff$scores)

## ----eval=TRUE------------------------------------------------------
head(group.OnOff$markers)

## ----eval=TRUE, echo=TRUE, results='hide'---------------------------
individ.OnOff <- OnOff(eset, cell.change, out.put="individual")

## ----eval=TRUE------------------------------------------------------
save(file="OnOffscore.RData", individ.OnOff, group.OnOff)

## ----eval=TRUE------------------------------------------------------
sel.transition <- individ.OnOff$scores$start == "FIB" &
                  individ.OnOff$scores$target == "ESC"
sel.esc <- grepl("embryonic stem cell", individ.OnOff$scores$cell_type)
onoff.sel <- individ.OnOff$scores[sel.esc & sel.transition,
                                  c(4,6,7,12,13,19,20,21)]

## ----eval=TRUE------------------------------------------------------
summary(onoff.sel$OnOffScore)

## ----eval=TRUE------------------------------------------------------
onoff.sel[order(onoff.sel$OnOffScore)[1:4],]

## ----label='barplotonoffcode', eval=FALSE, echo=FALSE---------------
#  barplot.out <- BarplotOnOff(eset, group.OnOff$scores)

## ----label='barplotonoffdata', eval=TRUE, echo=TRUE, fig.keep='none'----
barplot.out <- BarplotOnOff(eset, group.OnOff$scores)
barplot.out

## ----pyramidbarplot, eval=TRUE, echo=FALSE, fig.pos ='!ht', fig.keep='high', fig.align='center', out.width='\\linewidth', fig.width=6, fig.height=4, fig.cap='\\textbf{Pyramid barplot of on/off scores.} The horizontal barplot shows the fraction of donor markers lost (blue) and the fraction of target markers gained (green). Each bar shows the results of one cell transition (left margin) for a particular derived cell type (right margin). The longer the coloured bars, the more successful the transition.'----
barplot.out <- BarplotOnOff(eset, group.OnOff$scores)

## ----eval=FALSE-----------------------------------------------------
#  pdf(file="GroupOnOffScore_Barplot.pdf")
#  barplot.out <- BarplotOnOff(eset, group.OnOff$scores)
#  dev.off()

## ----eval=TRUE------------------------------------------------------
tmp.time <- system.time(cs <- CosineSimScore(eset, cell.change,
                                             iqr.cutoff=0.1))
tmp.time

## ----eval=TRUE, results='hide'--------------------------------------
PlotCosineSimHeatmap(cs$cosine.general.groups, "general groups",
                     width=20, height=20, x=-20, y=3)

## ----eval=TRUE, results='hide'--------------------------------------
PlotCosineSimHeatmap(cs$cosine.subgroups, "subgroups",
                     width=14, height=14, x=-14, y=3)

## ----eval=FALSE-----------------------------------------------------
#  PlotCosineSimHeatmap(cs$cosine.samples, "samples",
#                       width=50, height=50, x=-50, y=10)

## ----label='HeatmapDataPrep', eval=TRUE, echo=FALSE, results='hide'----
## Get the names (IDs) of the sample and their description
samples.cs <- colnames(cs$cosine.samples)
samples.eset <- sampleNames(eset)

## Select the samples of interest and their corresponding score
sel.ips <- eset$category == "test" &
            eset$sub_cell_type1 %in% c("piPS-FIB", "iPS-FIB")
sel <- samples.cs %in% c(c("FIB", "ESC"), samples.eset[sel.ips])
cs.sel <- cs$cosine.samples[sel, sel]

## Rename columns/rownames to more descriptive labels
## as cs.sel is a symetric matrix, these are identical
ids <- match(colnames(cs.sel), samples.eset)
ids.na <- is.na(ids)
ids.rest <- na.omit(ids)
new.colnames <- c(colnames(cs.sel)[ids.na],
                  paste(eset$sub_cell_type1[ids.rest],
                        eset$sample_id[ids.rest],
                        sep="_")
                  )
colnames(cs.sel) <- rownames(cs.sel) <- new.colnames

## ----label='HeatmapCode', eval=TRUE, echo=FALSE, results='hide', fig.keep='none'----
## Plot the heatmap
PlotCosineSimHeatmap(cs.sel, "piPS", width=10, height=10, x=-10, y=3)

## ----eval=TRUE, echo=FALSE, results='hide'--------------------------
heatmap.filename <- "CosineSimilarityHeatmap_piPS.pdf"
heatmap.newpath <- file.path(dir.figures, heatmap.filename)
system(paste("mv", heatmap.filename, heatmap.newpath))

## ----eval=FALSE, echo=TRUE------------------------------------------
#  ## Get the names (IDs) of the sample and their description
#  samples.cs <- colnames(cs$cosine.samples)
#  samples.eset <- sampleNames(eset)
#  
#  ## Select the samples of interest and their corresponding score
#  sel.ips <- eset$category == "test" &
#              eset$sub_cell_type1 %in% c("piPS-FIB", "iPS-FIB")
#  sel <- samples.cs %in% c(c("FIB", "ESC"), samples.eset[sel.ips])
#  cs.sel <- cs$cosine.samples[sel, sel]
#  
#  ## Rename columns/rownames to more descriptive labels
#  ## as cs.sel is a symetric matrix, these are identical
#  ids <- match(colnames(cs.sel), samples.eset)
#  ids.na <- is.na(ids)
#  ids.rest <- na.omit(ids)
#  new.colnames <- c(colnames(cs.sel)[ids.na],
#                    paste(eset$sub_cell_type1[ids.rest],
#                          eset$sample_id[ids.rest],
#                          sep="_")
#                    )
#  colnames(cs.sel) <- rownames(cs.sel) <- new.colnames
#  
#  ## Plot the heatmap
#  PlotCosineSimHeatmap(cs.sel, "piPS", width=10, height=10, x=-10, y=3)

## ----eval=FALSE-----------------------------------------------------
#  PcaStandards(cs$pdataSub$experiment_id, "Experiment ID",
#               cs$esetSub.IQR)

## ----label='PcaCode', eval=FALSE------------------------------------
#  PcaStandards(cs$pdataSub$general_cell_type, "General Cell Type",
#               cs$esetSub.IQR)

## ----label='Pca', eval=TRUE, echo=FALSE, results='hide', fig.keep='high', out.width='1.3\\linewidth', fig.align='center', fig.pos='!ht', fig.width=14, fig.height=7, fig.cap='\\textbf{Principal component analysis of the reference dataset.} The reference samples cover a wide range of tissues and cell types from many studies and are the basis for comparison in the CellScore method. The analysis was applied on the expression matrix corresponidng to the most variable genes (as defined by the IQR cutoff). The plot, based on the first two prinicipal components, shows that similar cell types tend to cluster together. Brain tissue clusters on the right side, pluripotent and multipotent stem cells at the bottom, and somatic cell types on the left and upper clusters.'----
PcaStandards(cs$pdataSub$general_cell_type, "General Cell Type",
             cs$esetSub.IQR)

## ----eval=FALSE-----------------------------------------------------
#  pdf(file="StandardSamples_PCA_Labels.pdf", width=28, height=14)
#  PcaStandards(cs$pdataSub$general_cell_type, "General Cell Type",
#               cs$esetSub.IQR)
#  PcaStandards(cs$pdataSub$general_cell_type, "General Cell Type",
#               cs$esetSub.IQR,
#               text.label=cs$pdataSub$general_cell_type)
#  
#  dev.off()

## ----eval=TRUE, echo=TRUE-------------------------------------------
cellscore <- CellScore(eset, cell.change, individ.OnOff$scores,
                       cs$cosine.samples)

## ----eval=FALSE-----------------------------------------------------
#  save(file="VignetteResults.RData",
#       eset,                       # the combined expression dataset
#       group.OnOff, individ.OnOff, # the on/off score values
#       cs,                         # the cosine similaritiy values
#       cellscore                   # the CellScore values
#       )

## ----label='ScatterplotCode', eval=FALSE----------------------------
#  ScatterplotCellScoreComponents(cellscore, cell.change, FALSE)

## ----label='Scatterplot', eval=TRUE, echo=FALSE, out.width='\\linewidth', fig.pos='!ht', fig.align='center', fig.width=14, fig.height=7, fig.cap='\\textbf{Scatter plot of donor-like and target-like scores.} Derived cell types with the most successful transition have low donor-like score and high target-like scores and should cluster in the upper left-hand corner. In this example, the partially reprogrammed iPS cells (piPS-FIB) show gradual transition to their desired target cell type. A few iPS-FIB samples retain high donor-like scores, indicating unusual properties of these lines.'----
ScatterplotCellScoreComponents(cellscore, cell.change, FALSE)

## ----eval=FALSE-----------------------------------------------------
#  ScatterplotCellScoreComponents(cellscore, cell.change[2,], FALSE)

## ----label='Boxplot', eval=TRUE, fig.pos='!ht', out.width='\\linewidth', fig.align='center', fig.width=9, fig.height=6, fig.cap='\\textbf{Boxplot of CellScore values by subgroups.}'----
BoxplotCellScore(cellscore, cell.change)

## ----label='RugplotCode', eval=FALSE, echo=FALSE--------------------
#  secondary.property <- "transition_induction_method"
#  RugplotCellScore(cellscore, cell.change, secondary.property)

## ----label='Rugplot', eval=TRUE, echo=FALSE, out.width='1.1\\linewidth', fig.pos='!ht', fig.align='center', fig.width=14, fig.height=7, fig.cap='\\textbf{Rug plot of CellScore values.} The CellScore values are plotted for each derived cell type (left margin) within a study (GEO accession numbers in left margin). Each test sample is represented by a vertical line and colored by its transition induction method.'----
secondary.property <- "transition_induction_method"
RugplotCellScore(cellscore, cell.change, secondary.property)

## ----eval=FALSE-----------------------------------------------------
#  pdf(file="CellScoreReport_PerTransition.pdf", width=7, height=11)
#  CellScoreReport(cellscore, cell.change, group.OnOff$markers, eset)
#  dev.off()

## ----label='ReportDataPerp', include=FALSE, echo=FALSE--------------
# get the test cellscores from valid transitions
# defined by cell.change table
plot.data <- extractTransitions(cellscore, cell.change)
# Extract CellScores of test that should be plotted
# on the same page into list
plotgroup <- paste(plot.data$experiment_id,
                   plot.data$cxkey.subcelltype,
                   sep="_")
temp <- data.frame(plot.data, plotgroup, stringsAsFactors=FALSE)
# sort the table, show when the list is made, everything
# is already in the right order:
#  o by target
#  o by donor
#  o by study
ind <- order(temp$target, temp$donor_tissue, temp$experiment_id)
plot.data.ordered <- temp[ind,]
tg <- unique(paste(plot.data.ordered$experiment_id,
                   plot.data.ordered$cxkey.subcelltype,
                   sep="_")
             )
# get test data from plotgroup
test.data <- plot.data.ordered[plot.data.ordered$plotgroup %in% tg[1], ]

## ----label='ReportFig1', echo=FALSE, out.width='1.2\\linewidth', fig.keep='high', fig.width=14, fig.height=7, fig.pos='!ht', fig.align='center', fig.cap='\\textbf{Scatter plot of CellScore components.} The first plot of the CellScore report shows a scatter plot of the donor-like and target-like scores of the donor standard (in this case, fibroblasts (FIB); red) and target standard (embryonic stem cell (ESC); blue), as well as the derived cell types (induced pluripotent stem cells from fibroblasts (iPS-FIB); brown crosses). The number of samples from each group is indicated in parentheses in the figure legend.'----
mat <- matrix(c(1,1,2,2), nrow=1, ncol=4, byrow=TRUE)
layout(mat)
scatterplotDonorTargetTest(test.data, cellscore, FALSE)

## ----label='ReportFig2', echo=FALSE, out.width='\\linewidth', fig.width=8, fig.height=8, fig.keep='high', fig.pos='!ht', fig.align='center', fig.cap='\\textbf{Density and rug plots of CellScore values.} The CellScore values are shown as a density plot (upper panel) and as a rug plot (lower panel). The rug plot displays the test cells (iPS-FIB) in relation to the desired target cell type (ESC). Vertical dashed lines indicate the number of standard deviations away from the mean CellScore (bold vertical dashed line) of the target cell type.'----
mat <- matrix(c(1,1,1,2,2,2), nrow=2, ncol=3, byrow=TRUE)
layout(mat)
rugplotDonorTargetTest(test.data, cellscore)

## ----label='ReportFig3', echo=FALSE, out.width='\\linewidth', fig.width=9, fig.height=12, fig.keep='high', fig.pos='!ht', fig.align='center', fig.cap='\\textbf{Heatmap of donor and target markers.} Fibroblasts (red) are the donor cells, and embryonic stem cells are the desired target cells (blue). The test cells are induced pluripotent stem cells from fibroblasts (iPS-FIB; green).'----
mat <- matrix(c(1,1,1,1,1,1), nrow=2, ncol=3, byrow=TRUE)
layout(mat)
calls <- assayDataElement(eset, "calls")
rownames(calls) <- fData(eset)[, "probe_id"]
heatmapOnOffMarkers(test.data, group.OnOff$markers, pData(eset), calls)

## ----eval=TRUE, echo=TRUE-------------------------------------------
sessionInfo()

