context("CellScore Constructor")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(Biobase)
BiocManager::install("hgu133plus2CellScore")
library(hgu133plus2CellScore)

test_that("We can create a CellScore object from data", {
                 ## Locate the external data files in the CellScore package
     rdata.path <- system.file("extdata", "eset48.RData", package = "CellScore")
     tsvdata.path <- system.file("extdata", "cell_change_test.tsv",
                                  package = "CellScore")
     
     if (file.exists(rdata.path) && file.exists(tsvdata.path)) {
     
        ## Load the expression set with normalized expressions of 48 test samples
        load(rdata.path)
     
        ## Import the cell change info for the loaded test samples
        cell.change <- read.delim(file= tsvdata.path, sep="\t",
                                  header=TRUE, stringsAsFactors=FALSE)
        ## Combine the standards and the test data
        eset <- combine(eset.std, eset48)
     
        ## Generate cosine similarity for the combined data
        ## NOTE: May take 1-2 minutes on the full eset object
        ## so we subset it for 4 cell types
        pdata <- pData(eset)
        sel.samples <- pdata$general_cell_type %in% c("ESC", "EC", "FIB", "KER")
        eset.sub <- eset[, sel.samples]
        cs <- CosineSimScore(eset.sub, cell.change, iqr.cutoff=0.1)
     
        ## Generate the on/off scores for the combined data
        individ.OnOff <- OnOff(eset.sub, cell.change, out.put="individual")
     
        ## Generate the CellScore values for all samples
        cellscore <- CellScore(eset.sub, cell.change, individ.OnOff$scores,
                               cs$cosine.samples)
        expect_named(cellscore, c('composite.ID', 'experiment_id', 'sample_id', 'platform_id', 'cell_type', 'disease_status', 'category', 'general_cell_type', 'donor_tissue', 'sub_cell_type1', 'transition_induction_method', 'donor_cell_body_location', 'start', 'target', 'markers.start', 'markers.target', 'start.mkrs.in.test', 'target.mkrs.in.test', 'loss.start.mkrs', 'gain.target.mkrs', 'OnOffScore', 'fraction.target', 'cosine.target', 'fraction.donor', 'cosine.donor', 'target.like', 'donor.like', 'CellScore', 'index'))
        expect_length(cellscore$platform_id, 430)
        expect_length(cellscore$CellScore, 430)
     }
        
})
