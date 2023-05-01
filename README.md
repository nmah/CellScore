# CellScore

## Local Development

To install CellScore locally from this repo for development or testing, you need to have a recent R (4 or later) installed.

Then clone the repo:

```sh
git clone https://github.com/flaviusb/CellScore.git
```

Then, run R inside the cloned directory. In order to load the package locally, you can run:

```R
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("SummarizedExperiment")
BiocManager::install("getDEE2")
devtools::load_all()
```

