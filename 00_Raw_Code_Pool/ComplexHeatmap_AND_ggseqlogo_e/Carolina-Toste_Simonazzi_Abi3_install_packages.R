#| label: load libraries
#| include: false

here::i_am("Rscripts/install_packages.R")

#Install and load packages

## Create vector of used packages for CRAN 

pkg <- c("usethis", "janitor", "tidyverse", "reshape2", "DT",  "kableExtra", "labelled", 
         "downlit", "xml2", "devtools", "corrgram", "RColorBrewer", "viridis", "pheatmap", 
         "ggsci", "cowplot", "ggstatsplot", "ggvenn", "ggpubr","gridExtra", "tidytext", "Matrix",
         "ggthemes", "anndata", "qs", "clustree", "ggnewscale", "scCustomize", "gprofiler2",
         "FactoMineR", "factoextra", "rrcov", "bigutilsr", "remotes", "reticulate", "circlize", 
         "rstatix", "textshape", "hdf5r", "clustermole", "future", "msigdbr", "tidygraph", 
         "visNetwork", "naniar", "visdat", "assertr", "assertthat", "mashr", "digest", "scRNAseq", 
         "DuoClustering2018", "sctransform", "writexl", "circlize", "xgboost")

## Create vector of used packages for Bioconductor 

pkg_bioconductor <- c("cellity", "scater", "SingleCellExperiment",  
                      "DESeq2", "edgeR", "BiocParallel", "org.Mm.eg.db", "pathview", 
                      "gage", "gageData", "glmGamPoi", "DropletUtils", "ensembldb", 
                      "AnnotationHub", "patchwork", "scran", "scuttle", "PCAtools",
                      "batchelor", "bluster", "igraph", "harmony", "celldex", "AUCell",
                      "EnhancedVolcano", "dorothea", "decoupleR", "zellkonverter", 
                      "SummarizedExperiment", "ComplexHeatmap", "GenomeInfoDb", 
                      "biovizBase", "GENIE3", "RcisTarget", "zoo", "mixtools",
                      "NMF", "R2HTML", "Rtsne", "doMC", "doRNG", "TFBSTools", 
                      "universalmotif", "motifmatchr", "ggseqlogo", "jsonlite", 
                      "BSgenome.Mmusculus.UCSC.mm10", "GenomicRanges", "clusterProfiler",
                      "enrichplot", "variancePartition", "dreamlet", "WGCNA", "GeneOverlap",
                      "ggrepel", "UCell", "OmnipathR", "biomaRt", "rrvgo", "GSEABase", 
                      "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "lme4", "S4Vectors",
                      "HDF5Array", "ggrastr", "Signac", "SingleR", "GEOquery", "scry", "MAST", "scDblFinder", "decontX", "Nebulosa")




## Create install and load function

check_packages <- function(packages, bioconductor = FALSE) { 
  new_pkg <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  if (length(new_pkg) > 0) { 
    if (bioconductor) { 
      if (!requireNamespace("BiocManager", quietly = TRUE)) { 
        install.packages("BiocManager", repos = "https://cran.rstudio.com") 
      } 
      BiocManager::install(new_pkg, ask = FALSE) 
    } else { 
      install.packages(new_pkg, repos = "https://cran.rstudio.com") 
    } 
  }
  
  ## Load all packages
  
  lapply(packages, require, character.only = TRUE)
}

##Call function check_packages on vectors of CRAN and Bioconductor packages

check_packages(pkg) 
check_packages(pkg_bioconductor, bioconductor = TRUE)


##Create vector of packages from github


github_packages <- list(
  "satijalab/seurat-wrappers"      = "SeuratWrappers",
  "NightingaleHealth/ggforestplot" = "ggforestplot",
  "satijalab/azimuth"              = "Azimuth",
  "satijalab/seurat-data"          = "SeuratData",
  "cole-trapnell-lab/monocle3"     = "monocle3",
  "aertslab/SCopeLoomR"            = "SCopeLoomR",
  "aertslab/SCENIC"                = "SCENIC",
  "willtownes/quminorm"            = "quminorm",
  "satijalab/seurat"               = "Seurat",
  "igrabski/sc-SHC"                = "scSHC", 
  "crazyhottommy/scclusteval"      = "scclusteval",
  "lcrawlab/recall"                = "recall", 
  "immunogenomics/presto"          = "presto",
  "corceslab/CHOIR"                = "CHOIR",
  "ZJU-UoE-CCW-LAB/scCDC"          = "scCDC")

## Create seurat5 repos vector

seurat5_repos <- c(
  "satijalab/seurat",
  "satijalab/seurat-data",
  "satijalab/azimuth",
  "satijalab/seurat-wrappers")

##Update Seurat key dependencies

update.packages(oldPkgs = c("withr", "rlang"), ask = FALSE, quiet=TRUE)

##Use remotes to install github packages

for (repo in names(github_packages)) {
  pkg_name <- github_packages[[repo]]
  
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    if (repo %in% seurat5_repos) {
      remotes::install_github(repo, ref = "seurat5", quiet = TRUE)
    } else {
      remotes::install_github(repo)
    }
  }
  
  suppressPackageStartupMessages(
    library(pkg_name, character.only = TRUE)
  )
}


##Create vector of assertive packages

bitbucket_packages <- list(
  "assertive.properties",
  "assertive.base",
  "assertive.types",
  "assertive.strings",
  "assertive.datetimes",
  "assertive.data",
  "assertive.data.uk",
  "assertive.data.us",
  "assertive.code",
  "assertive.numbers",
  "assertive.files",
  "assertive.sets",
  "assertive.matrices",
  "assertive.models",
  "assertive.reflection",
  "assertive"
)

##Use remotes to install assertive packages

for (pkg in bitbucket_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    ## Construct repo string with username
    
    repo <- paste0("richierocks/", pkg)
    remotes::install_bitbucket(repo)
  }
  
  ## Load packages
  
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
  
}

## Remove vectors of packages and function from environment

rm(pkg, pkg_bioconductor,pkg_name, repo, seurat5_repos,  github_packages, bitbucket_packages, check_packages)


##Add README.md file at root

if (interactive()) {
  usethis::use_readme_md()
}
