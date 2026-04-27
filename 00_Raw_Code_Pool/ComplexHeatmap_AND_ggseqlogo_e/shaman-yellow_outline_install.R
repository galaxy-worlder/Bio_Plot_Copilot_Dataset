
install.packages(c("usethis", "devtools", "officer", "gt"))

# install.packages(c("CliquePercolation", "qgraph", "Matrix"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BiocStyle", "ChemmineOB"))

BiocManager::install(c("rvest"))

# install.packages(c("igraph", "ggraph"))

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Cao-lab-zcmu/MCnebula2")

BiocManager::install(c("FELLA", "xcms", "ggtree"))
## For tools query chemical classification via ClassyFire API
remotes::install_github('aberHRML/classyfireR')
## For tools convert CID to KEGG ID
remotes::install_github('xia-lab/MetaboAnalystR')

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova",
    "Rgraphviz", "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph",
    "siggenes","BiocParallel", "MSnbase",
    "multtest","RBGL","edgeR","fgsea","devtools","crmn","httr","qs")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    message(c(new_pkgs, " packages added..."))
  }
  if((length(new_pkgs)<1)){
    message("No new packages added...")
  }
}

metanr_packages()

install.packages(c("pdftools", "officedown"))

BiocManager::install(c("limma", "biomaRt", "STRINGdb"))
remotes::install_github('YuLab-SMU/clusterProfiler')

install.packages("RSelenium")

install.packages(c("car", "rms"))
BiocManager::install("ropls")

install.packages(c("randomForest"))

install.packages(c("pROC"))

install.packages(c("DescTools"))

# install.packages(c("gmodels"))

BiocManager::install(c("meta"))
BiocManager::install(c("GEOquery"))

# BiocManager::install("GenomicAlignments")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

install.packages("ggVennDiagram")

BiocManager::install("RNAseq123")

BiocManager::install("WGCNA")

install.packages(c("simputation", "naniar", "agricolae"))

install.packages(c("pander"))

install.packages(c("survival", "survminer"))

## <https://satijalab.org/seurat/>
## <https://satijalab.org/seurat/articles/install.html>

BiocManager::install("SparseArray")
install.packages(c("sva"))
BiocManager::install(c('fastDummies', 'RcppHNSW', 'RSpectra'))
remotes::install_github("satijalab/seurat", "seurat5")

BiocManager::install(c('HDF5Array', 'beachmat'))
BiocManager::install("glmGamPoi")

reticulate::install_miniconda()
reticulate::py_install(packages = 'umap-learn')

BiocManager::install(c("SingleR"))
BiocManager::install(c("celldex"))

BiocManager::install(c('NMF', 'circlize', 'ComplexHeatmap', 'BiocNeighbors'))
remotes::install_github("sqjin/CellChat")

remotes::install_github('satijalab/seurat-wrappers')

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
    'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
    'SummarizedExperiment', 'batchelor', 'HDF5Array',
    'terra', 'ggrastr'))
BiocManager::install("rsample")
remotes::install_github('cole-trapnell-lab/monocle3')

BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
BiocManager::install("rly")
remotes::install_github("cole-trapnell-lab/garnett", ref = "monocle3")

## 
install.packages("ggupset")
install.packages("UpSetR")

## remotes::install_github("yixianfan/CRDscore")
## BiocManager::install("Rsamtools")

BiocManager::install("TCGAbiolinks")

install.packages(c("grImport"))
install.packages(c("EFS"))
install.packages(c("lazyWeave"))

system("pip3 install umap-learn")

install.packages("spiralize")

## https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html
## https://bioconductor.org/packages/release/bioc/vignettes/MicrobiotaProcess/inst/doc/MicrobiotaProcess.html
BiocManager::install("biomformat")
# BiocManager::install("MicrobiotaProcess")
remotes::install_github("YuLab-SMU/MicrobiotaProcess")
BiocManager::install("qiime2R")
install.packages(c("gghalves", "ggh4x", "corrr", "ggside"))

install.packages("gsalib")

install.packages("vcfR")
BiocManager::install("maftools")

BiocManager::install("org.Ss.eg.db")
# install.packages("queryup")
BiocManager::install("UniProt.ws")


# BiocManager::install('biodbUniprot')

install.packages("r3dmol")

install.packages("ggpval")

BiocManager::install("pathview")

## Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression
## [@EstimatingTheBecht2016]
## https://github.com/ebecht/MCPcounter

## xCell: digitally portraying the tissue cellular heterogeneity landscape
## https://github.com/dviraran/xCell

## MitoCarta3.0（http://www.broadinstitute.org/mitocarta）

## ST
BiocManager::install("hdf5r")

devtools::install_github("liuhong-jia/scAnno")
BiocManager::install("infercnv")

devtools::install_github("navinlabcode/copykat")

## GWAS: GAPIT <https://github.com/jiabowang/GAPIT>
## GWAS: PoPs <https://www.nature.com/articles/s41588-023-01443-6#code-availability>

## Transcriptomic deconvolution in cancer and other heterogeneous tissues remains challenging. Available methods lack the ability to estimate both component-specific proportions and expression profiles for individual samples
## https://wwylab.github.io/DeMixT/
## https://github.com/wwylab/TmS

## https://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html

## Robust integration of multiple single-cell RNA sequencing datasets using a single reference space
## https://github.com/bioinfoDZ/RISC
## https://github.com/bioinfoDZ/RISC/blob/master/GSE123813_Vignette_RISC_v1.6.pdf
remotes::install_github("cvarrichio/Matrix.utils")
devtools::install_github("https://github.com/bioinfoDZ/RISC.git")

BiocManager::install("DropletUtils")

## De novo detection of somatic mutations in high-throughput single-cell profiling data sets
## [@DeNovoDetectiMuyas2023]
## https://github.com/cortes-ciriano-lab/SComatic

## PhyloVelo - Phylogeny-based transcriptomic velocity of single cells
## https://github.com/kunwang34/PhyloVelo

install.packages("filesstrings")

## https://phylovelo.readthedocs.io/en/latest/notebook/cd8%2BTcells-Copy1.html

## https://github.com/velocyto-team/velocyto.R
remotes::install_github("velocyto-team/velocyto.R")

devtools::install_github("Hy4m/linkET")

# https://github.com/irrationone/cellassign

install.packages("bibliometrix")

# https://cran.r-project.org/web/packages/PantaRhei/vignettes/panta-rhei.html
install.packages("lubridate")
install.packages("openxlsx2")
# https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html

install.packages("randomForestSRC")
install.packages("ggthemes")

# install.packages("ecce")
install.packages("gtranslate")

remotes::install_github("omarwagih/ggseqlogo")

install.packages("showtext")

options(timeout = 100000)
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("GenomicFiles")

BiocManager::install("MungeSumstats")

install.packages("PubChemR")
# BiocManager::install("GEOmetadb")
install.packages("rscopus")
install.packages("foreign")

BiocManager::install("xcms")
# sudo apt install libgsl-dev
BiocManager::install("curatedMetagenomicData")
install.packages("readxlsb")

install.packages("tesseract")
