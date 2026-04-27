#'Computation of gene expression signature scores.
#'
#'Computes gene expression signature scores. Also draws boxplots
#'representing the average signature scores for each subtype.
#'
#'@param expr A \code{SummarizedExperiment} object or a matrix containig gene
#'expression profiles. If input is a \code{SummarizedExperiment}, the first
#'element in the assays list should be a matrix of gene expression.
#'Rows and columns of the gene expression matrix correspond to genes and
#'samples, respectively (rownames must be to gene symbols).
#'@param pred A vector of predicted consensus molecular subtypes.
#'@param rnaseq logical to determine if input data is
#'RNA-Seq gene expression profile. By default, it is FALSE.
#'@return A matrix of gene expression signature scores.
#'@details
#'
#'\code{computeGES} calculates the following 7 gene expression
#'signature scores:
#'\itemize{
#'\item EMT (epithelial-mesenchymal transition):
#'average of expression values of genes included in the
#'EMT signature published by \cite{Tan et al. (2014)}.
#'\item Stromal: stromal score representing the presence of
#'stromal cells in tumor tissues (computed using the ESTIMATE algorithm).
#'\item Immune: immune score representing the presence of
#'immune cells in tumor tissues (computed using the ESTIMATE algorithm).
#'\item Microenvironment: microenvironment score representing the sum
#'of all immune and stromal cell types (computed using xCell)
#'\item Stemness: stemness index computed using the method developed by
#'\cite{Malta et al. (2018)}.
#'\item Hormone: average of expression values of AR, ERBB2, ESR1, and PGR.
#'\item CIN (chromosomal instability): average of expression values of genes
#'included in the CIN70 signature published by \cite{Carter et al. (2006)}.
#'}
#'
#'@references Aran, D. et al. (2017). xCell: digitally portraying
#'the tissue cellular heterogeneity landscape. \emph{Genome biology}, 18, 220.
#'@references Carter, S.L. et al. (2006). A signature of chromosomal
#'instability inferred from gene expression profiles predicts
#'clinical outcome in multiple human cancers. \emph{Nature genetics}, 38, 1043.
#'@references Malta, T.M. et al. (2018). Machine learning
#'identifies stemness features associated with
#'oncogenic dedifferentiation. \emph{Cell}, 173, 338-354.
#'@references Tan, T.Z. et al. (2014). Epithelial-mesenchymal
#'transition spectrum quantification and its efficacy in deciphering
#'survival and drug responses of cancer
#'patients. \emph{EMBO molecular medicine}, 6, 1279-93.
#'@references Yoshihara, K. et al. (2013). Inferring
#'tumour purity and stromal and immune cell admixture
#'from expression data. \emph{Nature communications}, 4, 2612.
#'
#'@importFrom pheatmap pheatmap
#'@importFrom grDevices colorRampPalette
#'@importFrom RColorBrewer brewer.pal
#'@importFrom stats cor wilcox.test
#'@importFrom ggpubr ggboxplot
#'@importFrom methods is
#'@import ggplot2
#'@import grid
#'@import SummarizedExperiment
#'@export
#'@examples
#'# Load gene expression profiles of TNBC samples
#'data(GSE25055)
#'
#'# Predict consensus molecular subtypes of TNBC samples
#'prediction <- predictCMS(expr = GSE25055)
#'
#'# Compute gene expression signature scores
#'resultGES <- computeGES(expr = GSE25055, pred = prediction, rnaseq = FALSE)
computeGES <- function(expr, pred, rnaseq = FALSE){

  if (is(expr, "SummarizedExperiment")){
    exp.mat <- assays(expr)[[1]]
  } else{
    exp.mat <- expr
  }

  CMS.palette <- c("MSL" = "brown2", "IM" = "gold2",
                   "LAR" = "yellowgreen", "SL" = "midnightblue")

  pred <- factor(pred, levels = c("MSL", "IM", "LAR", "SL"))
  EMT.score <- colMeans(exp.mat[rownames(exp.mat) %in% EMT.geneset,])
  CIN.score <- colMeans(exp.mat[rownames(exp.mat) %in% CIN.geneset,])
  hormone.score <- colMeans(exp.mat[rownames(exp.mat) %in%
                                      c("AR", "ESR1", "ERBB2", "PGR"),])
  X <- exp.mat[rownames(exp.mat) %in% names(Stemness.geneset),]
  Stemness.geneset = Stemness.geneset[rownames(X)]
  s <- apply(X, 2, function(z) {cor(z, Stemness.geneset,
                                    method="sp", use="complete.obs")})
  s <- s - min(s)
  stemness.score <- s / max(s)
  ESTIMATE.score <- computeESTIMATEscore(exp.mat)
  stromal.score <- unlist(ESTIMATE.score["StromalSignature",])
  immune.score <- unlist(ESTIMATE.score["ImmuneSignature",])
  microenvironment.score <- computexCellScore(exp.mat, rnaseq)
  sig.mat <- rbind(EMT.score, stromal.score, immune.score,
                   microenvironment.score, stemness.score,
                   hormone.score, CIN.score)
  rownames(sig.mat) <- c("EMT", "Stromal", "Immune", "Microenvironment",
                         "Stemness", "Hormone", "CIN")

  pval1 <- pval2 <- pval3 <- pval4 <- NA
  if(length(unique(pred)) > 1){

    if("MSL" %in% as.character(pred)){
      pval1 <- wilcox.test(stromal.score ~
                             ifelse(pred == "MSL", 1, 0))$p.value
    }
    if("IM" %in% as.character(pred)){
      pval2 <- wilcox.test(immune.score ~
                             ifelse(pred == "IM", 1, 0))$p.value
    }
    if("LAR" %in% as.character(pred)){
      pval3 <- wilcox.test(hormone.score ~
                             ifelse(pred == "LAR", 1, 0))$p.value
    }
    if("SL" %in% as.character(pred)){
      pval4 <- wilcox.test(stemness.score ~
                             ifelse(pred == "SL", 1, 0))$p.value
    }

  }

  sub1 <- sub2 <- sub3 <- sub4 <- ""

  if(!is.na(pval1)){
    sub1 <- bquote(paste("Wilcoxon (MSL vs. others) ",
                         italic(p), " = ", .(format(pval1, digits = 2))))
  }

  if(!is.na(pval2)){
    sub2 <- bquote(paste("Wilcoxon (IM vs. others) ",
                         italic(p), " = ", .(format(pval2, digits = 2))))
  }

  if(!is.na(pval3)){
    sub3 <- bquote(paste("Wilcoxon (LAR vs. others) ",
                         italic(p), " = ", .(format(pval3, digits = 2))))
  }

  if(!is.na(pval4)){
    sub4 <- bquote(paste("Wilcoxon (SL vs. others) ",
                         italic(p), " = ", .(format(pval4, digits = 2))))
  }

  TITLE_SIZE <- 12
  SUBTITLE_SIZE <- 10

  sigdat <- data.frame(CMS = pred, Stromal = stromal.score,
                       Immune = immune.score, Hormone = hormone.score,
                       Stem = stemness.score)

  p1 <- ggboxplot(sigdat, x = "CMS", y = "Stromal",
                  fill = "CMS", palette = CMS.palette) +
    theme_bw() + labs(title = "Stromal", subtitle = sub1) +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = TITLE_SIZE,
                                    hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size=  SUBTITLE_SIZE, hjust = 0.5))

  p2 <- ggboxplot(sigdat, x = "CMS", y = "Immune",
                  fill = "CMS", palette = CMS.palette) +
    theme_bw() + labs(title = "Immune", subtitle = sub2) +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = TITLE_SIZE,
                                    hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size=  SUBTITLE_SIZE, hjust = 0.5))

  p3 <- ggboxplot(sigdat, x = "CMS", y = "Hormone",
                  fill = "CMS", palette = CMS.palette) +
    theme_bw() + labs(title = "Hormone", subtitle = sub3) +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = TITLE_SIZE,
                                    hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size=  SUBTITLE_SIZE, hjust = 0.5))

  p4 <- ggboxplot(sigdat, x = "CMS", y = "Stem",
                  fill = "CMS", palette = CMS.palette) +
    theme_bw() + labs(title = "Stem", subtitle = sub4) +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = TITLE_SIZE,
                                    hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(size=  SUBTITLE_SIZE, hjust = 0.5))

  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(p3)
  g4 <- ggplotGrob(p4)
  g <- cbind(rbind(g1, g3, size = "first"),
             rbind(g2,  g4, size = "first"), size = "first")
  grid.newpage()
  grid.draw(g)

  return(sig.mat)

}

#'Computation of microenvironment score
#'
#'Computes a microenvironment score. This function wraps around
#'the \code{xCellAnalysis} function of the \code{xCell} package to
#'compute a microenvironment score.
#'
#'@param mat A matrix of gene expression with genes in rows and
#'samples in columns (rownames corresopnding to gene symbol).
#'@param rnaseq logical to determine if input data is
#'RNA-Seq gene expression profile
#'@return A data frame containing stromal and immune scores
#'@importFrom GSVA gsva
#'@importFrom pracma lsqlincon
#'@importFrom stats aggregate
#'@import quadprog
#'@keywords internal
computexCellScore <- function(mat, rnaseq){

  signatures <- xCell.data$signatures
  genes <- xCell.data$genes
  if(rnaseq){
    spill <- xCell.data$spill
  } else{
    spill <- xCell.data$spill.array
  }

  #Compute enrichment scores
  shared.genes <- intersect(rownames(mat), genes)
  expr <- mat[shared.genes,]
  expr <- apply(expr, 2, rank)
  scores <- gsva(expr, signatures, method = "ssgsea", ssgsea.norm = FALSE,
                 parallel.sz = 1, verbose = FALSE)
  scores <- scores - apply(scores, 1, min)
  cell.types <- unlist(strsplit(rownames(scores), "%"))
  cell.types <- cell.types[seq(1, length(cell.types), 3)]
  agg <- aggregate(scores ~ cell.types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  #Transform scores
  rows <- rownames(scores)[rownames(scores) %in% rownames(spill$fv)]
  tscores <- scores[rows, ]
  minX <- apply(tscores, 1, min)
  A <- rownames(tscores)
  tscores <- (as.matrix(tscores) - minX)/5000
  tscores[tscores < 0] <- 0
  tscores <- (tscores^spill$fv[A,2])/(spill$fv[A,3]*2)

  #Adjust scores
  K <- spill$K * 0.5
  diag(K) <- 1
  rows <- rownames(tscores)[rownames(tscores) %in%
                              rownames(K)]
  ascores <- apply(tscores[rows, ], 2, function(x) lsqlincon(K[rows,rows],
                                                             x, lb = 0))
  ascores[ascores<0] <- 0
  rownames(ascores) <- rows

  #Compute microenvironment scores
  immune.score <- apply(ascores[c('B-cells','CD4+ T-cells',
                                  'CD8+ T-cells','DC','Eosinophils',
                                  'Macrophages','Monocytes','Mast cells',
                                  'Neutrophils','NK cells'),],2,sum)/1.5
  stromal.score <- apply(ascores[c('Adipocytes','Endothelial cells',
                                   'Fibroblasts'),],2,sum)/2
  microenvironment.score <- immune.score + stromal.score

  microenvironment.score
}

#'Computation of stromal and immune scores
#'
#'Computes stromal and immune scores. This function was borrowed from
#'the \code{estimate} package and changed to accept R object as input.
#'
#'@param mat A matrix of gene expression with genes in rows and
#'samples in columns (rownames corresopnding to gene symbol).
#'@return A data frame containing stromal and immune scores
#'@keywords internal
computeESTIMATEscore <- function(mat){

  m <- mat
  gene.names <- rownames(mat)
  sample.names <- colnames(mat)
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(SI.geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI.geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES,
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  score.data

}
