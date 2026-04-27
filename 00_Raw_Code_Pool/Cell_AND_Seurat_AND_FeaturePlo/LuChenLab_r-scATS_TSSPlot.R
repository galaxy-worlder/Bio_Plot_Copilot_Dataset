#' @title TSSPlot
#' @name TSSPlot
#' @description Find markers of given group
#' @importFrom SummarizedExperiment colData
#' @importFrom reshape2 melt
#' @importFrom SummarizedExperiment rowData
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply
#' @importFrom rtracklayer readGFFAsGRanges
#' @importFrom data.table as.data.table
#' @importFrom data.table setnames
#' @importFrom ggridges geom_density_ridges
#' @importFrom GenomicRanges gaps
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges ranges
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom ggbio autoplot
#' @importFrom cowplot plot_grid
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom ggrepel geom_text_repel
#' @import ggplot2
#'
#' @param object An scATSDataSet from TSS.
#' @param bam Path of bam file.
#' @param gtf Path of GTF file or GRanges file of GTF file.
#' @param gene Gene for plot.
#' @param ExonicReadsOnly Only the R2 of reads  within annotated exon region are remained or not.
#' @param CellBarcodeName The column name of cell barcode (CB) in Seurat object meta.data.
#' @param sep The character string to separate the CB (cell barcode) and sample name.
#' @param MaxSigma,MinProp,TSS Limitation of inferred TSS for plotting.
#' @param groupBy The column name of group in scATSDataSet object meta.data.
#' @param groups,groupLevels Only for given groups.
#' @param mapqFilter The minimum map quqlity.
#' @param isSupplementaryAlignment,isSecondaryAlignment parameters for ScanBamParam.
#' @param binwidth binwidth for histogram plot.
#' @param adjust adjust for density ridges plot.
#' @param base_size,FillSashimi,ColourTxName Plot parameters.
#' @param rel_heights Relative heights of final plots.
#'
#'
#' @concept visualization

TSSPlot <- function(object, bam, gtf, txdb = NULL, gene = NULL, Cells = NULL, ExonicReadsOnly = FALSE, CellBarcodeName = NULL, sep = "_",
                    MaxSigma = Inf, TSS = NULL, MinProp = 0, groupBy = NULL, groups = NULL, groupLevels = NULL, adjust = 1,
                    mapqFilter = 255, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = FALSE,
                    binwidth = 25, base_size = 15, FillSashimi = "#FEB24C", ColourTxName = "Red", rel_heights = c(1, 1, 5, 1)) {
  options(warn = -1)
  stopifnot(is(object, "scATSDataSet"))
  stopifnot(all(file.exists(bam)))
  stopifnot(all(file.exists(paste0(bam, ".bai"))))

  if(length(bam) > 1) {
    stopifnot(!is.null(sep))
    stopifnot(!is.na(sep))
    SampleName <- sort(unique(mapply(function(x) x[2], strsplit(colnames(object), split = sep))))
    stopifnot(identical(SampleName, sort(gsub(".bam", "", basename(bam)))))
  }

  if(is.null(txdb)) {
    txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtf)))
  } else {
    if(!is(txdb, "TxDb")) {
      txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtf)))
    }
  }
  EBG <- GenomicFeatures::exonsBy(txdb, by = "gene")

  if(!is(gtf, "GRanges")) {
    gtf <- rtracklayer::readGFFAsGRanges(gtf)
  } else {
    stopifnot(all(c("gene_name", "type") %in% names(S4Vectors::mcols(gtf))))
  }

  if(is.null(gene)) {
    stop("There is no gene for ATS identification!")
  } else {
    stopifnot(is.element(gene, SummarizedExperiment::rowData(object)$gene_name))
    stopifnot(is.element(gene, gtf$gene_name))
  }

  if(!is.null(groupBy)) {
    stopifnot(!is.null(CellBarcodeName))
  }

  stopifnot(mapqFilter >= 0)
  stopifnot(binwidth > 0)

  region <- gtf[with(data.frame(S4Vectors::mcols(gtf)), type == "gene" & gene_name == gene)]
  out <- rss(region = region, bam = bam, EBG = EBG, ExonicReadsOnly = ExonicReadsOnly, mapqFilter = mapqFilter, sep = sep)
  pos <- out[[2]]
  cov <- out[[1]]

  if(!is.null(groupBy)) {
    pos <- data.table::as.data.table(merge(pos, SummarizedExperiment::colData(object)[, c(groupBy, CellBarcodeName)], by.x = "CB", by.y = CellBarcodeName))
  }

  if(!is.null(Cells)) {
    pos <- pos[CB %in% Cells]
  }

  paras <- subset(SummarizedExperiment::rowData(object), gene_name == gene)
  paras <- subset(paras, sigma < MaxSigma & pro > MinProp)

  if(nrow(paras) > 0) {
    Xmin <- min(paras$mean); Xmax <- max(paras$mean)
  } else {
    Xmin <- pos[, min(pos)]; Xmax <- pos[, max(pos)]
  }

  if(!is.null(TSS)) paras <- paras[row.names(paras) %in% TSS, ]
  exb <- max(round(max(paras$sigma) * 4), 1000)

  dx <- seq.int(Xmin - exb, Xmax + exb, length.out = 10240)

  paras <- lapply(seq_len(nrow(paras)), function(i) {
    data.table::data.table(TSS = row.names(paras)[i], x = dx, y = dnorm(dx, mean = paras[i, "mean"], sd = paras[i, "sigma"]))
  })
  paras <- do.call(rbind, paras)

  ggplot(cov) +
    geom_ribbon(aes(x = pos, ymax = D, ymin = 0), fill = FillSashimi) +
    lims(x = c(Xmin - exb, Xmax + exb)) +
    labs(y = "Depth") +
    scale_y_continuous(n.breaks = 3) +
    theme_light(base_size = base_size) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> p1

  ggplot() +
    geom_histogram(data = pos, aes(pos, after_stat(density)), binwidth = binwidth, fill = "grey50") +
    geom_ribbon(data = paras, aes(x = x, ymax = y, ymin = 0, fill = TSS), alpha = 0.5) +
    ggrepel::geom_text_repel(data = paras[, .SD[which.max(y), ], TSS], aes(x, y, label = TSS), min.segment.length = 1, direction = "y", nudge_y = paras[y > 0, mean(y)]) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(y = "TSS") +
    lims(x = c(Xmin - exb, Xmax + exb)) +
    scale_y_continuous(breaks = paras[, max(y)]/2, labels = "Density") +
    theme_light(base_size = base_size) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.ticks = element_blank(),
          legend.background = element_blank(),
          # legend.position = c(Legend.x, Legend.y),
          legend.position = "none",
          legend.title = element_blank()) -> p2

  if(!is.null(groupBy)) {
    data.table::setnames(pos, groupBy, "groupBy")
    if(is.null(groups)) groups <- unique(pos[, groupBy])
    pos2 <- pos[groupBy %in% groups]
    if(!is.null(groupLevels)) pos2[, groupBy := factor(groupBy, levels = groupLevels)]
    ggplot(pos2, aes(x = pos, y = groupBy, fill = groupBy, height = ..density..)) +
      ggridges::geom_density_ridges(scale = 2, stat = "density", adjust = adjust, bw = binwidth) +
      lims(x = c(Xmin - exb, Xmax + exb)) +
      theme_light(base_size = base_size) +
      theme(legend.position = "none",
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()) -> p3
  } else {
    p3 <- NULL
  }

  ref1 <- gtf[with(data.frame(S4Vectors::mcols(gtf)), gene_name == gene & type == "exon")]
  ref1 <- split(ref1, ref1$transcript_name)
  ref1 <- ref1[mapply(ref1, FUN = function(x) length(IRanges::findOverlaps(IRanges::ranges(x), pos[, IRanges::IRanges(min(pos), max(pos))])) > 0)]
  ref1 <- lapply(ref1, function(x) {
    if(length(x) == 1) return(x)
    ins <- GenomicRanges::gaps(x)
    ins <- ins[GenomicRanges::start(ins) != 1]
    ins$type <- "gap" # "cds", "exon", "utr", "gap"
    sort(c(x, ins))
  })
  ref1 <- GenomicRanges::GRangesList(ref1)

  p4 <- lapply(seq_along(ref1), function(x) {
    p <- ggbio::autoplot(GenomicRanges::GRangesList(ref1[[x]]))@ggplot
    p + lims(x = c(Xmin - exb, Xmax + exb)) + theme_void(base_size = base_size) +
      geom_text(aes(x = Xmin - exb, y = 1, label = names(ref1)[x]), colour = ColourTxName)
  })
  p4 <- cowplot::plot_grid(plotlist = p4, ncol = 1, align = "v")

  if(is.null(groupBy)) {
    suppressWarnings(cowplot::plot_grid(p1, p2, p4, ncol = 1, align = "v", rel_heights = rel_heights[c(1, 2, 4)]))
  } else {
    suppressWarnings(cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, align = "v", rel_heights = rel_heights))
  }
}

#' @title FeaturePlot
#' @importFrom SummarizedExperiment assayNames colData
#' @importFrom reshape2 melt
#' @importFrom scales muted alpha
#'
#' @param object An scATSDataSet from TSS.
#' @param assay Which assay to pull expression data from?
#' @param cells,TSSs Which cells or TSSs for plot.
#' @param reduction Which reduction result are used for plot.
#' @param dim Which dimension of reduction result are used for plot.
#' @param groupBy The column name of group in scATSDataSet object meta.data.
#' @param groups Only for given groups.
#' @param plot Plot or return data frame.
#' @param shapeForGroup Given shape for group or not?
#' @param order Order the point by value or not?
#' @param point.size,plot.theme,color.low,color.high,color.na,midpoint,color.mid Parameters for ggplot.
#' @param ... Additional parameters for \code{scale_colour_gradient2}
#'
#' @export
#' @concept visualization


FeaturePlot <- function(object,
                        assay = "psi",
                        cells = NULL,
                        TSSs = NULL,
                        reduction = "umap",
                        dim = 1:2,
                        logScale = FALSE,
                        groupBy = NULL,
                        groups = NULL,
                        plot = TRUE,
                        shapeForGroup = FALSE,
                        point.size = 1,
                        order = TRUE,
                        plot.theme = NULL,
                        base_size = 15,
                        color.high = "#642594",
                        color.mid = NULL,
                        color.low = NULL,
                        midpoint = NULL,
                        color.na = "grey90", ...) {
  stopifnot(is(object, "scATSDataSet"))

  if(length(slot(object, "reductions")) == 0) {
    stop("There is no reduction results.")
  }

  stopifnot(is.logical(plot))
  stopifnot(is.logical(shapeForGroup))
  stopifnot(is.logical(order))
  stopifnot(is.logical(logScale))
  if(shapeForGroup) stopifnot(!is.null(groupBy))

  stopifnot(!is.null(color.high))
  color.high <- scales::muted(color.high)
  if(!is.null(color.mid)) color.mid <- scales::muted(color.mid)
  color.mid <- color.mid %||% scales::alpha(color.high, 0.5)
  if(!is.null(color.low)) color.low <- scales::muted(color.low)
  color.low <- color.low %||% scales::alpha(color.high, 0.05)

  if(!is.null(assay)) stopifnot(length(assay) == 1)
  if(!is.null(assay)) stopifnot(is.element(assay, SummarizedExperiment::assayNames(object)))
  assay <- assay %||% "psi"

  if(!is.null(cells)) stopifnot(all(is.element(cells, colnames(object))))
  cells <- cells %||% colnames(object)

  if(!is.null(TSSs)) stopifnot(all(is.element(TSSs, row.names(object))))
  if(is.null(TSSs)) stop("No TSSs"); if(length(TSSs) == 0) stop("No TSSs")

  stopifnot(length(dim) == 2)
  if(!is.null(groupBy)) stopifnot(length(groupBy) == 1)
  if(!is.null(groupBy)) stopifnot(groupBy %in% colnames(SummarizedExperiment::colData(object)))

  if(assay == "psi") {
    data <- psi(object, cells, TSSs)
  }
  if(assay == "counts") {
    data <- counts(object, cells, TSSs)
  }
  if(assay == "expected.count") {
    data <- expected.count(object, cells, TSSs)
  }

  if(!is.null(reduction)) stopifnot(length(reduction) == 1)
  if(!is.null(reduction)) {
    if(!is.element(reduction, names(slot(object, "reductions")))) {
      stop(paste("reduction should be one of", paste(names(slot(object, "reductions")), collapse = ", ")))
    }
  }
  reduction <- reduction %||% setdiff(names(slot(object, "reductions")), "pca")[1]
  reduct <- slot(object, "reductions")[[which(names(slot(object, "reductions")) == reduction)]]
  rds <- paste0(slot(reduct, "key"), dim)
  if(!all(is.element(rds, colnames(reduct@cell.embeddings)))) {
    stop("Given dim are not in the reduction results.")
  }
  reduct <- reduct@cell.embeddings[, rds]
  data <- merge(reshape2::melt(t(data)), reduct, by.x = "Var1", by.y = 0)

  if(!is.null(groupBy)) {
    data <- merge(data, SummarizedExperiment::colData(object)[groupBy], by.x = "Var1", by.y = 0)
  }
  data <- data.table::as.data.table(data)
  if(!is.null(groups)) {
    data[!data[, groupBy, with = F][[1]] %in% groups, value := NaN]
  }

  if(order) {
    data <- rbind(data[is.na(value)], data[!is.na(value)][order(value, decreasing = F)])
  }
  if(!plot) return(data)

  if(is.null(plot.theme)) {
    plot.theme <- theme_bw(base_size = base_size) +
      theme(panel.grid = element_blank(), strip.background = element_blank())
  } else {
    stopifnot(is(plot.theme, "theme"))
  }

  if(logScale) {
    midpoint <- midpoint %||% mean(c(data[, min(log1p(value), na.rm = TRUE)], data[, max(log1p(value), na.rm = TRUE)]))
    if(shapeForGroup) {
      ggplot(data, aes(x = eval(expr = parse(text = rds[1])), y = eval(expr = parse(text = rds[2])), colour = log1p(value))) +
        geom_point(mapping = aes(shape = eval(expr = parse(text = groupBy))), point.size = point.size) +
        labs(x = rds[1], y = rds[2]) +
        facet_wrap(~ Var2) +
        scale_colour_gradient2(low = color.low,
                               mid = color.mid,
                               high = color.high,
                               na.value = color.na,
                               midpoint = midpoint,
                               guide = guide_colorbar(title = paste0("log(", assay, ")"),
                                                      title.position = "left",
                                                      title.theme = element_text(angle = 90, hjust = 0.5)), ...) +
        scale_shape_manual(values = seq_along(unique(data[, groupBy, with = FALSE][[1]])),
                           guide = guide_legend(title = groupBy)) +
        plot.theme
    } else {
      ggplot(data, aes(x = eval(expr = parse(text = rds[1])), y = eval(expr = parse(text = rds[2])), colour = log1p(value))) +
        geom_point(point.size = point.size) +
        labs(x = rds[1], y = rds[2]) +
        facet_wrap(~ Var2) +
        scale_colour_gradient2(low = color.low,
                               mid = color.mid,
                               high = color.high,
                               na.value = color.na,
                               midpoint = midpoint,
                               guide = guide_colorbar(title = paste0("log(", assay, ")"),
                                                      title.position = "left",
                                                      title.theme = element_text(angle = 90, hjust = 0.5)), ...) +
        plot.theme
    }
  } else {
    midpoint <- midpoint %||% mean(c(data[, min(value, na.rm = TRUE)], data[, max(value, na.rm = TRUE)]))
    if(shapeForGroup) {
      ggplot(data, aes(x = eval(expr = parse(text = rds[1])), y = eval(expr = parse(text = rds[2])), colour = value)) +
        geom_point(mapping = aes(shape = eval(expr = parse(text = groupBy))), point.size = point.size) +
        labs(x = rds[1], y = rds[2]) +
        facet_wrap(~ Var2) +
        scale_colour_gradient2(low = color.low,
                               mid = color.mid,
                               high = color.high,
                               na.value = color.na,
                               midpoint = midpoint,
                               guide = guide_colorbar(title = assay), ...) +
        scale_shape_manual(values = seq_along(unique(data[, groupBy, with = FALSE][[1]])),
                           guide = guide_legend(title = groupBy)) +
        plot.theme
    } else {
      ggplot(data, aes(x = eval(expr = parse(text = rds[1])), y = eval(expr = parse(text = rds[2])), colour = value)) +
        geom_point(point.size = point.size) +
        labs(x = rds[1], y = rds[2]) +
        facet_wrap(~ Var2) +
        scale_colour_gradient2(low = color.low,
                               mid = color.mid,
                               high = color.high,
                               na.value = color.na,
                               midpoint = midpoint,
                               guide = guide_colorbar(title = assay), ...) +
        plot.theme
    }
  }
}


PlotGene <- function(gtffile = NULL, txdb = NULL, gtfrange = NULL, gene) {
  stopifnot(!is.null(gene))
  stopifnot(length(gene) == 1)
  if(!is.null(txdb)) stopifnot(is(txdb, "TxDb"))
  if(!is.null(gtfrange)) stopifnot(is(gtfrange, "GRanges"))
  if(is.null(txdb) | is.null(gtfrange)) stopifnot(!is.null(gtffile))
  if(is.null(txdb) | is.null(gtfrange)) stopifnot(file.exists(gtffile))
  if(is.null(txdb)) {
    txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtffile)))
  }
  if(is.null(gtfrange)) gtfrange <- rtracklayer::import(gtffile)
  stopifnot(is.element(gene, with(as.data.frame(gtfrange), gene_id)) | is.element(gene, with(as.data.frame(gtfrange), gene_name)))
  if(is.element(gene, with(as.data.frame(gtfrange), gene_name))) {
    gene <- with(as.data.frame(gtfrange), gene_id[type == "gene" & gene_name == gene])
  }

  gi <- biovizBase::crunch(txdb, GenomicFeatures::genes(txdb)[gene])
  gi <- GenomicRanges::GRangesList(split(gi, gi$tx_name))
  names(gi) <- plyr::mapvalues(names(gi), with(as.data.frame(gtf), transcript_id[type == "transcript"]), with(as.data.frame(gtf), transcript_name[type == "transcript"]), warn_missing = F)
  ggbio::autoplot(gi, space.skip = 0)@ggplot
}


PlotR1 <- function(bam, gene, gtffile = NULL, txdb = NULL, gtfrange = NULL, exonIntersect = TRUE, junctionIntersect = FALSE, binwidth = 10, index = NULL,
                   mapqFilter = 255, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = FALSE, sep = "_", ggtheme = NULL) {
  stopifnot(!is.null(gene))
  stopifnot(length(gene) == 1)
  if(!is.null(txdb)) stopifnot(is(txdb, "TxDb"))
  if(!is.null(gtfrange)) stopifnot(is(gtfrange, "GRanges"))
  if(is.null(txdb) | is.null(gtfrange)) stopifnot(!is.null(gtffile))
  if(is.null(txdb) | is.null(gtfrange)) stopifnot(file.exists(gtffile))
  if(is.null(txdb)) {
    txdb <- suppressMessages(suppressWarnings(txdbmaker::makeTxDbFromGFF(gtffile)))
  }
  if(is.null(gtfrange)) gtfrange <- rtracklayer::import(gtffile)
  stopifnot(is.element(gene, with(as.data.frame(gtfrange), gene_id)) | is.element(gene, with(as.data.frame(gtfrange), gene_name)))
  if(is.element(gene, with(as.data.frame(gtfrange), gene_name))) {
    gene <- with(as.data.frame(gtfrange), gene_id[type == "gene" & gene_name == gene])
  }
  if(is.null(ggtheme)) {
    ggtheme <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
  }
  stopifnot(is(ggtheme, "theme"))


  map0 <- loadbamgene(bam = bam, gene = gene, txdb = txdb, exonIntersect = exonIntersect, junctionIntersect = junctionIntersect,
                      isSupplementaryAlignment = isSupplementaryAlignment, isSecondaryAlignment = isSecondaryAlignment,
                      isDuplicate = isDuplicate, mapqFilter = mapqFilter, sep = sep)

  if(S4Vectors::runValue(GenomicAlignments::strand(map0)) == "+") {
    ggplot(data.table(x = min(GenomicAlignments::start(GenomicAlignments::first(map0))):max(GenomicAlignments::start(GenomicAlignments::first(map0))),
                      y = ecdf(GenomicAlignments::start(GenomicAlignments::first(map0)))(min(GenomicAlignments::start(GenomicAlignments::first(map0))):max(GenomicAlignments::start(GenomicAlignments::first(map0)))) * length(map0)),
           aes(x, y)) +
      geom_step() +
      labs(y = "Reads start CDF", x = "Genomic loci") +
      ggtheme -> p0

    ggplot(data.table(x = GenomicAlignments::start(GenomicAlignments::first(map0))), aes(x)) +
      geom_histogram(binwidth = binwidth) +
      labs(y = "Reads start count", x = "Genomic loci") +
      ggtheme -> p1
  } else {
    ggplot(data.table(x = min(GenomicAlignments::end(GenomicAlignments::first(map0))):max(GenomicAlignments::end(GenomicAlignments::first(map0))),
                      y = (1 - ecdf(GenomicAlignments::end(GenomicAlignments::first(map0)))(min(GenomicAlignments::end(GenomicAlignments::first(map0))):max(GenomicAlignments::end(GenomicAlignments::first(map0))))) * length(map0)),
           aes(x, y)) +
      geom_step() +
      labs(y = "Reads start CDF", x = "Genomic loci") +
      ggtheme -> p0

    ggplot(data.table(x = GenomicAlignments::end(GenomicAlignments::first(map0))), aes(x)) +
      geom_histogram(binwidth = binwidth) +
      labs(y = "Reads start count", x = "Genomic loci") +
      ggtheme -> p1
  }

  test <- GenomicAlignments::coverage(GenomicAlignments::first(map0))[[unique(as.character(GenomeInfoDb::seqnames(map0)))]]
  test2 <- test[IRanges::start(range(IRanges::IRanges(test > 0))):IRanges::end(range(IRanges::IRanges(test > 0)))]

  ggplot(data.table(x = IRanges::start(range(IRanges::IRanges(test > 0))):IRanges::end(range(IRanges::IRanges(test > 0))),
                    y = as.numeric(test2)), aes(x = x, y = y)) +
    geom_step() +
    labs(y = "Reads coverage", x = "Genomic loci") +
    ggtheme -> p2

  p3 <- PlotGene(gtffile = gtffile, txdb = txdb, gtfrange = gtfrange, gene = gene) + ggtheme +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

  p0 <- p0 + scale_x_continuous(limits = layer_scales(p3)$x$range$range)
  p1 <- p1 + scale_x_continuous(limits = layer_scales(p3)$x$range$range)
  p2 <- p2 + scale_x_continuous(limits = layer_scales(p3)$x$range$range)

  if(!is.null(index)) {
    index <- index[index >= min(layer_scales(p3)$x$range$range) & index <= max(layer_scales(p3)$x$range$range)]
    if(length(index) > 0) {
      p0 <- p0 + geom_vline(xintercept = index, lty = 2, col = 2)
      p3 <- p3 + geom_vline(xintercept = index, lty = 2, col = 2)
    }
  }
  cowplot::plot_grid(p0, p1, p2, p3, ncol = 1, align = "v")
}

