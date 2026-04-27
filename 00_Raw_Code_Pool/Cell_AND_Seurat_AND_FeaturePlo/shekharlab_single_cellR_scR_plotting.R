#require(ggjoy)

#' Visualize features on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#'
#' @param object scR object
#' @param features.plot Vector of features to plot (typically genes)
#' @param min.cutoff Vector of minimum cutoff values for each feature, may specify quantile in the form of 'q##' where '##' is the quantile (eg, 1, 10)
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param cols.use The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param pch.use Pch for plotting
#' @param overlay Plot two features overlayed one on top of the other
#' @param do.hover Enable hovering over points to view information
#' @param data.hover Data to add to the hover, pass a character vector of features to add. Defaults to cell name and identity. Pass 'NULL' to remove extra data.
#' @param do.identify Opens a locator session to identify clusters of cells
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "tsne", can also be "pca", or "ica", assuming these are precomputed.
#' @param use.imputed Use imputed values for gene expression (default is FALSE)
#' @param nCol Number of columns to use when plotting multiple features.
#' @param no.axes Remove axis labels
#' @param no.legend Remove legend from the graph. Default is TRUE.
#' @param dark.theme Plot in a dark theme
#' @param do.return return the ggplot2 object
#'
#' @importFrom RColorBrewer brewer.pal.info
#'
#' @return No return value, only a graphical output
#'
#' @export
#'
#' @examples
#' FeaturePlot(object = pbmc_small, features.plot = 'PC1')
#'
FeaturePlot <- function(
  object,
  features.plot,
  min.cutoff = NA,
  max.cutoff = NA,
  dim.1 = 1,
  dim.2 = 2,
  cells.use = NULL,
  pt.size = 1,
  cols.use = c("yellow", "red"),
  pch.use = 16,
  overlay = FALSE,
  do.hover = FALSE,
  data.hover = 'ident',
  do.identify = FALSE,
  reduction.use = "tsne",
  use.imputed = FALSE,
  nCol = NULL,
  no.axes = FALSE,
  no.legend = TRUE,
  dark.theme = FALSE,
  do.return = FALSE,
  use.raw = FALSE,
  use.count = FALSE,
  max.cells = NA
) {
  cells.use <- set.ifnull(cells.use, colnames(object@data))
  if (!is.na(max.cells)){
    cells.use = sample(cells.use, max.cells)
    }
  if (is.null(nCol)) {
    nCol <- 2
    if (length(features.plot) == 1) {
      nCol <- 1
    }
    if (length(features.plot) > 6) {
      nCol <- 3
    }
    if (length(features.plot) > 9) {
      nCol <- 4
    }
  }
  
  num.row <- floor(x = length(features.plot) / nCol - 1e-5) + 1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimCode(reduction.use)
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  print(1)
  if (reduction.use == "tsne") dim.codes = c("tsne_1","tsne_2")
  data.plot <- as.data.frame(GetCellEmbeddings(
    object = object,
    reduction.use = reduction.use,
    dims.use = c(dim.1, dim.2),
    cells.use = cells.use
  ))

  if (reduction.use == "tsne") colnames(data.plot) = c("tsne_1","tsne_2")
  data.plot$x <- data.plot[, dim.codes[1]]
  data.plot$y <- data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  #names(x = data.plot) <- c('x', 'y')

  data.use=t(fetch.data(
    object,
    features.plot,
    cells.use = cells.use, 
    use.raw=use.raw, 
    use.count=use.count))
  
  #   Check mins and maxes
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      ifelse(
        test = is.na(x = cutoff),
        yes = min(data.use[feature, ]),
        no = cutoff
      )
    },
    cutoff = min.cutoff,
    feature = features.plot
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      ifelse(
        test = is.na(x = cutoff),
        yes = max(data.use[feature, ]),
        no = cutoff
      )
    },
    cutoff = max.cutoff,
    feature = features.plot
  )
  check_lengths = unique(x = vapply(
    X = list(features.plot, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check_lengths) != 1) {
    stop('There must be the same number of minimum and maximum cuttoffs as there are features')
  }
  if (overlay) {
    stop("Overlay is not yet implemented")
    #   Wrap as a list for MutiPlotList
    # pList <- list(
    #   BlendPlot(
    #     data.use = data.use,
    #     features.plot = features.plot,
    #     data.plot = data.plot,
    #     pt.size = pt.size,
    #     pch.use = pch.use,
    #     cols.use = cols.use,
    #     dim.codes = dim.codes,
    #     min.cutoff = min.cutoff,
    #     max.cutoff = max.cutoff,
    #     no.axes = no.axes,
    #     no.legend = no.legend,
    #     dark.theme = dark.theme
    #   )
    # )
  } else {
    #   Use mapply instead of lapply for multiple iterative variables.
    pList <- mapply(
      FUN = SingleFeaturePlot,
      feature = features.plot,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff,
      MoreArgs = list( # Arguments that are not being repeated
        data.use = data.use,
        data.plot = data.plot,
        pt.size = pt.size,
        pch.use = pch.use,
        cols.use = cols.use,
        dim.codes = dim.codes,
        no.axes = no.axes,
        no.legend = no.legend,
        dark.theme = dark.theme
      ),
      SIMPLIFY = FALSE # Get list, not matrix
    )
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot") # Instead just skipt this step
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    } else {
      features.info <- fetch.data(object = object, vars.all = data.hover)
    }
    #   Use pList[[1]] to properly extract the ggplot out of the plot list
    return(HoverLocator(
      plot = pList[[1]],
      data.plot = data.plot,
      features.info = features.info,
      dark.theme = dark.theme,
      title = features.plot
    ))
    # invisible(readline(prompt = 'Press <Enter> to continue\n'))
  } else if (do.identify) {
    require(SDMTools)
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    #   Use pList[[1]] to properly extract the ggplot out of the plot list
    return(FeatureLocator(
      plot = pList[[1]],
      data.plot = data.plot,
      dark.theme = dark.theme
    ))
  } else {
    print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
  }
  ResetPar()
  if (do.return){
    return(pList)
  }
}


#' Plot a single feature
#'
#' @param data.use The data regarding the feature
#' @param feature The feature to plot
#' @param data.plot The data to be plotted
#' @param pt.size Size of each point
#' @param pch.use Shape of each point
#' @param cols.use Colors to plot
#' @param dim.codes Codes for the dimensions to plot in
#' @param min.cutoff Minimum cutoff for data
#' @param max.cutoff Maximum cutoff for data
#' @param no.axes Remove axes from plot
#' @param no.legend Remove legend from plot
#' @param dark.theme Plot in dark theme
#
# @return A ggplot2 scatterplot
#
#
SingleFeaturePlot <- function(
  data.use,
  feature,
  data.plot,
  pt.size,
  pch.use,
  cols.use,
  dim.codes,
  min.cutoff,
  max.cutoff,
  no.axes,
  no.legend,
  dark.theme
) {
  data.gene <- na.omit(object = data.frame(data.use[feature, ]))
  #   Check for quantiles
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.gene)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.gene)
  #   Mask any values below the minimum and above the maximum values
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x < min.cutoff, yes = min.cutoff, no = x))
    }
  )
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x > max.cutoff, yes = max.cutoff, no = x))
    }
  )
  data.plot$gene <- data.gene
  #   Stuff for break points
  if (length(x = cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  } else {
    brewer.gran <- length(x = cols.use)
  }
  #   Cut points
  if (all(data.gene == 0)) {
    data.cut <- 0
  } else {
    data.cut <- as.numeric(x = as.factor(x = cut(
      x = as.numeric(x = data.gene),
      breaks = brewer.gran
    )))
  }
  data.plot$col <- as.factor(x = data.cut)
  #   Start plotting
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  if (brewer.gran != 2) {
    if (length(x = cols.use) == 1) {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_brewer(palette = cols.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_manual(values = cols.use)
    }
  } else {
    if (all(data.plot$gene == data.plot$gene[1])) {
      warning(paste0("All cells have the same value of ", feature, "."))
      p <- p + geom_point(color = cols.use[1], size = pt.size, shape = pch.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = gene),
        size = pt.size,
        shape = pch.use
      ) + scale_color_gradientn(
        colors = cols.use,
        guide = guide_colorbar(title = feature)
      )
    }
  }
  if (no.axes) {
    p <- p + labs(title = feature, x ="", y="") + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  } else {
    p <- p + labs(title = feature, x = dim.codes[1], y = dim.codes[2])
  }
  if (no.legend) {
    p <- p + theme(legend.position = 'none')
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  return(p)
}

#' FeaturePlot for two features
BiFeaturePlot <- function(
  object,
  features.plot,
  quantile.cutoff=NA,
  dim.1 = 1,
  dim.2 = 2,
  cells.use = NULL,
  breaks = NULL,
  pt.size = 1,
  cols.use = c("red", "blue"),
  pch.use = 16,
  overlay = FALSE,
  do.hover = FALSE,
  data.hover = 'ident',
  do.identify = FALSE,
  reduction.use = "tsne",
  use.imputed = FALSE,
  nCol = NULL,
  no.axes = FALSE,
  no.legend = TRUE,
  dark.theme = FALSE,
  do.return = FALSE,
  use.raw = FALSE,
  use.count = FALSE
) {
  require(colorspace)
  if (!is.list(features.plot)){
    stop("Error: features.plot must be a list with each element being a character vector of length 2")
  }
  all.features = unique(unlist(features.plot)) 
  cells.use <- set.ifnull(cells.use, colnames(object@data))
  if (is.null(nCol)) {
    nCol <- 2
    if (length(features.plot) == 1) {
      nCol <- 1
    }
    if (length(features.plot) > 6) {
      nCol <- 3
    }
    if (length(features.plot) > 9) {
      nCol <- 4
    }
  }
  
  num.row <- floor(x = length(features.plot) / nCol - 1e-5) + 1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimCode(reduction.use)
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  print(1)
  if (reduction.use == "tsne") dim.codes = c("tsne_1","tsne_2")
  data.plot <- as.data.frame(GetCellEmbeddings(
    object = object,
    reduction.use = reduction.use,
    dims.use = c(dim.1, dim.2),
    cells.use = cells.use
  ))
  
  if (reduction.use == "tsne") colnames(data.plot) = c("tsne_1","tsne_2")
  data.plot$x <- data.plot[, dim.codes[1]]
  data.plot$y <- data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  #names(x = data.plot) <- c('x', 'y')
  
  data.use=t(fetch.data(
    object,
    all.features,
    cells.use = cells.use, 
    use.raw=use.raw, 
    use.count=use.count))
  
  # Normalize data
  print("Normalizing each feature from min to max")
  for (f in all.features){
    n = length(cells.use)
    if (is.na(quantile.cutoff) | quantile.cutoff > 1){
      cutoff = quantile(data.use[f,], max(0.99, 1-100/n))
    } else {
      cutoff = quantile(data.use[f,], quantile.cutoff)
    }
    data.use[f,] = data.use[f,]/cutoff
  }
 
  if (overlay) {
    stop("Overlay is not yet implemented")
    #   Wrap as a list for MutiPlotList
    # pList <- list(
    #   BlendPlot(
    #     data.use = data.use,
    #     features.plot = features.plot,
    #     data.plot = data.plot,
    #     pt.size = pt.size,
    #     pch.use = pch.use,
    #     cols.use = cols.use,
    #     dim.codes = dim.codes,
    #     min.cutoff = min.cutoff,
    #     max.cutoff = max.cutoff,
    #     no.axes = no.axes,
    #     no.legend = no.legend,
    #     dark.theme = dark.theme
    #   )
    # )
  } else {
    #   Use mapply instead of lapply for multiple iterative variables.
    pList <- mapply(
      FUN = DoubleFeaturePlot,
      features = features.plot,
      MoreArgs = list( # Arguments that are not being repeated
        data.use = data.use,
        data.plot = data.plot,
        breaks = breaks,
        pt.size = pt.size,
        pch.use = pch.use,
        cols.use = cols.use,
        dim.codes = dim.codes,
        no.axes = no.axes,
        no.legend = no.legend,
        dark.theme = dark.theme
      ),
      SIMPLIFY = FALSE # Get list, not matrix
    )
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot") # Instead just skipt this step
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    } else {
      features.info <- fetch.data(object = object, vars.all = data.hover)
    }
    #   Use pList[[1]] to properly extract the ggplot out of the plot list
    return(HoverLocator(
      plot = pList[[1]],
      data.plot = data.plot,
      features.info = features.info,
      dark.theme = dark.theme,
      title = features.plot
    ))
    # invisible(readline(prompt = 'Press <Enter> to continue\n'))
  } else if (do.identify) {
    require(SDMTools)
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    #   Use pList[[1]] to properly extract the ggplot out of the plot list
    return(FeatureLocator(
      plot = pList[[1]],
      data.plot = data.plot,
      dark.theme = dark.theme
    ))
  } else {
    print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
  }
  ResetPar()
  if (do.return){
    return(pList)
  }
}

#' Single Feature plot for two features

DoubleFeaturePlot = function(
  data.use,
  features,
  data.plot,
  breaks,
  pt.size,
  pch.use,
  cols.use,
  dim.codes,
  no.axes,
  no.legend,
  dark.theme
) {
  if (length(features) != 2){
    stop("Error: DoubleFeaturePlot works only for two features. Please supply a list where each element is a vector of two features")
  }
  data.gene <- na.omit(object = data.frame(data.use[features, ]))
  #   Check for quantiles
  #   Mask any values below the minimum and above the maximum values

  #   Stuff for break points
  if (is.null(breaks)){
    brewer.gran = 2
  } else {
    brewer.gran = breaks
  }
  if (!(brewer.gran == round(brewer.gran)) | (brewer.gran < 0)) {
    stop("breaks must be a positive integer")
  }
  print(paste0("using ", brewer.gran, " breaks" ))
  
  #   Cut points
  l=1
  for (f in features){
    if (all(data.gene == 0)) {
      data.cut <- 0
    } else {
      data.cut <- as.numeric(x = as.factor(x = cut(
        x = as.numeric(x = data.gene[f,]),
        breaks = brewer.gran
      )))
    }
    colx = cols.use[l]
    colramp = colorRampPalette(c("white",colx))(breaks)
    
    data.plot[,paste0("col",l)] <- sapply(data.cut, function(x) colramp[x] )
    l=l+1
  }
  
  data.plot[,features[1]] = as.numeric(data.gene[features[1],])
  data.plot[,features[2]] = as.numeric(data.gene[features[2],])
  
  data.plot$col = sapply(1:nrow(data.plot), function(x){
    x1 = as.numeric(col2rgb(data.plot[x,"col1"]))/255
    x2 = as.numeric(col2rgb(data.plot[x,"col2"]))/255
    m = as.numeric(colorspace::coords(mixcolor(0.5, RGB(x1[1],x1[2],x1[3]), RGB(x2[1],x2[2],x2[3]))))
    return(rgb(m[1],m[2],m[3]))
  })
  
  t1 = paste0(features[1], "(", cols.use[1], ")")
  t2 = paste0(features[2], "(", cols.use[2], ")")
  title.use = paste0(t1,", ", t2)

    #   Start plotting
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  
  p <- p + geom_point(
    size = pt.size,
    shape = pch.use,
    color = data.plot$col
  ) + scale_color_identity()
  if (no.axes) {
    p <- p + labs(title = title.use, x ="", y="") + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  } else {
    p <- p + labs(title = title.use, x = dim.codes[1], y = dim.codes[2])
  }
  if (no.legend) {
    p <- p + theme(legend.position = 'none')
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  return(p)
}

#' Plot tSNE map
#'
#' Graphs the output of a tSNE analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object scR object
#' @param do.label FALSE by default. If TRUE, plots an alternate view where the center of each
#' cluster is labeled
#' @param pt.size Set the point size
#' @param label.size Set the size of the text labels
#' @param cells.use Vector of cell names to use in the plot.
#' @param colors.use Manually set the color palette to use for the points
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#'
#' @seealso DimPlot
#'
#' @export
#'
#' @examples
#' TSNEPlot(object = pbmc_small)
#'
TSNEPlot <- function(
  object,
  do.label = TRUE,
  pt.size=1,
  label.size=6,
  cells.use = NULL,
  cols.use = NULL,
  group.by = "ident",
  do.return = FALSE,
  ...
) {
  return(DimPlot(
    object = object,
    reduction.use = "tsne",
    cells.use = cells.use,
    pt.size = pt.size,
    do.label = do.label,
    label.size = label.size,
    cols.use = cols.use,
    do.return = do.return,
    ...
  ))
}

#' Plot PCA map
#'
#' Graphs the output of a PCA analysis
#' Cells are colored by their identity class.
#'
#' This function is a wrapper for DimPlot. See ?DimPlot for a full list of possible
#' arguments which can be passed in here.
#'
#' @param object scR object
#' @param \dots Additional parameters to DimPlot, for example, which dimensions to plot.
#'
#' @export
#'
#' @examples
#' PCAPlot(object = pbmc_small)
#'
PCAPlot <- function(
  object, 
  dim.1=1, dim.2=2, 
  pt.size=1,
  group.by = "ident",
  ...) {
  return(DimPlot(object = object, reduction.use = "pca", 
                 dim.1=dim.1, dim.2=dim.2, 
                 pt.size=pt.size, 
                 group.by = group.by,
                 label.size = 4, ...))
}

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique (PCA by default).
#' Cells are colored by their identity class.
#'
#' @param object scR object
#' @param reduction.use Which dimensionality reduction to use. Default is
#' "pca", can also be "tsne", or "ica", assuming these are precomputed.
#' @param dim.1 Dimension for x-axis (default 1)
#' @param dim.2 Dimension for y-axis (default 2)
#' @param cells.use Vector of cells to plot (default is all cells)
#' @param pt.size Adjust point size for plotting
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param do.bare Do only minimal formatting (default : FALSE)
#' @param cols.use Vector of colors, each color corresponds to an identity
#' class. By default, ggplot assigns colors.
#' @param group.by Group (color) cells in different ways (for example, orig)
#' @param pt.shape If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with fetch.data) allowing for both different colors and
#' different shapes on cells.
#' @param do.hover Enable hovering over points to view information
#' @param data.hover Data to add to the hover, pass a character vector of features to add. Defaults to cell name and ident. Pass 'NULL' to clear extra information.
#' @param do.identify Opens a locator session to identify clusters of cells.
#' @param do.label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param no.legend Setting to TRUE will remove the legend
#' @param no.axes Setting to TRUE will remove the axes
#' @param dark.theme Use a dark theme for the plot
#' @param ... Extra parameters to FeatureLocator for do.identify = TRUE
#'
#' @return If do.return==TRUE, returns a ggplot2 object. Otherwise, only
#' graphical output.
#'
#' @seealso \code{FeatureLocator}
#'
#' @import SDMTools
#' @importFrom stats median
#' @importFrom dplyr summarize group_by
#'
#' @export
#'
#' @examples
#' DimPlot(object = pbmc_small)
#'
DimPlot <- function(
  object,
  reduction.use = "pca",
  dim.1 = 1,
  dim.2 = 2,
  cells.use = NULL,
  pt.size = 1,
  do.return = FALSE,
  do.bare = FALSE,
  cols.use = NULL,
  group.by = "ident",
  pt.shape = NULL,
  do.hover = FALSE,
  data.hover = 'ident',
  do.identify = FALSE,
  do.label = FALSE,
  label.size = 4,
  no.legend = FALSE,
  no.axes = FALSE,
  dark.theme = FALSE,
  max.cells = NA,
  ...
) {
  
  cells.use <- set.ifnull(cells.use, colnames(object@data))
  embeddings.use = GetCellEmbeddings(object = object, reduction.use = reduction.use, dims.use = c(dim.1, dim.2),
                                     cells.use=cells.use)
  if (reduction.use == "tsne") colnames(embeddings.use) = c("tsne_1","tsne_2")
  if (length( embeddings.use) == 0) {
    stop(paste(reduction.use, "has not been run for this object yet."))
  }
  dim.code <- GetDimCode(reduction.use = reduction.use)
  if (dim.code == "tsne") dim.code = paste0(dim.code,"_")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(embeddings.use)
  ident.use <- as.factor(object@ident[cells.use])
  if (group.by != "ident") {
    ident.use <- as.factor(fetch.data(
      object = object,
      vars.all = group.by
    )[cells.use, 1])
    names(ident.use) = cells.use
  }

  data.plot$ident <- ident.use
  data.plot$x <- data.plot[, dim.codes[1]]
  data.plot$y <- data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  
  # Downsample if necessary
  if (!is.na(max.cells)){
    data.plot = data.plot[sample(rownames(data.plot),max.cells),]   
  } else {
      data.plot=data.plot[sample(rownames(data.plot)),]
  }
  
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
    geom_point(mapping = aes(colour = factor(x = ident)), size = pt.size)
  if (! is.null(pt.shape)) {
    shape.val <- fetch.data(object = object, vars.all = pt.shape)[cells.use, 1]; names(shape.val) = cells.use
    if (is.numeric(shape.val)) {
      shape.val <- cut(x = shape.val, breaks = 5)
    }
    data.plot[, "pt.shape"] <- shape.val
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
      geom_point(
        mapping = aes(colour = factor(x = ident), shape = factor(x = pt.shape)),
        size = pt.size
      )
  }
  if (! is.null(cols.use)) {
    p <- p + scale_colour_manual(values = cols.use)
  } else{
    cols.use=set.ifnull(rainbow(length(levels(data.plot$ident)))); cols.use[1]="lightgrey"
    p <- p + scale_colour_manual(values = cols.use)
  }
  p2 <- p +
    xlab(label = dim.codes[[1]]) +
    ylab(label = dim.codes[[2]]) +
    scale_size(range = c(pt.size, pt.size))
  p3 <- p2 +
    SetXAxisGG() +
    SetYAxisGG() +
    SetLegendPointsGG(x = 6) +
    SetLegendTextGG(x = 12) +
    no.legend.title +
    theme_bw() +
    nogrid
  p3 <- p3 + theme(legend.title = element_blank())
  if (do.label) {
    flag=0
    if("plyr" %in% (.packages())){
      detach(package:Rmisc)
      detach(package:plyr)
      flag=1
    }
    data.plot %>%
      dplyr::group_by(ident) %>%
      summarise(x = median(x = x), y = median(x = y)) -> centers
    centers <- as.data.frame(centers)
    #return(centers)
    p3 <- p3 +
      geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0) +
      geom_text(data = centers, mapping = aes(label = ident), size = label.size)

    if(flag==1){
      library(Rmisc)
      library(plyr)
    }
    
  }
  if (dark.theme) {
    p <- p + DarkTheme()
    p3 <- p3 + DarkTheme()
  }
  if (no.legend) {
    p3 <- p3 + theme(legend.position = "none")
  }
  if (no.axes) {
    p3 <- p3 + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
  }
  if (do.identify || do.hover) {
    if (do.bare) {
      plot.use <- p
    } else {
      plot.use <- p3
    }
    if (do.hover) {
      if (is.null(x = data.hover)) {
        features.info <- NULL
      } else {
        features.info <- fetch.data(object = object, vars.all = data.hover)
      }
      return(HoverLocator(
        plot = plot.use,
        data.plot = data.plot,
        features.info = features.info,
        dark.theme = dark.theme
      ))
    } else if (do.identify) {
      return(FeatureLocator(
        plot = plot.use,
        data.plot = data.plot,
        dark.theme = dark.theme,
        ...
      ))
    }
  }
  if (do.return) {
    if (do.bare) {
      return(p)
    } else {
      return(p3)
    }
  }
  if (do.bare) {
    print(p)
  } else {
    print(p3)
  }
}

######################## JOY PLOT ###################################

#' Single cell joy plot
#'
#' Draws a joy plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object scR object
#' @param features.plot Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by fetch.data)
#' @param ident.include Which classes to include in the plot (default is all)
#' @param id.use which id to use (default, cluster id.)
#' @param nCol Number of columns if multiple plots are displayed
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param size.x.use X axis title font size
#' @param size.y.use Y axis title font size
#' @param size.title.use Main title font size
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param x.lab.rot Rotate x-axis labels
#' @param y.lab.rot Rotate y-axis labels
#' @param legend.position Position the legend for the plot
#' @param single.legend Consolidate legend the legend for all plots
#' @param remove.legend Remove the legend from the plot
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param return.plotlist Return the list of individual plots instead of compiled plot.
#' @param \dots additional parameters to pass to FetchData (for example, use.imputed, use.scaled, use.raw)
#'
#' @import ggplot2
#' @importFrom cowplot get_legend
#' @importFrom ggjoy geom_joy theme_joy
#' @importFrom cowplot plot_grid
#'
#' @return By default, no return, only graphical output. If do.return=TRUE,
#' returns a list of ggplot objects.
#'
#' @export
#'
#' @examples
#' JoyPlot(object = pbmc_small, features.plot = 'PC1')
#'

JoyPlot <- function(
  object,
  features.plot,
  ident.include = NULL,
  id.use = NULL,
  nCol = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  size.x.use = 16,
  size.y.use = 16,
  size.title.use = 20,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  x.lab.rot = FALSE,
  y.lab.rot = FALSE,
  legend.position = "right",
  single.legend = TRUE,
  remove.legend = FALSE,
  do.return = FALSE,
  return.plotlist = FALSE,
  ...
) {
  if (is.null(nCol)) {
    if (length(features.plot) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(features.plot), 3)
    }
  }
  data.use <- data.frame(fetch.data(object = object, vars.all = features.plot, ...),
                         check.names = F)
  require(ggjoy)
  if (is.null(id.use)){
    if ("m" %in% colnames(object@data.info)){
      id.use = "m"
    } else {
      id.use = "orig"
    }
  } 
  
  if (!(id.use %in% colnames(object@data.info))){
     stop("Error : id.use is not included in data.info")
  }
  
  if (is.null(ident.include)) {
    cells.to.include <- object@cell.names
  } else {
    cells.to.include <- which.cells(object = object, ident = ident.include, id=id.use)
    if (length(cells.include) == 0) stop("Erorr: must specify valid ident.include and id.use")
  }
  
  data.use <- data.use[cells.to.include, ,drop = FALSE]
  
  if (!is.null(group.by)) {
    ident.use <- as.factor(fetch.data(
      object = object,
      vars.all = group.by
    )[cells.to.include, 1])
  } else {
    ident.use <- object@ident[cells.to.include]
  }
  gene.names <- colnames(data.use)[colnames(data.use) %in% rownames(object@data)]
  if (single.legend) {
    remove.legend <- TRUE
  }
  if (same.y.lims && is.null(y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(
    X = features.plot,
    FUN = function(x) {
      return(SingleJoyPlot(
        feature = x,
        data = data.use[, x, drop = FALSE],
        cell.ident = ident.use,
        do.sort = do.sort, y.max = y.max,
        size.x.use = size.x.use,
        size.y.use = size.y.use,
        size.title.use = size.title.use,
        cols.use = cols.use,
        gene.names = gene.names,
        y.log = y.log,
        x.lab.rot = x.lab.rot,
        y.lab.rot = y.lab.rot,
        legend.position = legend.position,
        remove.legend = remove.legend
      ))
    }
  )
  if (length(features.plot) > 1) {
    plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
    if (single.legend && !remove.legend) {
      legend <- get_legend(
        plot = plots[[1]] + theme(legend.position = legend.position)
      )
      if (legend.position == "bottom") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          ncol = 1,
          rel_heights = c(1, .2)
        )
      } else if (legend.position == "right") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          rel_widths = c(3, .3)
        )
      } else {
        warning("Shared legends must be at the bottom or right of the plot")
      }
    }
  } else {
    plots.combined <- plots[[1]]
  }
  if (do.return) {
    if (return.plotlist) {
      return(plots)
    } else {
      return(plots.combined)
    }
  } else {
    if (length(x = plots.combined) > 1) {
      plots.combined
    }
    else {
      invisible(x = lapply(X = plots.combined, FUN = print))
    }
  }
}

# Plot a single feature on a joy plot
#
# @param feature Feature to plot
# @param data Data to plot
# @param cell.ident Idents to use
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param size.title.use Main title font size
# @param cols.use Colors to use for plotting
# @param gene.names
# @param y.log plot Y axis on log scale
# @param x.lab.rot Rotate x-axis labels
# @param y.lab.rot Rotate y-axis labels
# @param legend.position Position the legend for the plot
# @param remove.legend Remove the legend from the plot
#
# @return A ggplot-based violin plot
#
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggjoy geom_joy theme_joy
#
SingleJoyPlot <- function(
  feature,
  data,
  cell.ident,
  do.sort,
  y.max,
  size.x.use,
  size.y.use,
  size.title.use,
  cols.use,
  gene.names,
  y.log,
  x.lab.rot,
  y.lab.rot,
  legend.position,
  remove.legend
) {
  set.seed(seed = 42)
  feature.name <- colnames(data)
  colnames(data) <- "feature"
  feature <- "feature"
  data$ident <- cell.ident
  if (do.sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(x = tapply(
        X = data[, feature],
        INDEX = data$ident,
        FUN = mean
      ))))
    )
  }
  if (y.log) {
    noise <- rnorm(n = length(data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(data[, feature])) / 100000
  }
  data[, feature] <- data[, feature] + noise
  y.max <- set.ifnull(y.max, max(data[, feature]))
  plot <- ggplot(
    data = data,
    mapping = aes(
      x = feature,
      y = factor(ident)
    )
  ) +
    geom_joy(scale = 4, mapping = aes(fill = factor(x = ident))) + theme_joy() +
    scale_y_discrete(expand = c(0.01, 0)) +   # will generally have to set the `expand` option
    scale_x_continuous(expand = c(0, 0))      # for both axes to remove unneeded padding
  plot <- plot+theme(
    legend.position = legend.position,
    axis.title.x = element_text(
      face = "bold",
      colour = "#990000",
      size = size.x.use
    ),
    axis.title.y = element_text(
      face = "bold",
      colour = "#990000",
      size = size.y.use
    )
  ) +
    guides(fill = guide_legend(title = NULL)) +
    #geom_jitter(height = 0, size = point.size.use) +
    ylab("Identity") +
    nogrid +
    ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"))
  plot <- plot + ggtitle(feature.name)
  if (y.log) {
    plot <- plot + scale_x_log10()
  } else {
    #plot <- plot + xlim(min(data[, feature]), y.max)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + xlab(label = "Log Expression level")
    } else {
      plot <- plot + xlab(label = "Expression level")
    }
  } else {
    plot <- plot + xlab(label = "")
  }
  if (! is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}

######################## HEATMAP #####################################
#' Gene expression heatmap (a la Seurat)
#'
#' Draws a heatmap of single cell gene expression using ggplot2.
#'
#' @param object scR object
#' @param data.use Option to pass in data to use in the heatmap. Default will pick from
#' object@data 
#' @param use.scaled Whether to scale the data (default FALSE)
#' @param cells.use Cells to include in the heatmap (default is all cells)
#' @param genes.use Genes to include in the heatmap (ordered)
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped)
#' @param group.by Groups cells by this variable. Default is object@@ident
#' @param draw.line Draw vertical lines delineating different groups
#' @param col.low Color for lowest expression value
#' @param col.mid Color for mid expression value
#' @param col.high Color for highest expression value
#' @param slim.col.label display only the identity class name once for each group
#' @param remove.key Removes the color key from the plot.
#' @param rotate.key Rotate color scale horizantally
#' @param title Title for plot
#' @param cex.col Controls size of column labels (cells)
#' @param cex.row Controls size of row labels (genes)
#' @param group.label.loc Place group labels on bottom or top of plot.
#' @param group.label.rot Whether to rotate the group label.
#' @param group.cex Size of group label text
#' @param group.spacing Controls amount of space between columns.
### @param assay.type to plot heatmap for (default is RNA)
#' @param do.plot Whether to display the plot.
#'
#' @return Returns a ggplot2 plot object
#'
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#'
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  data.use = NULL,
  use.scaled = FALSE,
  use.centered = FALSE,
  cells.use = NULL,
  genes.use = NULL,
  disp.min = NULL,
  disp.max = NULL,
  group.by = "ident",
  draw.line = TRUE,
  col.low = "#FF00FF",
  col.mid = "#000000",
  col.high = "#FFFF00",
  slim.col.label = FALSE,
  remove.key = FALSE,
  rotate.key = FALSE,
  title = NULL,
  cex.col = 10,
  cex.row = 10,
  midpoint.use=NULL,
  group.label.loc = "bottom",
  group.label.rot = FALSE,
  group.cex = 15,
  group.spacing = 0.15,
  do.plot = TRUE,
  max.cells.per.ident=500,
  downsample.perc = 1,
  ident.order=NULL,
  gene.ident = NULL,
  do.dendro.row=FALSE,
  label.repel = FALSE
) {
  
  
  genes.use = ainb(a=genes.use, b=rownames(object@data))
  if(length(genes.use) ==0) stop("ERROR: Specify valid genes.use parameter")
  
  if (is.null(data.use)) {
    # note: data.use should have cells as column names, genes as row names
    cells.use <- set.ifnull(cells.use, object@cell.names)
    cells.use <- intersect(cells.use, colnames(object@data))
    if (length(cells.use) == 0) {
      stop("No cells given to cells.use present in object")
    }
    data.use = object@data[genes.use, cells.use]
  }

  if (is.null(group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  } else {
    cells.ident <- factor(fetch.data(
      object = object,
      cells.use = cells.use,
      vars.all = group.by
    )[, 1])
    names(cells.ident) <- cells.use
  }
  cells.ident <- factor(
    cells.ident,
    labels = intersect(levels(cells.ident), cells.ident)
  )
  
  if (max.cells.per.ident < max(table(cells.ident))){
    if (downsample.perc == 1){
      cells.subsample = c()
      for (i in levels(cells.ident)){
        cells.in.clust = names(cells.ident)[cells.ident==i]
        cells.in.clust = sample(cells.in.clust, min(max.cells.per.ident, length(cells.in.clust)))
        cells.in.clust = sample(cells.in.clust)
        cells.subsample = c(cells.subsample, cells.in.clust)
      }
      cells.ident = cells.ident[cells.subsample]
      cells.use = cells.subsample
    } else {
      cells.subsample = subsample_by_ident(cells.ident, downsample.perc)
      cells.ident = cells.ident[cells.subsample]
      cells.use = cells.subsample
    }
  
  }
  
  if (!is.null(ident.order)){
    if (sum(ident.order %in% levels(cells.ident))){
      cells.ident = factor(cells.ident, levels = ident.order)
      cells.use2 = c()
      for (i in ident.order){
        cells.use2 = c(cells.use2, names(cells.ident)[cells.ident == i])
      }
      cells.use = cells.use2
      cells.ident = cells.ident[cells.use]
    } else {
      print("Error : ident.order does not match with the data")
    }
    
  }
  
  data.use = as.matrix(data.use[,cells.use])
  if (use.scaled & use.centered) {
    print("z-scoring the genes")
    data.use <- t(scale(t(data.use)))
  }
  
  if (use.scaled & !use.centered){
    data.use = t(scale(t(data.use), center=FALSE, scale=rowMeans(data.use)))
  }
  
  
  if (do.dendro.row){
    clust.tree = hclust(dist(data.use))
    genes.use = rownames(data.use)[clust.tree$order]
    data.use = data.use[genes.use,]
  }
  
  
    #if (disp.max==2.5) disp.max = 10;
  if (is.null(disp.min)){
    disp.min = quantile(data.use,0.2)
  }
  if (is.null(disp.max)){
    disp.max = quantile(data.use,0.9)
  }
  data.use <- minmax(data = data.use, min = disp.min, max = disp.max)
  
  if (is.null(midpoint.use)){
    midpoint.use = quantile(data.use,0.7)
  }
  data.use <- as.data.frame(t(data.use))
  data.use$cell <- rownames(data.use)
  colnames(data.use) <- make.unique(names = colnames(data.use))
  data.use %>% melt(id.vars = "cell") -> data.use
  names(data.use)[names(data.use) == 'variable'] <- 'gene'
  names(data.use)[names(data.use) == 'value'] <- 'expression'
  data.use$ident <- cells.ident[data.use$cell]
  

  breaks <- seq(
    from = min(data.use$expression),
    to = max(data.use$expression),
    length = length(PurpleAndYellow()) + 1
  )
  
  if (!is.null(gene.ident)){
    names(gene.ident) = genes.use
    data.use$gene.ident = gene.ident[data.use$gene]
    data.use$gene.ident = factor(data.use$gene.ident, levels = unique(data.use$gene.ident))
  }
  
  data.use$gene = factor(data.use$gene, levels = rev(levels(data.use$gene)))
  
  # For geom_repel
  data.use$label = ""
  first.cell = data.use[1,]$cell
  data.use[data.use$cell == first.cell,"label"] = as.character(data.use[data.use$cell == first.cell,"gene"])

  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  } else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  #print(head(data.use))
  print(head(data.use))
  hmap <- ggplot(
    data = data.use,
    mapping = aes(x = cell, y = gene, fill = expression)
  ) +
    geom_tile() +
    #scale_fill_viridis() +
    scale_fill_gradient2(
      low = col.low,
      mid = col.mid,
      high = col.high,
      midpoint = midpoint.use,
      name= "Expression",
      guide = guide_colorbar(
        direction = key.direction,
        title.position = key.title.pos
      )
    ) +
    #scale_y_discrete(position = "right", labels = levels(data.use$gene)) +
    theme(
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(size = group.cex),
      axis.text.y = element_text(size = cex.row),
      axis.text.x = element_text(size = cex.col),
      axis.title.x = element_blank()
    )
  print("Initial rendering done")
  if (slim.col.label) {
    hmap <- hmap +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  } else {
    hmap <- hmap + theme(axis.text.x = element_text(angle = 90))
  }
  if (! is.null(group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
    } else {
      switch <- 'x'
    }
    if (is.null(gene.ident)){
      hmap <- hmap +
        facet_grid(~ident,
                   drop = TRUE,
                   space = "free",
                   scales = "free",
                   switch = switch) +
        theme_bw() + theme(panel.spacing = unit(0.1, "lines"), strip.text.y = element_blank(), axis.text.x = element_blank())
      
    } else {
      hmap <- hmap +
        facet_grid(gene.ident~ident,
                   drop = TRUE,
                   space = "free",
                   scales = "free",
                   switch = switch) +
        theme_bw() + theme(panel.spacing = unit(0.1, "lines"), strip.text.y = element_blank(),axis.text.x = element_blank())
    }

      #scale_x_discrete(expand = c(0, 0), drop = TRUE) + theme(panel.border=element_rect(colour="black",size=1))
    if (draw.line) {
      panel.spacing <- unit(group.spacing, units = 'lines')
      # heatmap <- heatmap + theme(strip.background = element_blank(), panel.spacing = unit(group.spacing, "lines"))
    } else {
      panel.spacing <- unit(x = 0, units = 'lines')
      #
    }
    #heatmap <- heatmap +
    #  theme(strip.background = element_blank(), panel.spacing = panel.spacing)
    if (group.label.rot) {
      hmap <- hmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key) {
    hmap <- hmap + theme(legend.position = "none")
  }
  if (! is.null(x = title)) {
    hmap <- hmap + labs(title = title)
  }
  
  if (!label.repel){
    if (do.plot) {
      hmap
    } else {
      return(hmap)
    }
  } else {
    
    if (do.plot) {
       hmap + geom_text_repel(aes(label=label))
    } else {
      return(hmap + geom_text_repel(aes(label=label)))
    }
    
  }
  
  
  
}

#' Plotting heatmap from data matrix
DataHeatmap <- function(data = NULL,
                        xlab.use = "x",
                        ylab.use = "y",
                        val.lab = "value",
                        col.low = "#FF00FF",
                        col.mid = "#000000",
                        col.high = "#FFFF00",
                        col.tile.border = "white",
                        quantile.thresh = c(0.1,0.95),
                        min.break = NULL,
                        max.break = NULL,
                        cex.col = 10,
                        cex.row = 10,
                        group.label.loc = "bottom",
                        group.label.rot = FALSE,
                        group.cex = 15,
                        group.spacing = 0.15,
                        xlab.angle=90,
                        do.plot=TRUE,
                        unname=FALSE,
                        genes.label=NULL,
                        do.repel=FALSE, 
                        return.plot=FALSE){
  
  
  data = as.matrix(data)
  if (unname){
    rownames(data) = NULL; colnames(data) = NULL
  }
  data = minmax(data, quantile(data,quantile.thresh[1]), quantile(data,quantile.thresh[2]))
  data = as.data.frame(data)
  data$x = rownames(data)
  data %>% melt(id.vars = "x") -> data.use
  
  names(data.use)[names(data.use) == 'x'] <- xlab.use
  names(data.use)[names(data.use) == 'variable'] <- ylab.use
  names(data.use)[names(data.use) == 'value'] <- val.lab

  if (is.null(min.break)) min.break = min(data.use[,val.lab])
  if (is.null(max.break)) max.break = max(data.use[,val.lab])
  
  print(min.break)
  print(max.break)
  breaks <- seq(
    from = min.break,
    to = max.break,
    length = length(PurpleAndYellow()) + 1
  )
  
  data.use[,xlab.use] = factor(data.use[,xlab.use], levels = rownames(data))
  data.use[,ylab.use] = factor(data.use[,ylab.use], levels = rev(setdiff(colnames(data),"x")))
  
  key.direction <- "vertical"
  key.title.pos <- "left"
  
  if (!do.repel){
    heatmap <- ggplot(
      data = data.use,
      mapping = aes_string(x = xlab.use, y = ylab.use, fill = val.lab)
    ) + geom_tile(color=col.tile.border) +
      scale_fill_gradient2(
        low = col.low,
        mid = col.mid,
        high = col.high,
        name= val.lab,
        guide = guide_colorbar(
          direction = key.direction,
          title.position = key.title.pos
        )) + 
      scale_y_discrete(position = "left", labels = levels(data.use[,ylab.use])) +
      theme(
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = group.cex),
        axis.text.y = element_text(size = cex.row, face="italic"),
        axis.text.x = element_text(size = cex.col, face="italic"),
        axis.title.x = element_blank()
      ) + theme(axis.text.x=element_text(size=12, face="italic", angle=xlab.angle, hjust=1))
    
    if (do.plot){
      heatmap
    }
  } else {
    if (return.plot) return(heatmap)
    
    heatmap <- ggplot(
      data = data.use,
      mapping = aes_string(x = xlab.use, y = ylab.use, fill = val.lab)
    ) + geom_tile(color=col.tile.border) +
      scale_x_discrete(expand = c(0, 0), drop = TRUE)  + 
      scale_fill_gradient2(
        low = col.low,
        mid = col.mid,
        high = col.high,
        name= val.lab,
        guide = guide_colorbar(
          direction = key.direction,
          title.position = key.title.pos
        )) + 
      scale_y_discrete(position = "right", labels = levels(data.use[,ylab.use])) +
      theme(axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.text.x = element_text(size = group.cex),
            axis.text.x = element_text(size = cex.col, face="italic"),
            legend.justification = "left",
            plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5, "pt"))
    library(ggrepel)
    df1 = data.frame(y=1:(ncol(data)-1), gene = rev(colnames(data)[-ncol(data)]))
    df1$gene = as.character(df1$gene)
    if (!is.null(genes.label)) df1$gene[!(df1$gene %in% genes.label)] = ""
    axis1 <- ggplot(df1, aes(x=0,y=y)) + 
      geom_text_repel(aes(label=gene), min.segment.length = grid::unit(1, "pt")) + scale_x_continuous(limits = c(-3, 0), expand = c(0, 0),
                             breaks = NULL, labels = NULL, name = NULL) +
      scale_y_continuous(limits = c(0.5, ncol(data)+0.5 ), expand = c(0, 0),
                         breaks = NULL, labels = NULL, name = NULL) +
      theme(panel.background = element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
    
    aligned <- cowplot::align_plots(axis1,heatmap + theme(legend.position = "none"), align = "h", axis = "tb")
    aligned <- append(aligned, list(get_legend(heatmap)))
    if (do.plot) plot_grid(plotlist = aligned, nrow = 1, rel_widths = c(3, 7, 4))
    
  }
  
 
 
  return(heatmap)
  
}

######################## DIAGNOSTICS ###################################
#' Get/Plot Silhouette coefficients for a particular clustering
GetSilhouette = function(object,
                         reduction.use = "raw",
                         ident.use = "m",
                         genes.use=NULL,
                         dims.use=c(1:40),
                         do.scale=FALSE,
                         cells.use=NULL,
                         max.cells.per.ident=NULL,
                         do.plot=TRUE,
                         return.avg = FALSE,
                         ...){
  
  
  if (is.null(genes.use)){
    if (length(object@var.genes) > 0){
      genes.use = object@var.genes
    } else {
      genes.use = rownames(object@data)
    }
  }
  genes.use = ainb(genes.use, rownames(object@data))
  cells.use = set.ifnull(cells.use, colnames(object@data))
  if (!is.null(max.cells.per.ident)){
    orig_ident = object@ident
    object = set.all.ident(object, id = ident.use)
    cells.use=SampleCellsPerIdent(object, max.cells.per.ident = max.cells.per.ident,id=ident.use)
    object@ident=orig_ident
    
  }
  
  if (reduction.use == "raw"){
    data.use = object@data[genes.use, cells.use]
    if (do.scale){ data.use = t(scale(t(data.use)))}
    data.use=t(data.use)
  } else if(reduction.use == "regress"){
    data.use = object@regress.data[genes.use, cells.use]
    if (do.scale){ data.use = t(scale(t(data.use)))}
    data.use=t(data.use)
  } else {
    data.use = GetCellEmbeddings(object = object,
                                 reduction.use = reduction.use, 
                                 dims.use = dims.use, 
                                 cells.use = cells.use)

    print(paste0("Using data in ", reduction.use, " space across dimensions ", min(dims.use), " - ", max(dims.use)))
  }
  data.use=as.matrix(data.use)
  if (is.null(data.use)) stop("Invalid Data Space")
  
  orig_ident = object@ident
  object = set.all.ident(object, id = ident.use)
  clust.labels = as.numeric(object@ident[cells.use]); names(clust.labels) = cells.use
  object@ident=orig_ident
  
  #dissE = vectorized_pdist(A = data.use)
  dissE = fastPdist2(data.use, data.use)
  sk <- silhouette(clust.labels, dmatrix = dissE)
  
  
  if (do.plot){
    df = data.frame(cluster=sk[,1], sil_width=sk[,3])
    df_summary = summarySE(data=df,measurevar="sil_width", groupvars="cluster")
    p <- ggplot(df_summary, aes(x=cluster, y=sil_width)) + 
      geom_bar(stat="identity") + 
      geom_errorbar(aes(ymin=sil_width-ci, ymax=sil_width+ci), width=.2) +
      xlab("Cluster ID") + ylab("Average Silhouette Width") + ggtitle(paste0("Median Silhouette Width = ", round(summary(sk)$avg.width,2)))
    print(p)
  }
  
  if (return.avg){
    return(df_summary)
  } else {
    return(sk)
  }
  
  
}


######################## NEW PLOTTING UTILITIES #######################
#remove legend title
no.legend.title <- theme(legend.title = element_blank())

#set legend text
SetLegendTextGG <- function(x = 12, y = "bold") {
  return(theme(legend.text = element_text(size = x, face = y)))
}

#set legend point size
SetLegendPointsGG <- function(x = 6) {
  return(guides(colour = guide_legend(override.aes = list(size = x))))
}

#set x axis features
SetXAxisGG <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.x = element_text(face = z, colour = y, size = x),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}

#set y axis features
SetYAxisGG <- function(x = 16, y = "#990000", z = "bold", x2 = 12) {
  return(theme(
    axis.title.y = element_text(face = z, colour = y, size = x),
    axis.text.y = element_text(angle = 90, vjust = 0.5, size = x2)
  ))
}


#' Dark Theme
#'
#' Add a dark theme to ggplot objects
#'
#' @param ... Extra parameters to be passed to theme()
#' @import ggplot2
#' @return A ggplot2 theme object
#' @seealso \code{theme}
#' @import ggplot2
#' @export
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + DarkTheme(legend.position = 'none')
#'
DarkTheme <- function(...) {
  #   Some constants for easier changing in the future
  black.background <- element_rect(fill = 'black')
  black.background.no.border <- element_rect(fill = 'black', size = 0)
  font.margin <- 4
  white.text <- element_text(
    colour = 'white'
  )
  # margin = margin(
  #   t = font.margin,
  #   r = font.margin,
  #   b = font.margin,
  #   l = font.margin
  # )
  white.line <- element_line(colour = 'white', size = 1)
  no.line <- element_line(size = 0)
  #   Create the dark theme
  dark.theme <- theme(
    #   Set background colors
    plot.background = black.background,
    panel.background = black.background,
    legend.background = black.background,
    legend.box.background = black.background.no.border,
    legend.key = black.background.no.border,
    #   Set text colors
    plot.title = white.text,
    plot.subtitle = white.text,
    axis.title = white.text,
    axis.text = white.text,
    legend.title = white.text,
    legend.text = white.text,
    #   Set line colors
    axis.line.x = white.line,
    axis.line.y = white.line,
    panel.grid = no.line,
    panel.grid.minor = no.line,
    #   Make this a complete theme and validate it
    complete = TRUE,
    validate = TRUE,
    #   Extra parameters
    ...
  )
  return(dark.theme)
}

#' Hover Locator
#'
#' Get quick information from a scatterplot by hovering over points
#'
#' @param plot A ggplot2 plot
#' @param data.plot The oridinal data that went into the ggplot2 plot
#' @param features.info An optional dataframe or matrix of extra information to be displayed on hover
#' @param dark.theme Plot using a dark theme?
#' @param ... Extra parameters to be passed to plotly::layout
#'
#' @seealso \code{plotly::layout}
#' @seealso \code{ggplot2::ggplot_build}
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' HoverLocator(plot = p, data.plot = df)
#' }
#'
HoverLocator <- function(
  plot,
  data.plot,
  features.info = NULL,
  dark.theme = FALSE,
  ...
) {
  #   Use GGpointToBase because we already have ggplot objects
  #   with colors (which are annoying in plotly)
  plot.build <- GGpointToBase(plot = plot, do.plot = FALSE)
  rownames(x = plot.build) <- rownames(data.plot)
  #   Reset the names to 'x' and 'y'
  names(x = plot.build) <- c(
    'x',
    'y',
    names(x = plot.build)[3:length(x = plot.build)]
  )
  #   Add the names we're looking for (eg. cell name, gene name)
  if (is.null(x = features.info)) {
    plot.build$feature <- rownames(x = data.plot)
  } else {
    info <- apply(
      X = features.info,
      MARGIN = 1,
      FUN = function(x, names) {
        return(paste0(names, ': ', x, collapse = '<br>'))
      },
      names = colnames(x = features.info)
    )
    data.info <- data.frame(
      feature = paste(rownames(x = features.info), info, sep = '<br>'),
      row.names = rownames(x = features.info)
    )
    plot.build <- merge(x = plot.build, y = data.info, by = 0)
  }
  #   Set up axis labels here
  #   Also, a bunch of stuff to get axis lines done properly
  xaxis <- list(
    title = names(x = data.plot)[1],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  yaxis <- list(
    title = names(x = data.plot)[2],
    showgrid = FALSE,
    zeroline = FALSE,
    showline = TRUE
  )
  #   Check for dark theme
  if (dark.theme) {
    title <- list(color = 'white')
    xaxis <- c(xaxis, color = 'white')
    yaxis <- c(yaxis, color = 'white')
    plotbg <- 'black'
  } else {
    title = list(color = 'black')
    plotbg = 'white'
  }
  #   Start plotly and pipe it into layout for axis modifications
  #   The `~' means pull from the data passed (this is why we reset the names)
  #   Use I() to get plotly to accept the colors from the data as is
  #   Set hoverinfo to 'text' to override the default hover information
  #   rather than append to it
  plotly::plot_ly(
    data = plot.build,
    x = ~x,
    y = ~y,
    type = 'scatter',
    mode = 'markers',
    color = ~I(color),
    hoverinfo = 'text',
    text = ~feature
  ) %>% plotly::layout(
    xaxis = xaxis,
    yaxis = yaxis,
    titlefont = title,
    paper_bgcolor = plotbg,
    plot_bgcolor = plotbg,
    ...
  )
}

# Convert a ggplot2 scatterplot to base R graphics
#
# @param plot A ggplot2 scatterplot
# @param do.plot Create the plot with base R graphics
# @param ... Extra parameters passed to PlotBuild
#
# @return A dataframe with the data that created the ggplot2 scatterplot
#
GGpointToBase <- function(plot, do.plot = TRUE, ...) {
  plot.build <- ggplot2::ggplot_build(plot = plot)
  build.data <- plot.build$data[[1]]
  plot.data <- build.data[, c('x', 'y', 'colour', 'shape', 'size')]
  names(x = plot.data) <- c(
    plot.build$plot$labels$x,
    plot.build$plot$labels$y,
    'color',
    'pch',
    'cex'
  )
  if (do.plot) {
    PlotBuild(plot.data = plot.data, ...)
  }
  return(plot.data)
}

# Create a scatterplot with data from a ggplot2 scatterplot
#
# @param plot.data The original ggplot2 scatterplot data
# This is taken from ggplot2::ggplot_build
# @param dark.theme Plot using a dark theme
# @param smooth Use a smooth scatterplot instead of a standard scatterplot
# @param ... Extra parameters passed to graphics::plot or graphics::smoothScatter
#
#' @importFrom graphics axis
#
PlotBuild <- function(plot.data, dark.theme = FALSE, smooth = FALSE, ...) {
  #   Do we use a smooth scatterplot?
  #   Take advantage of functions as first class objects
  #   to dynamically choose normal vs smooth scatterplot
  if (smooth) {
    myplot <- smoothScatter
  } else {
    myplot <- plot
  }
  if (dark.theme) {
    par(bg = 'black')
    axes = FALSE
    col.lab = 'white'
  } else {
    axes = 'TRUE'
    col.lab = 'black'
  }
  myplot(
    plot.data[, c(1, 2)],
    col = plot.data$color,
    pch = plot.data$pch,
    cex = vapply(
      X = plot.data$cex,
      FUN = function(x) {
        return(max(x / 2, 0.5))
      },
      FUN.VALUE = numeric(1)
    ),
    axes = axes,
    col.lab = col.lab,
    col.main = col.lab,
    ...
  )
  if (dark.theme) {
    axis(
      side = 1,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
    axis(
      side = 2,
      at = NULL,
      labels = TRUE,
      col.axis = col.lab,
      col = col.lab
    )
  }
}



#' Feature Locator
#'
#' Select points on a scatterplot and get information about them
#'
#' @param plot A ggplot2 plot
#' @param data.plot The oridinal data that went into the ggplot2 plot
#' @param ... Extra parameters, such as dark.theme, recolor, or smooth for using a dark theme,
#' recoloring based on selected cells, or using a smooth scatterplot, respectively
#'
#' @return The names of the points selected
#'
#' @seealso \code{locator}
#' @seealso \code{ggplot2::ggplot_build}
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' FeatureLocator(plot = p, data.plot = df)
#' }
#'
FeatureLocator <- function(plot, data.plot, ...) {
  points.located <- PointLocator(plot = plot, ...)
  #   The rownames for points.located correspond to the row indecies
  #   of data.plot thanks to the way the ggplot object was made
  selected <- data.plot[as.numeric(x = rownames(x = points.located)), ]
  return(rownames(x = selected))
}

# Locate points on a plot and return them
#
# @param plot A ggplot2 plot
# @param recolor Do we recolor the plot to highlight selected points?
# @param dark.theme Plot using a dark theme
# @param ... Exptra parameters to PlotBuild
#
# @return A dataframe of x and y coordinates for points selected
#
#' @importFrom graphics locator
#
PointLocator <- function(plot, recolor=TRUE, dark.theme = FALSE, ...) {
  #   Convert the ggplot object to a data.frame
  plot.data <- GGpointToBase(plot = plot, dark.theme = dark.theme, ...)
  npoints <- nrow(x = plot.data)
  cat("Click around the cluster of points you wish to select\n")
  cat("ie. select the vertecies of a shape around the cluster you\n")
  cat("are interested in. Press <Esc> when finished (right click for R-terminal users)\n\n")
  polygon <- locator(n = npoints, type = 'l')
  polygon <- data.frame(polygon)
  #   pnt.in.poly returns a data.frame of points
  points.all <- SDMTools::pnt.in.poly(
    pnts = plot.data[, c(1, 2)],
    poly.pnts = polygon
  )
  #   Find the located points
  points.located <- points.all[which(x = points.all$pip == 1), ]
  #   If we're recoloring, do the recolor
  if(recolor) {
    if (dark.theme) {
      no = 'white'
    } else {
      no = 'black'
    }
    points.all$color <- ifelse(test = points.all$pip == 1, yes = 'red', no = no)
    plot.data$color <- points.all$color
    PlotBuild(plot.data = plot.data, dark.theme = dark.theme, ...)
  }
  return(points.located[, c(1, 2)])
}

#' Create a custom color palette
#'
#' Creates a custom color palette based on low, middle, and high color values
#'
#' @param low low color
#' @param high high color
#' @param mid middle color. Optional.
#' @param k number of steps (colors levels) to include between low and high values
#'
#' @return A color palette for plotting
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @export
#'
#' @examples
#' myPalette <- CustomPalette()
#' myPalette
#'
CustomPalette <- function(
  low = "white",
  high = "red",
  mid = NULL,
  k = 50
) {
  low <- col2rgb(col = low) / 255
  high <- col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2],len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(rgb(red = r, green = g, blue = b))
}

#' A black and white color palette
#'
#' @param ... Extra parameters to CustomPalette
#'
#' @return A color palette
#'
#' @seealso \code{CustomPalette}
#'
#' @export
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlackAndWhite())
#'
BlackAndWhite <- function(...) {
  return(CustomPalette(low = "white", high="black", ...))
}

#' A purple and yellow color palette
#'
#' @param ... Extra parameters to CustomPalette
#'
#' @return A color palette
#'
#' @seealso \code{CustomPalette}
#'
#' @export
#'
#' @examples
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' plot(df, col = BlackAndWhite())
#'
PurpleAndYellow <- function(...) {
  return(CustomPalette(low = "magenta", high = "yellow", mid = "black", ...))
}


################################# OLD ############################

setGeneric("feature.plot", function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors,pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL, cex.main=1.2, scale.range=NULL, do.scale=FALSE, xlab.use=NULL, ylab.use=NULL, plot.kcenters=FALSE) standardGeneric("feature.plot"))
setMethod("feature.plot", "scR", 
          function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors,pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL, cex.main=1.2, scale.range=NULL, do.scale=FALSE, xlab.use=NULL, ylab.use=NULL, plot.kcenters=FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)==1) nCol=1
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }         
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            data.plot = NULL
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
              if (plot.kcenters){
                k.centers=t(sapply(levels(object@ident),function(x) apply(object@pca.rot[which.cells(object,x),],2,mean)))
              }
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
              if (plot.kcenters){
                k.centers=t(sapply(levels(object@ident),function(x) apply(object@tsne.rot[which.cells(object,x),],2,mean)))
              }
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
              if (plot.kcenters){
                k.centers=t(sapply(levels(object@ident),function(x) apply(object@ica.rot[which.cells(object,x),],2,mean)))
              }
            }
            
            if (is.null(data.plot)){
              
              data.plot <- as.data.frame(GetCellEmbeddings(
                object = object,
                reduction.use = reduction.use,
                dims.use = c(pc.1, pc.2),
                cells.use = cells.use
              ))
              dim.code = GetDimCode(reduction.use)
              if (plot.kcenters){
                k.centers=t(sapply(levels(object@ident),function(x) apply(data.plot[which.cells(object,x),],2,mean)))
              }
            }
            
            ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use
            xlab.use=set.ifnull(xlab.use,paste(dim.code,pc.1,sep="")); ylab.use=set.ifnull(ylab.use,paste(dim.code,pc.2,sep=""))
            x1 = paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use, use.raw=use.raw, use.count=use.count)))
            if (do.scale){
              data.use = t(scale(t(data.use)))
              data.use = as.data.frame(data.use)
            }
            for(i in features.plot) {
              data.gene=as.numeric(na.omit(data.frame(data.use[i,])))
              
              if (is.null(scale.range)){  
                data.range = c(min(data.gene), max(data.gene)) 
              } else {
                data.range = c(scale.range[1], scale.range[2])
              }
              data.gene[data.gene>=scale.range[2]] = scale.range[2]
              data.gene[data.gene<=scale.range[1]] = scale.range[1]
              
              
              if(data.range[2] > 5){
                b=1
              }else{
                b=(data.range[2] - data.range[1]) /5
              }
              breaks = seq(data.range[1], data.range[2], by=b)
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = breaks)))
              data.cut[is.na(data.cut)]=1
              data.col=rev(cols.use(length(breaks)))[data.cut]
              print(rev(cols.use(length(breaks))))
              plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=xlab.use,ylab=ylab.use, xlim=xlim, ylim=ylim, cex.main=cex.main)
              if (plot.kcenters){
                points(k.centers[,pc.1],k.centers[,pc.2],cex=1.3,col="white",pch=16); 
                text(k.centers[,pc.1],k.centers[,pc.2],levels(object@ident),cex=1)
              }
            }
            rp()
          }
)

setGeneric("feature.plot.3d", function(object,features.plot,pc.1=1,pc.2=2,pc.3=3,cells.use=NULL,pt.size=0.5,cols.use=heat.colors,pch.use=16,reduction.use="pca",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL, zlim=NULL, cex.main=1.2, scale.range=NULL) standardGeneric("feature.plot.3d"))
setMethod("feature.plot.3d", "scR", 
          function(object,features.plot,pc.1=1,pc.2=2,pc.3=3,cells.use=NULL,pt.size=0.5,cols.use=heat.colors,pch.use=16,reduction.use="pca",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL,zlim=NULL, cex.main=1.2,scale.range=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)==1) nCol=1
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }         
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            
            ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep=""); x3=paste(dim.code,pc.3,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]; data.plot$z=data.plot[,x3]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use, use.raw=use.raw, use.count=use.count)))
            
            for(i in features.plot) {
              data.gene=as.numeric(na.omit(data.frame(data.use[i,])))
              
              if (is.null(scale.range)){  
                data.range = c(min(data.gene), max(data.gene)) 
              } else {
                data.range = c(scale.range[1], scale.range[2])
              }
              data.gene[data.gene>=scale.range[2]] = scale.range[2]
              data.gene[data.gene<=scale.range[1]] = scale.range[1]
              
              if(data.range[2] > 5){
                b=1
              }else{
                b=(data.range[2] - data.range[1]) /5
              }
              
              breaks = seq(data.range[1], data.range[2], by=b)
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = breaks)))
              data.cut[is.na(data.cut)]=1
              data.col=rev(cols.use(length(breaks)))[data.cut]
              clear3d("all")
              par3d(windowRect = c(0, 0, 700, 700)) # make the window large
              par3d(zoom = 1.1)
              light3d()	
              aspect3d(1,1,1)
              plot3d(data.plot$x, data.plot$y, data.plot$z, col=data.col,size=pt.size, zlab=x3,xlab=x1,ylab=x2,xlim=xlim, ylim=ylim, zlim=zlim, type="s", box=TRUE, cex.main=cex.main)
            }
            rp()
          }
)


setGeneric("gg.feature", function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors,pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL, cex.main=1.2, max.val=NULL, scale.range=NULL) standardGeneric("gg.feature"))
setMethod("gg.feature", "scR", 
          function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors,pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE, use.count=FALSE, xlim=NULL, ylim=NULL, cex.main=1.2,max.val=NULL, scale.range=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }         
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            
            ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use, use.raw=use.raw, use.count=use.count)))
            data.plot = cbind(data.plot,t(data.use))
            p=list()
            for (i in c(1:length(features.plot)))
              p[[i]] = ggplot(data.plot, aes_string(x1,x2)) + geom_point(aes_string(color=features.plot[[i]])) + geom_point(aes(alpha=1)) + scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                                                                                                                                                                  high="red", space ="Lab" )
            for(i in features.plot) {
              data.gene=as.numeric(na.omit(data.frame(data.use[i,])))
              
              if (is.null(scale.range)){  
                data.range = c(min(data.gene), max(data.gene)) 
              } else {
                data.range = c(scale.range[1], scale.range[2])
              }
              data.gene[data.gene>=scale.range[2]] = scale.range[2]
              data.gene[data.gene<=scale.range[1]] = scale.range[1]
              if (!is.null(max.val)) data.gene[data.gene>=max.val] = max.val
              
              
              if(data.range[2] > 5){
                b=1
              }else{
                b=(data.range[2] - data.range[1]) /5
              }
              breaks = seq(data.range[1], data.range[2], by=b)
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = breaks)))
              data.cut[is.na(data.cut)]=1
              data.col=rev(cols.use(length(breaks)))[data.cut]
              print(rev(cols.use(length(breaks))))
              #plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2, xlim=xlim, ylim=ylim, cex.main=cex.main)
            }
            rp()
          }
)


setGeneric("feature.plot.scale", function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=redblue(10),pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE,use.count=FALSE, show.labels=TRUE, ylim=NULL, xlim=NULL, max.val=NULL, scale.range=NULL) standardGeneric("feature.plot.scale"))
setMethod("feature.plot.scale", "scR", 
          function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=redblue(10),pch.use=16,reduction.use="tsne",nCol=NULL, use.raw=FALSE,use.count=FALSE, show.labels=TRUE, ylim=NULL, xlim=NULL, max.val=NULL,scale.range=NULL ) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }         
            num.row=floor(length(features.plot)/nCol-1e-5)+1
            par(mfrow=c(num.row,nCol))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            print(features.plot)
            data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use, use.raw=use.raw, use.count=use.count)))
            print(1)
            for(i in features.plot) {
              
              # get data for the gene
              data.gene=as.numeric(na.omit(data.frame(data.use[i,])))
              
              if (is.null(scale.range)){  
                data.range = c(min(data.gene), max(data.gene)) 
              } else {
                data.range = c(scale.range[1], scale.range[2])
              }
              
              data.gene[data.gene>=scale.range[2]] = scale.range[2]
              data.gene[data.gene<=scale.range[1]] = scale.range[1]
              if (!is.null(max.val)) data.gene[data.gene>=max.val] = max.val
              
              
              if(data.range[2] > 5){
                b=1
              }else{
                b=(data.range[2] - data.range[1]) / 5
              }
              breaks = seq(data.range[1], data.range[2], by=b)
              
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = breaks)))
              data.cut[is.na(data.cut)]=1
              data.col=rev(cols.use(length(breaks)))[data.cut]
              par("mai"=c(1.5, 1.5, 1.5, 1.5))
              plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2, xlim=xlim, ylim=ylim)
              #plot(data.plot$x,data.plot$y,bg=data.col, col="gray",cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2, xlim=xlim, ylim=ylim)
              
              
              if(data.range[2] > 5){
                b=1
              }else{
                b=(data.range[2] - data.range[1]) / 5
              }
              
              
              
              
              print(1)
              digit=0;
              if (data.range[2] < 2) digit=1
              if(show.labels)
              {
                shape::colorlegend(col=rev(cols.use(length(breaks))), zlim=data.range, zval=breaks,  left=F, posy=c(0.1,0.8), posx=c(0.85,0.88), main = "log(Exp)", digit=digit) #main="Log2(count)",
              }else{
                shape::colorlegend(col=rev(cols.use(length(breaks))), zlim=NULL, zval=NULL, posy=c(0,0.8), main = "log(Exp)", digit=digit) #main="Log2(count)",
              }
              
            }
            par(mfrow=c(1,1))
            rp()
          }
)

setGeneric("plot.cluster", function(object,clusters.plot=NULL,cells.plot=NULL,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=redblue(3),pch.use=16,reduction.use="tsne",nCol=NULL, id.use=NULL, ident.use=NULL, main.use=NULL) standardGeneric("plot.cluster"))
setMethod("plot.cluster", "scR", 
          function(object,clusters.plot=NULL,cells.plot=NULL,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=redblue(3),pch.use=16,reduction.use="tsne",nCol=NULL, id.use=NULL, ident.use=NULL, main.use=NULL) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            if (is.null(nCol)) {
              nCol=2
              if (length(clusters.plot)>6) nCol=3
              if (length(clusters.plot)>9) nCol=4
            }         
            
            if (is.null(clusters.plot)){
              num.row = 1
            } else {
              num.row=floor(length(clusters.plot)/nCol-1e-5)+1
            }
            par(mfrow=c(num.row,nCol))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot[cells.use,]
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            
            id.use = set.ifnull(id.use, "orig")
            ident.use = object@data.info[cells.use,id.use]
            names(ident.use) = cells.use
            #ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use = data.frame()
            if (!is.null(clusters.plot)){
              for (i in clusters.plot){
                vec = as.numeric(ident.use[cells.use] == i)
                data.use=rbind(data.use, vec)
              }
              rownames(data.use) = clusters.plot
            } else {
              data.use=rbind(data.use, as.numeric(names(ident.use) %in% cells.plot))
              rownames(data.use) = "cells"
            }
            colnames(data.use) = cells.use
            
            for(i in rownames(data.use)) {
              data.gene=na.omit(data.frame(data.use[i,]))
              data.cut=as.numeric(as.factor(cut(as.numeric(data.gene),breaks = length(cols.use))))
              data.col=rev(cols.use)[data.cut]
              maintext = paste0(paste0(id.use, "_", i), " (n = ", length(which.cells(object, i, id=id.use)), ")")
              if (!is.null(main.use)) maintext = main.use
              plot(data.plot$x[order(data.gene)],data.plot$y[order(data.gene)],col=data.col[order(data.gene)],cex=pt.size,pch=pch.use,main=maintext,xlab=x1,ylab=x2)
            }
            rp()
          }
)


setGeneric("feature.heatmap", function(object,features.plot,pc.1=1,pc.2=2,idents.use=NULL,pt.size=2,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") standardGeneric("feature.heatmap"))
setMethod("feature.heatmap", "scR", 
          function(object,features.plot,pc.1=1,pc.2=2,idents.use=NULL,pt.size=1,cols.use=rev(heat.colors(10)),pch.use=16,reduction.use="tsne") {
            idents.use=set.ifnull(idents.use,sort(unique(object@ident)))
            dim.code="PC"
            par(mfrow=c(length(features.plot),length(idents.use)))
            if (reduction.use=="pca") {
              data.plot=object@pca.rot
              dim.code="PC"
            }
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot
              dim.code="IC"
            }
            
            
            ident.use=as.factor(object@ident)
            data.plot$ident=ident.use
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
            data.use=data.frame(t(fetch.data(object,features.plot)))
            data.scale=apply(t(data.use),2,function(x)(factor(cut(x,breaks=length(cols.use),labels=1:length(cols.use),ordered=TRUE))))
            data.plot.all=cbind(data.plot,data.scale)
            data.reshape=melt(data.plot.all,id = colnames(data.plot))
            data.reshape=data.reshape[data.reshape$ident%in%idents.use,]
            p <- ggplot(data.reshape, aes(x,y)) + geom_point(aes(colour=reorder(value,1:length(cols.use)),size=pt.size)) + scale_colour_manual(values=cols.use)
            p=p + facet_grid(variable~ident) + scale_size(range = c(pt.size, pt.size))
            p2=p+gg.xax()+gg.yax()+gg.legend.pts(6)+ggplot.legend.text(12)+no.legend.title+theme_bw()+nogrid+theme(legend.title=element_blank())
            print(p2)
          }
)



setGeneric("pca.plot", function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,reduction.use="pca", plot.kcenters=FALSE, kid.use="m", size.kcenters=7) standardGeneric("pca.plot"))
setMethod("pca.plot", "scR", 
          function(object,pc.1=1,pc.2=2,cells.use=NULL,pt.size=3,do.return=FALSE,do.bare=FALSE,cols.use=NULL,reduction.use="pca", plot.kcenters=FALSE, kid.use="m",size.kcenters=7) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code="PC"
            order = sample(c(1:length(cells.use)), replace=FALSE)
            if (reduction.use=="pca") data.plot=object@pca.rot[cells.use,]
            if (reduction.use=="tsne") {
              data.plot=object@tsne.rot[cells.use,]
              dim.code="tsne_"
            }
            if (reduction.use=="ica") {
              data.plot=object@ica.rot[cells.use,]
              dim.code="IC"
            }
            ident.use=as.factor(object@ident[cells.use])
            data.plot$ident=ident.use[order]
            x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
            data.plot$x=data.plot[order,x1]; data.plot$y=data.plot[order,x2]
            data.plot$pt.size=pt.size
            rownames(data.plot) = cells.use
            
            p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident), fill=factor(ident),shape=factor(ident),size=pt.size)) + scale_color_hue(l=55) + 
              scale_shape_manual(values=1:nlevels(data.plot$ident))
            #p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident),size=pt.size)) + scale_color_hue(l=55) 
            if (plot.kcenters){
              orig_id = object@ident
              object = set.all.ident(object, id=kid.use)
              k.centers=t(sapply(levels(object@ident),function(x) apply(data.plot[intersect(which.cells(object,x), rownames(data.plot)),c(1,2)],2,median)))
              object@ident = orig_id
              for (i in 1:nrow(k.centers)){
                p = p + annotate(geom="text", x=k.centers[i,1], y=k.centers[i,2], label=i,
                                 color="black", size=size.kcenters)
              }
            }  
            
            if (!is.null(cols.use)) {
              p=p+scale_colour_manual(values=cols.use) 
            }
            
            p2=p+xlab(x1)+ylab(x2)+scale_size(range = c(pt.size, pt.size))
            p3=p2+gg.xax()+gg.yax()+gg.legend.pts(6)+ggplot.legend.text(12)+no.legend.title+theme_bw()+nogrid
            
            if (do.return) {
              if (do.bare) return(p)
              return(p3)
            }
            
            if (do.bare){ 
              print(p)
              
            }else { print(p3)}
          }
)

#####################HEATMAP FUNCTIONS########################

setGeneric("doHeatMap", function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,ident.order = NULL, col.use=pyCols,do.Rowv=FALSE,row.ident=NULL,rowsep.use=NULL,cex.row.use=0.3, cex.col.use=0.3, window.size=30, do.smooth=FALSE,do.scale=TRUE,...) standardGeneric("doHeatMap"))

setMethod("doHeatMap","scR",
          function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,ident.order = NULL,col.use=pyCols,do.Rowv=FALSE,row.ident=NULL,rowsep.use=NULL,cex.row.use=0.3, cex.col.use=0.3, window.size=30,do.smooth=FALSE,do.scale=TRUE,...) {
            cells.use=set.ifnull(cells.use,object@cell.names)
            genes.use=ainb(genes.use,rownames(object@data))
            cells.use=ainb(cells.use,object@cell.names)
            if (order.by.ident) {
              if (is.null(ident.order)){
                cells.ident=object@ident[cells.use]
                cells.use=cells.use[order(cells.ident)]
              } else {
                cells.ident = c()
                for (i in ident.order){
                  cells.i = which.cells(object, i); 
                  cells.i = ainb(cells.i, cells.use)
                  cells.ident = c(cells.ident, cells.i)
                }
              }
            }
            
            cells.use = names(cells.ident)
            data.use=object@data[rev(genes.use),cells.use]
            #data.use=as.matrix(object@count.data[genes.use,cells.use])
            vline.use=NULL;
            if (draw.line) {
              colsep.use=cumsum(table(cells.ident))
            }
            if (!is.null(row.ident)){
              rowsep.use=cumsum(table(row.ident))
            }
            
            
            lab.col.use = rep("", ncol(data.use))
            
            #lab.col.use = rep("", ncol(data.plot))
            pdf("abc.pdf",w=20,h=20)
            heatmap.2(as.matrix(data.use),Rowv=do.Rowv,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,margins = c(15,5),rowsep=rowsep.use)
            dev.off()
            print(dim(data.use))
            if (max(c(nrow(data.use), ncol(data.use))) < 500){
              if (do.smooth){
                print(paste0("Smoothing rows. Using window size = ", window.size))
                data.use1 = t(caTools::runmean(t(data.use), k=window.size))
                rownames(data.use1) = rownames(data.use)
                colnames(data.use1) = colnames(data.use)
                data.use = data.use1
              }
              if (do.scale) data.use = t(scale(t(data.use)))
              data.use=minmax(data.use,min=disp.min,max=disp.max)
              heatmap.2(data.use,Rowv=do.Rowv,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=rowsep.use,...)
            } else {
              if (do.smooth){
                print(paste0("Smoothing rows. Using window size = ", window.size))
                data.use1 = t(caTools::runmean(t(data.use), k=window.size))
                rownames(data.use1) = rownames(data.use)
                colnames(data.use1) = colnames(data.use)
                data.use = data.use1
              }
              if (do.scale) data.use = t(scale(t(data.use)))
              data.use=minmax(data.use,min=disp.min,max=disp.max)
              dim(data.use)
              pdf("abc.pdf", w=20,h=20)
              heatmap3(data.use,Colv=NA, Rowv=do.Rowv, scale="none", cexCol=cex.col.use, balanceColor = T, cexRow=cex.row.use, labCol=lab.col.use)
              dev.off()
            }
            if (do.return) {
              return(data.use)
            }
          }
)

setGeneric("doCHmap.by.group", function(object,genes.use=NULL,cells.use=NULL, ident.use=NULL,fxn.x=expMean, do.scale=TRUE,
                                        ylab.use=TRUE,xlab.use=TRUE,color.range=NULL,rcut.frac=0.8, ccut.frac=0.8,return.plot=TRUE,
                                        write.clusters=FALSE, clusters.file = "Cluster_genes.txt",leg.title="Intensity",do.cluster=TRUE,...) standardGeneric("doCHmap.by.group"))
setMethod("doCHmap.by.group", "scR", 
          function(object,genes.use=NULL,cells.use=NULL,ident.use=NULL, fxn.x=expMean, do.scale=TRUE, 
                   ylab.use=TRUE,xlab.use=TRUE, color.range=NULL,rcut.frac=0.8, ccut.frac=0.8,return.plot=TRUE,
                   write.clusters=FALSE, clusters.file = "Cluster_genes.txt",leg.title="Intensity",do.cluster=TRUE,...) {
            
            genes.use = set.ifnull(genes.use, rownames(object@data))
            genes.use = ainb(genes.use, rownames(object@data))
            cells.use = set.ifnull(cells.use, colnames(object@data))
            cells.use = ainb(cells.use, colnames(object@data))
            ident.use = set.ifnull(ident.use,object@ident[cells.use])
            ident.levels = levels(ident.use)[levels(ident.use) %in% unique(as.character(ident.use))]
            #Average expression by cluster
            ExpMat = matrix(0, nrow=length(genes.use), ncol = length(ident.levels))
            rownames(ExpMat) = genes.use; colnames(ExpMat) = ident.levels
            
            for (i in ident.levels){
              cells.in.cluster = cells.use[which(ident.use== i)]
              vec.exp = apply(object@data[genes.use, cells.in.cluster], 1, function(x) fxn.x(x)) 
              ExpMat[, i] = vec.exp
            }
            
            if (do.scale){
              ExpMat = t(scale(t(ExpMat)))
              color.range = set.ifnull(color.range, c(-2,2))
              if (do.cluster){
                p = create.gg.hmap.w.barps(hmap.dat=ExpMat, cut.frac.r.h=rcut.frac, cut.frac.c.h=ccut.frac, cap=color.range, no.y.labels=!ylab.use,no.x.labels=!xlab.use,leg.title=leg.title,...)
              } else {
                p = gg.hmap(dat = ExpMat, cap=color.range, no.y.labels=!ylab.use,no.x.labels=!xlab.use,...)
              }
            } else {
              color.range = set.ifnull(color.range, c(0,4))
              if (do.cluster){ 
                p = create.gg.hmap.w.barps(hmap.dat=ExpMat, cut.frac.r.h=rcut.frac, cut.frac.c.h=ccut.frac, cap=color.range, no.y.labels=!ylab.use,no.x.labels=!xlab.use,leg.title=leg.title,...)
              } else {
                p = gg.hmap(dat = ExpMat, cap=color.range, no.y.labels=!ylab.use,no.x.labels=!xlab.use,...)
              }
            }
            
            if (do.cluster & write.clusters){
              file.create(clusters.file)
              genes.ord = genes.use[p$p.hm.l$ord.r]
              genes.cluster = list()
              splits = p$p.hm.l$r.splits
              len.list = length(splits)-1
              #Store cluster genes in a list
              for (i in 1:len.list){
                write(paste0("\n \n Cluster ", i," \n\n"), file=clusters.file, append=TRUE)
                genes.cluster[[i]] = genes.ord[splits[(len.list-i+1)]:splits[(len.list-i+2)]]
                write(genes.cluster[[i]], file=clusters.file, ncolumns=10, append=TRUE, sep="\t")
              }
              
            }
            
            if (!return.plot){
              if(do.cluster){
                print(p$p.hm.l)
              } else { print(p)}
            } else {
              return(p)
            }
            
            
          }
)


#############
######Violin Plots##########################

#' Single cell violin plot (from Seurat)
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object scR object
#' @param features.plot Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param ident.include Which classes to include in the plot (default is all)
#' @param nCol Number of columns if multiple plots are displayed
#' @param do.sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param size.x.use X axis title font size
#' @param size.y.use Y axis title font size
#' @param size.title.use Main title font size
#' @param adjust.use Adjust parameter for geom_violin
#' @param point.size.use Point size for geom_violin
#' @param cols.use Colors to use for plotting
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.log plot Y axis on log scale
#' @param x.lab.rot Rotate x-axis labels
#' @param y.lab.rot Rotate y-axis labels
#' @param legend.position Position the legend for the plot
#' @param single.legend Consolidate legend the legend for all plots
#' @param remove.legend Remove the legend from the plot
#' @param do.return Return a ggplot2 object (default : FALSE)
#' @param return.plotlist Return the list of individual plots instead of compiled plot.
#' @param \dots additional parameters to pass to FetchData (for example, use.imputed, use.scaled, use.raw)
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid get_legend
#'
#' @return By default, no return, only graphical output. If do.return=TRUE,
#' returns a list of ggplot objects.
#'
#' @export
#'
#' @examples
#' VlnPlot(object = pbmc_small, features.plot = 'PC1')
#'
VlnPlot <- function(
  object,
  features.plot,
  ident.include = NULL,
  nCol = NULL,
  do.sort = FALSE,
  y.max = NULL,
  same.y.lims = FALSE,
  size.x.use = 16,
  size.y.use = 16,
  size.title.use = 20,
  adjust.use = 1,
  point.size.use = 1,
  cols.use = NULL,
  group.by = NULL,
  y.log = FALSE,
  x.lab.rot = FALSE,
  y.lab.rot = FALSE,
  legend.position = "right",
  single.legend = TRUE,
  remove.legend = FALSE,
  do.return = FALSE,
  return.plotlist = FALSE,
  do.jitter=FALSE,
  print.frac=TRUE,
  do.boxplot=TRUE,
  use.regress=FALSE,
  ...
) {
  if (is.null(x = nCol)) {
    if (length(x = features.plot) > 9) {
      nCol <- 4
    } else {
      nCol <- min(length(x = features.plot), 3)
    }
  }
  data.use <- data.frame(fetch.data(object = object, vars.all = features.plot, use.regress=use.regress,...),
                         check.names = F)
  if (is.null(x = ident.include)) {
    cells.to.include <- object@cell.names
  } else {
    cells.to.include <- which.cells(object = object, ident.include)
  }
  data.use <- data.use[cells.to.include, ,drop = FALSE]
  if (!is.null(x = group.by)) {
    ident.use <- as.factor(x=fetch.data(
      object = object,
      vars.all = group.by
    )[cells.to.include, 1])
  } else {
    ident.use <- object@ident[cells.to.include]
  }
  
  if (!is.null(x = ident.include)){
    ident.use = factor(ident.use, levels = ident.include)
  }
  
  gene.names <- colnames(x = data.use)[colnames(x = data.use) %in% rownames(x = object@data)]
  if (single.legend) {
    remove.legend <- TRUE
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data.use)
  }
  plots <- lapply(
    X = features.plot,
    FUN = function(x) {
      return(SingleVlnPlot(
        feature = x,
        data = data.use[, x, drop = FALSE],
        cell.ident = ident.use,
        do.sort = do.sort, y.max = y.max,
        size.x.use = size.x.use,
        size.y.use = size.y.use,
        size.title.use = size.title.use,
        adjust.use = adjust.use,
        point.size.use = point.size.use,
        cols.use = cols.use,
        gene.names = gene.names,
        y.log = y.log,
        x.lab.rot = x.lab.rot,
        y.lab.rot = y.lab.rot,
        legend.position = legend.position,
        remove.legend = remove.legend,
        do.jitter = do.jitter,
        do.boxplot = do.boxplot,
        print.frac=print.frac,
        ...
      ))
    }
  )
  if (length(x = features.plot) > 1) {
    plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
    if (single.legend && !remove.legend) {
      legend <- get_legend(
        plot = plots[[1]] + theme(legend.position = legend.position)
      )
      if (legend.position == "bottom") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          ncol = 1,
          rel_heights = c(1, .2)
        )
      } else if (legend.position == "right") {
        plots.combined <- plot_grid(
          plots.combined,
          legend,
          rel_widths = c(3, .3)
        )
      } else {
        warning("Shared legends must be at the bottom or right of the plot")
      }
    }
  } else {
    plots.combined <- plots[[1]]
  }
  if (do.return) {
    if (return.plotlist) {
      return(plots)
    } else {
      return(plots.combined)
    }
  } else {
    if (length(x = plots.combined) > 1) {
      plots.combined
    }
    else {
      invisible(x = lapply(X = plots.combined, FUN = print))
    }
  }
}

# Plot a single feature on a violin plot
#
# @param feature Feature to plot
# @param data Data to plot
# @param cell.ident Idents to use
# @param do.sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param size.x.use X axis title font size
# @param size.y.use Y axis title font size
# @param size.title.use Main title font size
# @param adjust.use Adjust parameter for geom_violin
# @param point.size.use Point size for geom_violin
# @param cols.use Colors to use for plotting
# @param gene.names
# @param y.log plot Y axis on log scale
# @param x.lab.rot Rotate x-axis labels
# @param y.lab.rot Rotate y-axis labels
# @param legend.position Position the legend for the plot
# @param remove.legend Remove the legend from the plot
#
# @return A ggplot-based violin plot
#
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#
SingleVlnPlot <- function(
  feature,
  data,
  cell.ident,
  do.sort,
  y.max,
  size.x.use,
  size.y.use,
  size.title.use,
  adjust.use,
  point.size.use,
  cols.use,
  gene.names,
  y.log,
  x.lab.rot,
  y.lab.rot,
  legend.position,
  remove.legend,
  do.jitter=FALSE,
  do.boxplot = TRUE,
  print.frac=TRUE,
  is.expr=0
) {
  feature.name <- colnames(data)
  colnames(data) <- "feature"
  feature <- "feature"
  set.seed(seed = 42)
  data$ident <- cell.ident
  if (print.frac){
    perc.pos = c()
    for (i in sort(unique(cell.ident))){
      perc.pos = c(perc.pos, 100*sum(data[data$ident==i, "feature"] > is.expr) / sum(data$ident == i) )
    }
    #print(perc.pos)
    
  }
  
  if (do.sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(x = tapply(
        X = data[, feature],
        INDEX = data$ident,
        FUN = mean
      ))))
    )
  }
  if (y.log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  y.max <- set.ifnull(y.max, max(data[, feature]))
  plot <- ggplot(
    data = data,
    mapping = aes(
      x = factor(x = ident),
      y = feature
    )
  ) +
    geom_violin(
      scale = "width",
      adjust = adjust.use,
      trim = TRUE,
      mapping = aes(fill = factor(x = ident))
    ) +
    theme(
      legend.position = legend.position,
      axis.title.x = element_text(
        face = "bold",
        colour = "#990000",
        size = size.x.use
      ),
      axis.title.y = element_text(
        face = "bold",
        colour = "#990000",
        size = size.y.use
      )
    ) +
    guides(fill = guide_legend(title = NULL)) +
    xlab("Identity") +
    nogrid +
    ggtitle(feature) +
    theme(plot.title = element_text(size = size.title.use, face = "bold"))
  if (do.boxplot){
    plot <- plot + geom_boxplot(mapping = aes(fill = factor(x = ident)), width=0.1)
  }
  if (do.jitter){
    plot <- plot + geom_jitter(height = 0, size = point.size.use) 
  }
  
  plot <- plot + ggtitle(feature.name)
  if (y.log) {
    plot <- plot + scale_y_log10()
  } else {
    plot <- plot + ylim(min(data[, feature]), y.max*1.2)
  }
  if (feature %in% gene.names) {
    if (y.log) {
      plot <- plot + ylab(label = "Log Expression level")
    } else {
      plot <- plot + ylab(label = "Expression level")
    }
  } else {
    plot <- plot + ylab(label = "")
  }
  if (! is.null(x = cols.use)) {
    plot <- plot + scale_fill_manual(values = cols.use)
  }
  if (x.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  if (y.lab.rot) {
    plot <- plot + theme(axis.text.x = element_text(angle = 90))
  }
  if (remove.legend) {
    plot <- plot + theme(legend.position = "none")
  }
  if (print.frac){
    plot <- plot + annotate("text",x=1:length(unique(cell.ident)),y=max(data$feature)*1.05, label=as.character(round(perc.pos,0)))
  }
  return(plot)
}

setGeneric("vlnPlot", function(object,features.plot,cells.use=NULL,nCol=NULL,ylab.max=12,ymax.use=NULL,ymin.use=0,do.ret=FALSE,do.sort=FALSE,
                               size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,use.raw=FALSE, use.count=FALSE, do.log=FALSE,
                               ylab.use="Expression level (log TPM)",xlab.use="",print.frac=TRUE,...)  standardGeneric("vlnPlot"))
setMethod("vlnPlot","scR",
          function(object,features.plot,cells.use=NULL,nCol=NULL,ylab.max=12,ymax.use=NULL,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,use.raw=FALSE,use.count=FALSE,do.log=FALSE,ylab.use="Expression level (log TPM)",xlab.use="",print.frac=TRUE,...) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            cells.use = set.ifnull(cells.use, object@cell.names)
            data.use = data.frame(t(fetch.data(object,features.plot,cells.use = cells.use,use.imputed=use.imputed, use.raw=use.raw, use.count=use.count)));
            if (use.count) data.use = as.data.frame(log(data.use+1))
            #data.use = data.use[,colnames(object@data)]
            data.use[setdiff(features.plot, rownames(data.use)),] = 0
            data.use = data.use[features.plot,]
            #print(head(data.use))
            ident.use=drop.levels(object@ident[cells.use])
            
            pList=lapply(features.plot,function(x) plot.Vln(x,data.use[x,],ident.use,ylab.max,ymax.use,ymin.use,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use, do.log,ylab.use,xlab.use,print.frac))
            
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)  

plot.Vln=function(gene,data,cell.ident,ylab.max=12,ymax.use=NULL,ymin.use=0,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,adjust.use=1,size.use=1,cols.use=NULL, 
                  do.log=FALSE, ylab.use="Expression level (log TPM)", xlab.use="", print.frac=TRUE) {
  data$gene=as.character(rownames(data))
  data.use=data.frame(data[gene,])
  cell.ident = factor(cell.ident, levels=sort(unique(cell.ident)))
  if (length(gene)==1) {
    data.melt=data.frame(rep(gene,length(cell.ident))); colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data[1,1:length(cell.ident)])
    data.melt$id=names(data)[1:length(cell.ident)]
    
    if (print.frac){
      perc.pos = c()
      for (i in sort(unique(cell.ident))){
        perc.pos = c(perc.pos, 100*sum(data[gene, cell.ident==i] > 0) / length(data[gene, cell.ident==i]) )
      }
      
      #print(perc.pos)
      
    }
    
    
  }
  
  
  if (length(gene)>1) data.melt=melt(data.use,id="gene")
  data.melt$ident=cell.ident
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  if(do.sort) {
    data.melt$ident=factor(data.melt$ident,levels=names(rev(sort(tapply(data.melt$value,data.melt$ident,mean)))))
  }
  p=ggplot(data.melt,aes(factor(ident),value))
  p2=p + geom_violin(scale="width",adjust=adjust.use,trim=TRUE,aes(fill=factor(ident))) + ylab(ylab.use) + xlab(xlab.use)
  if (!is.null(cols.use)) {
    p2=p2+scale_fill_manual(values=cols.use)
  }
  
  # if (do.log){
  #  p2 = p2 + scale_y_log10()
  #}
  
  if (!is.null(ymax.use)){
    p2=p2+ylim(ymin.use, ymax.use)
  } else {
    #p2=p2+ylim(0, max(data.melt$value)*1.5)
  }
  
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0,size=size.use)
  p4=p3+theme_bw()+nogrid
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=90, vjust=0.5, size=15),
               axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=45, vjust=0.5, size=15)))
  p5= p5 +ggtitle(gene)+theme(plot.title = element_text(size=size.title.use, face="bold", hjust=0.5))

  if (print.frac){
    p5 = p5 + annotate("text",x=1:length(unique(cell.ident)),y=-0.5, label=as.character(round(perc.pos,0)))
  }
  return(p5)
  if(do.ret==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}



setGeneric("vlnPlot2", function(object,features.plot,nCol=NULL,ylab.max=12,ymax.use=NULL,ymin.use=0,do.ret=TRUE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,
                                do.mean=FALSE,rev.factor=FALSE,do.jitter=FALSE,use.raw=FALSE,vln.vert=FALSE,ident.use=NULL,...)  standardGeneric("vlnPlot2"))
setMethod("vlnPlot2","scR",
          function(object,features.plot,nCol=NULL,ylab.max=12,ymax.use=NULL,ymin.use=0,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,do.mean=FALSE,rev.factor=FALSE,do.jitter=FALSE,use.raw=FALSE,vln.vert=FALSE,ident.use=NULL,...) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            idents=set.ifnull(ident.use,levels(object@ident))
            cells.use=which.cells(object, idents)
            data.use=data.frame(t(fetch.data(object,features.plot,use.imputed=use.imputed, use.raw=use.raw, cells.use=cells.use)))
            data.use[setdiff(features.plot, rownames(data.use)),] = 0
            data.use = data.use[features.plot,]
            #print(head(data.use))
            ident.use=object@ident[cells.use]
            if (rev.factor) ident.use = factor(ident.use, rev(levels(ident.use)))
            pList=lapply(features.plot,function(x) plot.Vln2(x,data.use[x,],ident.use,ylab.max,ymax.use,ymin.use,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use, do.mean=do.mean, do.jitter=do.jitter, vln.vert=vln.vert))
            
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)

plot.Vln2=function(gene,data,cell.ident,ylab.max=12,ymax.use=NULL,ymin.use=0,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,
                   adjust.use=1,size.use=1,cols.use=NULL, do.mean=FALSE, do.jitter=FALSE, vln.vert = FALSE) {
  data$gene=as.character(rownames(data))
  data.use=data.frame(data[gene,])
  if (length(gene)==1) {
    data.melt=data.frame(rep(gene,length(cell.ident))); colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data[1,1:length(cell.ident)])
    data.melt$id=names(data)[1:length(cell.ident)]
  }
  #print(head(data.melt))
  
  if (length(gene)>1) data.melt=melt(data.use,id="gene")
  data.melt$ident=cell.ident
  
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  if(do.sort) {
    data.melt$ident=factor(data.melt$ident,levels=names(rev(sort(tapply(data.melt$value,data.melt$ident,mean)))))
  }
  p=ggplot(data.melt,aes(factor(ident),value))
  p2=p + geom_violin(scale="width",adjust=adjust.use,trim=TRUE,aes(fill=factor(ident))) + ylab("Expression level (log TPM)")
  if (do.mean) p2 = p2 + stat_summary(fun.y="mean", geom="point", size=2)
  if (!is.null(cols.use)) {
    p2=p2+scale_fill_manual(values=cols.use)
  }
  
  if (!is.null(ymax.use)){
    p2=p2+ylim(ymin.use, ymax.use)
  } else {
    #p2=p2+ylim(0, max(data.melt$value)*1.5)
  }
  
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+xlab("Cell Type")
  if (do.jitter) p3 = p3 + geom_jitter(height=0,size=size.use)
  p4=p3+ theme(axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=90, vjust=0.5, size=12))+theme_bw()+nogrid
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=90, vjust=0.5, size=12))+ggtitle(gene)+theme(plot.title = element_text(size=size.title.use, face="bold", hjust=0.5)))
  
  if (vln.vert) p5 = p5 + ylim(0, max(data.melt$value+0.2))
  
  if(do.ret==TRUE) {
    return(p5)
  }
  else {
    print(p5)
  }
}

setGeneric("vlnPlot.by.ident", function(object,features.plot, ident.use=NULL,id.use=NULL, nCol=NULL,ylab.max=12,ymax.use=NULL,ymin.use=NULL,do.ret=TRUE,do.sort=FALSE,
                                        size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,cols.use=NULL,use.raw=FALSE,do.log=FALSE,
                                        ylab.use="Expression level (log TPM)",xlab.use="Genes", print.frac=TRUE,...)  standardGeneric("vlnPlot.by.ident"))
setMethod("vlnPlot.by.ident","scR",
          function(object,features.plot,ident.use,id.use=NULL, nCol=NULL,ylab.max=12,ymax.use=NULL,ymin.use=NULL,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,use.raw=FALSE,do.log=FALSE,ylab.use="Expression level (log TPM)",xlab.use="Genes",print.frac=TRUE,...) {
            if (is.null(nCol)) {
              nCol=2
              if (length(features.plot)>6) nCol=3
              if (length(features.plot)>9) nCol=4
            }
            pList = list()
            
            id.use = set.ifnull(id.use, "orig")
            ident.use = set.ifnull(ident.use, levels(object@data.info[,id.use]))
            l=1;
            for (ident in ident.use){
              
              cells.use = which.cells(object, value=ident, id=id.use)     
              data.use = data.frame(t(fetch.data(object,vars.all=features.plot,cells.use=cells.use, use.imputed=use.imputed, use.raw=use.raw)));
              #Add rows that are not present
              data.use[setdiff(features.plot, rownames(data.use)),] = 0
              a = rownames(data.use)
              
              if (grepl("[|]", a[1])){
                rownames(data.use) = unlist(strsplit(a,"[|]"))[seq(1,2*length(a),2)]
                features.plot1 = unlist(strsplit(features.plot,"[|]"))[seq(1,2*length(a),2)]
              } else {
                features.plot1 = features.plot
              }
              
              df = as.data.frame(matrix(nrow=1, ncol=0))
              rownames(df) = ident
              df.ident = c()
              print(features.plot1)
              for (j in features.plot1){
                df = cbind(df, data.use[j, ])
                df.ident = c(df.ident, rep(j, length(data.use[j,])))
              }
              df.ident = factor(df.ident, levels=features.plot)
              pList[[l]] = plot.Vln(ident,df,df.ident,ylab.max,ymax.use,ymin.use,TRUE,do.sort,size.x.use,size.y.use,size.title.use,adjust.use,size.use,cols.use, do.log,ylab.use,xlab.use,print.frac)
              l=l+1;
            }
            
            
            
            if(do.ret) {
              return(pList)
            }
            else {
              multiplotList(pList,cols = nCol)
              rp()
            }
          }
)  

setGeneric("vlnPlot.vertical", function(object,features.plot,ylab.max=12,ymax.use=NULL,do.ret=FALSE,do.sort=FALSE,
                                        size.x.use=16,size.y.use=16,size.title.use=20, use.imputed=FALSE,adjust.use=1,size.use=1,
                                        cols.use=NULL,do.jitter=FALSE,use.raw=FALSE,ident.include=NULL,...)  standardGeneric("vlnPlot.vertical"))
setMethod("vlnPlot.vertical","scR",
          function(object,features.plot,ylab.max=12,ymax.use=NULL,do.ret=FALSE,do.sort=FALSE,size.x.use=16,size.y.use=16,size.title.use=20,use.imputed=FALSE,adjust.use=1,
                   size.use=1,cols.use=NULL,do.jitter=FALSE,use.raw=FALSE,ident.include=NULL,...) {
            
            ident.include = set.ifnull(ident.include, levels(object@ident))
            ident.include = rev(ident.include)
            plots = VlnPlot(object, features.plot, return.plotlist = TRUE, do.return=TRUE, do.boxplot = FALSE, print.frac = FALSE, ident.include=ident.include)
            plots[[1]]=plots[[1]]+coord_flip()+ scale_y_continuous(expand = c(0,0.1))+ labs(x=NULL, y=NULL)+
              theme(line = element_blank(),axis.text.x = element_blank(), axis.text.y = element_text(angle=0),
                    legend.position="none",plot.margin = unit(c(0,0,0,0),"cm"), panel.border=element_rect(size=2), panel.margin=unit(c(0,0,0,0),"cm"))
            for(i in 2:length(features.plot)){
              plots[[i]]=plots[[i]]+coord_flip()+ scale_y_continuous(expand = c(0,0.1))+ labs(x=NULL, y=NULL)+
                theme(line = element_blank(),axis.text = element_blank(),legend.position="none",plot.margin = unit(c(0,0,0,0),"cm"),
                      panel.border=element_rect(size=2), panel.spacing=unit(c(0,0,0,0),"cm"))
            }
            do.call("grid.arrange",c(plots, ncol=length(plots)))
            
          }
)  


############################### SCATTER PLOTS ####################################################

setGeneric("genePlot", function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                                pch.use=16,cex.use=2,use.imputed=FALSE,do.ident=FALSE,...)  standardGeneric("genePlot"))
setMethod("genePlot","scR",
          function(object, gene1, gene2, cell.ids=NULL,col.use=NULL,
                   pch.use=16,cex.use=2,use.imputed=FALSE,do.ident=FALSE,...) {
            cell.ids=set.ifnull(cell.ids,object@cell.names)
            data.use=t(fetch.data(object,c(gene1,gene2),cells.use = cell.ids,use.imputed=use.imputed))
            corner(data.use)
            g1=as.numeric(data.use[gene1,cell.ids])
            g2=as.numeric(data.use[gene2,cell.ids])
            ident.use=as.factor(object@ident[cell.ids])
            if (length(col.use)>1) {
              col.use=col.use[as.numeric(ident.use)]
            }
            else {
              col.use=set.ifnull(col.use,as.numeric(ident.use))
            }
            gene.cor=round(cor(g1,g2),2)
            plot(g1,g2,xlab=gene1,ylab=gene2,col=col.use,cex=cex.use,main=gene.cor,pch=pch.use,...)
            if (do.ident) {
              return(identify(g1,g2,labels = cell.ids))
            }
          }
)


setGeneric("cellPlot", function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE, plot.return=FALSE, xlab.use=NULL, ylab.use=NULL,use.count=FALSE,...)  standardGeneric("cellPlot"))
setMethod("cellPlot","scR",
          function(object, cell1, cell2, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,plot.return=FALSE, xlab.use=NULL, ylab.use=NULL,use.count=FALSE,...) {
            gene.ids=set.ifnull(gene.ids,rownames(object@data))
            if (length(cell1)==1){
              c1=as.numeric(object@data[gene.ids,cell1])
              if (use.count) c1 = as.numeric(object@count.data[gene.ids,cell1])
              xlab.use = set.ifnull(xlab.use, cell1)
            } else {
              if (use.count){
                c1 = apply(object@count.data[gene.ids, cell1], 1, function(x) mean(x))
              } else {
                c1 = apply(object@data[gene.ids, cell1], 1, function(x) expMean(x))
              }
              xlab.use = set.ifnull(xlab.use, unique(object@ident[cell1]))
            }
            if (length(cell2)==1){
              c2=as.numeric(object@data[gene.ids,cell2])
              if (use.count) c2 = as.numeric(object@count.data[gene.ids,cell2])
              ylab.use = set.ifnull(ylab.use, cell2)
            } else {
              if (use.count) {
                c2 = apply(object@count.data[gene.ids, cell2], 1, function(x) mean(x))
              } else{
                c2 = apply(object@data[gene.ids, cell2], 1, function(x) expMean(x))
              }
              ylab.use = set.ifnull(ylab.use, unique(object@ident[cell2]))
            }
            gene.cor=round(cor(c1,c2),2)
            
            if (!plot.return){
              smoothScatter(c1,c2,xlab=xlab.use,ylab=ylab.use,col=col.use,nrpoints=nrpoints.use,pch=pch.use,cex=cex.use,main=paste0("r = ",gene.cor))
              if (do.ident) {
                identify(c1,c2,labels = gene.ids)
              }
            } else {
              smoothScatter(c1,c2,xlab=xlab.use,ylab=ylab.use,col=col.use,nrpoints=nrpoints.use,pch=pch.use,cex=cex.use,main=paste0("r = ",gene.cor))
              p = recordPlot()
              return(p)
            }
            
          }
          
)

#Remove
setGeneric("cellPlot_AvgVsPop", function(object, pop, sing.cells, xlab.use=NULL, ylab.use=NULL, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...)  standardGeneric("cellPlot_AvgVsPop"))
setMethod("cellPlot_AvgVsPop","scR",
          function(object,pop, sing.cells, xlab.use=NULL, ylab.use=NULL, gene.ids=NULL,col.use="black",nrpoints.use=Inf,pch.use=16,cex.use=0.5,do.ident=FALSE,...) {
            gene.ids=set.ifnull(gene.ids,rownames(object@data))
            c2=apply(object@data[gene.ids,sing.cells],1, expMean)
            xlab.use=set.ifnull(xlab.use, paste("[", pop, "]"))
            ylab.use=set.ifnull(ylab.use, paste("Avg. [", sing.cells, "]"))
            if (length(pop) == 1){
              c1=as.numeric(object@data[gene.ids,pop])
            } else {
              c1 = apply(object@data[gene.ids, pop],1, expMean)
            }
            gene.cor=round(cor(c1,c2),2)
            smoothScatter(c1,c2,xlab=xlab.use,ylab=ylab.use,col=col.use,nrpoints=nrpoints.use,pch=pch.use,cex=cex.use,main=gene.cor)
            if (do.ident) {
              identify(c1,c2,labels = gene.ids)
            }
          }
)

#Remove
setGeneric("pairwiseCellCorrPlot", function(object, cell.names=NULL, gene.ids=NULL,col.use="black",point.size=1,...)  standardGeneric("pairwiseCellCorrPlot"))
setMethod("pairwiseCellCorrPlot","scR",
          function(object, cell.names=NULL, gene.ids=NULL,col.use="black",point.size=1,...) {
            gene.ids=set.ifnull(gene.ids,rownames(object@data))
            cell.names = set.ifnull(cell.names, colnames(object@data)[1:5])
            data.expression = object@data[gene.ids, cell.names]
            
            corrplot <- ggpairs(data.expression) 
            for (i in 1:length(cell.names)){
              for (j in (i+1):length(cell.names)){
                
                p <- ggplot(data.expression, aes_string(x=cell.names[i],y=cell.names[j])) + geom_point(size=point.size)
                corrplot <- putPlot(corrplot,p,j,i)
                
              }
            }
            print(corrplot)
            
          }
          
)


##########SATURATION PLOTS TO DIAGNOSE SEQUENCING DEPTH###########################

setGeneric("doSaturationPlots", function(object,ident.use=NULL,id =NULL,UMI.data=TRUE,cells.use=NULL,size.x.use=16,size.y.use=16,size.title.use=20,downsamp.levels=seq(0.1,1,length=10),high.thres=10,...) standardGeneric("doSaturationPlots"))
setMethod("doSaturationPlots", "scR", 
          function(object,ident.use=NULL, id=NULL, UMI.data=TRUE,cells.use=NULL,size.x.use=16,size.y.use=16,size.title.use=20,downsamp.levels=seq(0.1,1,length=10),high.thres=10,...) {
            
            ident.use = set.ifnull(ident.use, levels(object@ident)[1])
            cells.use = set.ifnull(cells.use, names(object@ident[which(object@ident %in% ident.use)]))
            
            pList=list()
            
            #Reads vs. # Genes
            p = geneSaturationPlot(object, ident.use=ident.use, id=id, cells.use=cells.use,size.x.use=size.x.use,
                                   size.y.use=size.y.use,size.title.use=size.title.use,downsamp.levels=downsamp.levels)
            #Reads vs. # Transcripts
            if (UMI.data){
              pList[[1]] = p
              pList[[2]] = TranscriptSaturationPlot(object, ident.use=ident.use, id=id, cells.use=cells.use,size.x.use=size.x.use,
                                                    size.y.use=size.y.use,size.title.use=size.title.use,downsamp.levels=downsamp.levels)
              #Reads vs. Transcripts
              reads = as.numeric(as.matrix(object@reads.data[,cells.use]))
              names(reads) = rep(rownames(object@reads.data), ncol(object@reads.data[,cells.use]))
              UMIs = as.numeric(as.matrix(object@count.data[,cells.use]))
              names(UMIs) = rep(rownames(object@count.data), ncol(object@count.data[,cells.use]))
              ind = (reads != 0); 
              reads = reads[ind]; UMIs = UMIs[ind]
              
              high_genes = which(reads/UMIs > high.thres)
              df = data.frame(reads=reads, UMIs=UMIs)
              samp_size = min(5e5, nrow(df))
              samples = union(sample(nrow(df), samp_size), c(order(df$reads, decreasing=TRUE)[1:100], order(df$UMIs, decreasing=TRUE)[1:100]))
              #samples = sample(nrow(df), samp_size)
              p = ggplot(df[samples,], aes(reads,UMIs)) + geom_point() + scale_x_log10() + scale_y_log10()
              p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=0, vjust=0.5, size=12))
              p <- p+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) 
              p <- p+ ggtitle(ident.use)+ theme(plot.title = element_text(size=size.title.use, face="bold"), legend.position="none")
              p <- p + xlab("# Reads") + ylab("# Transcripts (UMIs)") + annotate("text", x=reads[high_genes],y=UMIs[high_genes],
                                                                                 label=names(high_genes), size=3)
              pList[[3]] = p
              
              multiplotList(pList, cols=3)
              rp()
            } else {
              p
            }
            
          }
)


setGeneric("geneSaturationPlot", function(object,ident.use=NULL,id =NULL, cells.use=NULL,downsamp.levels=seq(0.1,1,length=10),size.x.use=16,size.y.use=16,size.title.use=20,
                                          adjust.use=1,size.use=1,...) standardGeneric("geneSaturationPlot"))
setMethod("geneSaturationPlot", "scR", 
          function(object,ident.use=NULL, cells.use=NULL,downsamp.levels=seq(0.1,1,length=10),size.x.use=16,size.y.use=16,size.title.use=20,
                   adjust.use=1,size.use=1,...) {
            
            ident.use = set.ifnull(ident.use, levels(object@ident)[1])
            cells.use = set.ifnull(cells.use, names(object@ident[which(object@ident %in% ident.use)]))
            data.use = object@reads.data[,cells.use]
            df = data.frame(matrix(nrow=length(cells.use), ncol=0))
            rownames(df) = cells.use
            for (i in downsamp.levels){
              print(i)
              num.genes = nGene_downsample(data.use, downsample.level = i,nboot=50)
              df = cbind(df, data.frame(num.genes))
            }
            
            df = t(df)
            rownames(df) = as.character(downsamp.levels)
            df = cbind(df, data.frame(downsamp.frac=downsamp.levels))
            
            df <- melt(df ,  id = 'downsamp.frac', variable_name = 'cell')
            p <- ggplot(df, aes(downsamp.frac, value)) + geom_line(aes(color = cell)) + theme_bw() + nogrid 
            p = p + ylab("# genes detected")
            p = p + xlab("Downsampling fraction (Reads)")
            p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=0, vjust=0.5, size=12))
            p <- p+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) 
            p <- p+ ggtitle(ident.use)+ theme(plot.title = element_text(size=size.title.use, face="bold"), legend.position="none") + xlim(c(min(downsamp.levels), max(downsamp.levels)))
            
            return(p)
          }
)

setGeneric("TranscriptSaturationPlot", function(object,ident.use=NULL,id =NULL, cells.use=NULL,downsamp.levels=seq(0.1,1,length=10),size.x.use=16,size.y.use=16,size.title.use=20,
                                                adjust.use=1,size.use=1,...) standardGeneric("TranscriptSaturationPlot"))
setMethod("TranscriptSaturationPlot", "scR", 
          function(object,ident.use=NULL, cells.use=NULL,downsamp.levels=seq(0.1,1,length=10),size.x.use=16,size.y.use=16,size.title.use=20,
                   adjust.use=1,size.use=1,...) {
            
            trans.data = object@count.data[,cells.use]
            reads.data = object@reads.data[,cells.use]
            df = data.frame(matrix(nrow=length(cells.use), ncol=0))
            rownames(df) = cells.use
            for (i in downsamp.levels){
              num.trans = nTrans_downsample(reads.data, trans.data, downsample.level = i,nboot=50)
              df = cbind(df, data.frame(num.trans))
            }
            
            df = t(df)
            rownames(df) = as.character(downsamp.levels)
            df = cbind(df, data.frame(downsamp.frac=downsamp.levels))
            
            df <- melt(df ,  id = 'downsamp.frac', variable_name = 'cell')
            p <- ggplot(df, aes(downsamp.frac, value)) + geom_line(aes(color = cell)) + theme_bw() + nogrid 
            
            p = p + ylab("# transcripts detected")
            p = p + xlab("Downsampling fraction (Reads)")
            p <- p + theme(axis.title.x = element_text(face="bold", colour="#990000", size=size.x.use), axis.text.x  = element_text(angle=0, vjust=0.5, size=12))
            p <- p+theme(axis.title.y = element_text(face="bold", colour="#990000", size=size.y.use), axis.text.y  = element_text(angle=0, vjust=0.5, size=12)) 
            p <- p+ ggtitle(ident.use)+ theme(plot.title = element_text(size=size.title.use, face="bold"), legend.position="none")
            
            return(p)
          }
)


# setGeneric("Perc.pos.by.ident", function(object,features.use=NULL, ident.use=NULL, thresh.use=NULL,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,return.plot=FALSE,...) standardGeneric("Perc.pos.by.ident"))
# setMethod("Perc.pos.by.ident", "scR", 
#           function(object,features.use=NULL,ident.use=NULL, thresh.use=NULL,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,return.plot=FALSE,...) {
#             
#             features.use=set.ifnull(features.use, object@var.genes[1:20])
#             features.use=features.use[features.use %in% rownames(object@data)]
#             ident.use=set.ifnull(ident.use, levels(object@ident))
#             thresh.use=set.ifnull(thresh.use,object@is.expr)
#             #Matrix of percent expressing cells
#             PercMat = matrix(0, nrow=length(features.use), ncol = 0)
#             rownames(PercMat) = features.use; 
#             
#             #Matrix of average transcript levels
#             ExpMat = PercMat;
#             #Count mat
#             Count.mat = object@count.data[features.use, colnames(object@data)]
#             
#             
#             for (i in ident.use){
#               cells.in.cluster = names(object@ident)[which(object@ident== i)]
#               vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
#               PercMat = cbind(PercMat,vec.exp)
#               
#               vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
#               ExpMat = cbind(ExpMat, vec.exp)
#             }
#             colnames(ExpMat) = ident.use
#             colnames(PercMat) = ident.use
#             
#             rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
#             PercMat = PercMat[rows.use,]
#             ExpMat = ExpMat[rows.use,]
#             features.use = rows.use
#             if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
#             if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
#             
#             if (!is.null(norm.exp)){
#               if (norm.exp < 0 | norm.exp > 1){
#                 print("Warning: norm.exp should be a value between (0,1). Skipping normalization")
#                 next
#               } else{
#                 quant.vals = apply(ExpMat,1, function(x) quantile(x, norm.exp))
#                 ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
#               }
#             }
#             
#             
#             
#             
#             if (do.plot){
#               
#               ExpVal = melt(ExpMat)
#               PercVal = melt(PercMat)
#               colnames(ExpVal) = c("gene","cluster","nTrans")
#               ExpVal$percExp = PercVal$value*100
#               
#               if (!do.transpose){
#                 ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
#                 ExpVal$cluster = factor(ExpVal$cluster, levels= ident.use)
#                 p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
#                   scale_color_gradient(low ="blue",   high = "red", limits=c(min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
#                 p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
#                   theme(axis.text.y=element_text(size=12, face="italic"))  
#                 if (return.plot){
#                   return(p)
#                 } else{
#                   print(p)
#                 }
#                 
#               } else {
#                 ExpVal$gene = factor(ExpVal$gene, levels=features.use)
#                 ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
#                 p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
#                   scale_color_gradient(low ="blue",   high = "red", limits=c( min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
#                 p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
#                   theme(axis.text.y=element_text(size=12, face="italic"))
#                 if (return.plot){
#                   return(p)
#                 } else{
#                   print(p)
#                 }
#                 
#               }
#               
#             }else {
#               to.return=list()
#               to.return$ExpMat = ExpMat;
#               to.return$PercMat = PercMat;
#               return(to.return)
#             }
#           }
# )

Perc.pos.by.ident = function(object,features.use=NULL,features.sign = NULL, Count.mat = NULL, use.raw=FALSE, ident.use=NULL, thresh.use=0,do.plot=TRUE,do.transpose=FALSE,
                             max.val.perc=NULL, max.val.exp=NULL,col.low="white", col.high="red",max.size=10,min.perc=0,norm.exp=NULL,
                             x.lab.rot=FALSE,title.use=NULL,name.replace.vec=NULL,make.clust.diag=FALSE,return.ident.order=FALSE,
                             place.sham.value=FALSE,sham.val.perc=75,
                             return.plot=FALSE,alpha.use=1,use.median=FALSE, convertNaNs = FALSE,...) {
  
  features.use=set.ifnull(features.use, object@var.genes[1:20])
  if (use.raw){
    data.use = object@count.data[,object@cell.names]
    Count.mat = object@count.data[,object@cell.names]
  } else {
    data.use = object@data
  }
  
  if (is.character(features.use)){
    features.use=features.use[features.use %in% rownames(data.use)]
    all_features = features.use
    feature.rownames = all_features
  } else if(is.list(features.use)){
    if (is.null(features.sign)){
      features.sign = list()
      for (item in 1:length(features.use)){
        features.sign[[item]] = rep("+",length(features.use[[item]]))
      }
    } else {
      if (length(features.sign) != length(features.use)){
        stop("Error: features.sign must have the same length as features.use")
      }
        for (item in 1:length(features.use)){
          if (features.sign[[item]] == ""){
            features.sign[[item]] = rep("+",length(features.use[[item]]))
          } else {
            if (length(features.sign[[item]]) != length(features.use[[item]]) ){
              stop(paste0("Error: entry ", item, " in features.sign is not of the correct length as the corresponding item in features.use"))
            }
          }
        }
      }
    all_features = c(); feature.rownames = c();
    for (l in 1:length(features.use)){
      #if (length(features.use[[l]]) > 2){
      #  stop("Error: Must specify gene sets of length 2")
      #}
      features.sign[[l]] = features.sign[[l]][features.use[[l]] %in% rownames(data.use)]
      features.use[[l]]=features.use[[l]][features.use[[l]] %in% rownames(data.use)]
      
      feature.rownames = c(feature.rownames, paste0(features.use[[l]], features.sign[[l]], collapse=""))
      all_features = union(all_features, features.use[[l]])
    }
  }

  ident.use=set.ifnull(ident.use, levels(object@ident))
  
  #Matrix of percent expressing cells
  PercMat = matrix(0, nrow=length(features.use), ncol = 0)
  rownames(PercMat) = feature.rownames; 
  
  #Matrix of average transcript levels
  ExpMat = PercMat;
  
  #Count mat
  Count.mat = set.ifnull(Count.mat, exp(data.use[all_features, colnames(data.use)]) - 1)
  
  
  for (i in ident.use){
    cells.in.cluster = names(object@ident)[which(object@ident== i)]
    if (is.character(features.use)){
      vec.perc = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
      if (!use.median){
        vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {0})
      } else {
        vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ median(x[x>0]) } else {0})
      }
    } else {
      vec.perc = c(); vec.exp=c()
      for (l in 1:length(features.use)){
        genes = features.use[[l]]
        if (length(genes) > 1){
          x1 = as.matrix(Count.mat[genes, cells.in.cluster])
          tempMat = matrix(0, nrow=nrow(x1), ncol = ncol(x1))
          for (rowi in 1:nrow(x1)){
            if (features.sign[[l]][rowi] == "+"){
              tempMat[rowi,] = as.numeric(x1[rowi,] > thresh.use)
            }
            
            if (features.sign[[l]][rowi] == "-"){
              tempMat[rowi,] = as.numeric(x1[rowi,] <= thresh.use)
            }
            
          }
          temp = apply(tempMat, 2, prod)
          vec.perc = c(vec.perc, sum(temp)/length(temp))
          cells.all = cells.in.cluster[temp==1]
          
          if (length(cells.all)<2){
            vec.exp = c(vec.exp, 0)
          } else {
            vec.exp = c(vec.exp,mean(x1[features.use[[l]][features.sign[[l]] == "+"],cells.all]))
          }
        } else {
          
          x1 = Count.mat[genes, cells.in.cluster]
          vec.perc = c(vec.perc, sum(x1 > thresh.use)/length(x1))
          vec.exp = c(vec.exp, if (sum(x1>0) > 1) { mean(x1[x1>0]) } else {0} ) 
          
        }
      }
      names(vec.perc) = feature.rownames
      names(vec.exp) = feature.rownames
    }
    
    PercMat = cbind(PercMat,vec.perc)
    ExpMat = cbind(ExpMat, vec.exp)
  }
  colnames(ExpMat) = ident.use
  colnames(PercMat) = ident.use

  
  if (!is.null(name.replace.vec)){
    rownames(ExpMat) = name.replace.vec[rownames(ExpMat)]
    rownames(PercMat) = name.replace.vec[rownames(PercMat)]
  }
  
  if (!is.null(norm.exp)){
    if (norm.exp < 0 | norm.exp > 1){
      print("Warning: norm.exp should be a value between (0,1). Skipping normalization")
      next
    } else{
      quant.vals = apply(ExpMat,1, function(x) quantile(x, norm.exp))
      ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
    }
  }
  
 
  rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
  PercMat = PercMat[rows.use,]
  ExpMat = ExpMat[rows.use,]
  feature.rownames = rows.use

  if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
  if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
  
  if (make.clust.diag){
    factor.levels = c()
    for (i1 in rownames(PercMat)){
      if (min(PercMat[i1,]) > 0.35 | max(PercMat[i1,]) < 0.15) next
      ind.sort = colnames(PercMat)[order(PercMat[i1,], decreasing=TRUE)]
      ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
      factor.levels = c(factor.levels, ind.sort[1])
    }
    factor.levels = c(factor.levels, setdiff(colnames(PercMat), factor.levels))
    factor.levels = factor.levels[!is.na(factor.levels)]
    ident.use=factor.levels
  } 
  
 
  
  if (return.ident.order){
    return(ident.use)
  }
  
  if (place.sham.value){
    ind = which(PercMat == 0)[1]
    PercMat[ind] = sham.val.perc/100
    ExpMat[ind] = 1e-10
  }
  
  if (convertNaNs){
    PercMat[is.nan(PercMat)] = 0
    ExpMat[is.nan(ExpMat)] = 0
  }
  
  
  if (do.plot){
    
    ExpVal = melt(ExpMat)
    PercVal = melt(PercMat)
    colnames(ExpVal) = c("gene","cluster","logExp")
    ExpVal$percExp = PercVal$value*100

    if (!do.transpose){
      ExpVal$gene = factor(ExpVal$gene, levels=rev(feature.rownames))
      ExpVal$cluster = factor(ExpVal$cluster, levels= ident.use)
      p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = logExp,  size =percExp), alpha=alpha.use) + 
       scale_color_gradient(low =col.low,   high = col.high, limits=c(0, max(ExpVal$logExp) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid + 
        #scale_fill_brewer(palette = "Blues")+scale_size(range = c(1, max.size))+   theme_bw() +nogrid + 
        ggtitle(title.use) + theme(plot.title = element_text(hjust = 0.5, size=15))
      if (x.lab.rot) {
        p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
          theme(axis.text.y=element_text(size=12, face="italic"))  
      } else { p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
        theme(axis.text.y=element_text(size=12, face="italic"))  
      }
    } else {
      ExpVal$gene = factor(ExpVal$gene, levels=feature.rownames)
      ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
      p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = logExp,  size =percExp), alpha=alpha.use) + 
        scale_color_gradient(low =col.low,   high = col.high, limits=c( 0, max(ExpVal$logExp) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid +
        #scale_fill_brewer(palette = "Blues")+scale_size(range = c(1, max.size))+   theme_bw() +nogrid + 
        ggtitle(title.use) + theme(plot.title = element_text(hjust = 0.5, size=15))
      if (x.lab.rot) {
        p = p + ylab("Cluster") + xlab("Gene") +  theme(axis.text.y=element_text(size=12, face="italic")) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1))
      } else {
        p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
          theme(axis.text.y=element_text(size=12, face="italic"))
      }
      
      if (return.plot){
        return(p)
      } else {
        print(p)
      }
      
      
    }
    
  }else {
    to.return=list()
    to.return$ExpMat = ExpMat;
    to.return$PercMat = PercMat;
    return(to.return)
  }
}


setGeneric("Perc.pos.by.ident.in.clust", function(object,features.use=NULL, clust.use=NULL, ident.use="orig", thresh.use=NULL,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,return.plot=FALSE,...) standardGeneric("Perc.pos.by.ident.in.clust"))
setMethod("Perc.pos.by.ident.in.clust", "scR", 
          function(object,features.use=NULL,clust.use=NULL, ident.use="orig", thresh.use=NULL,do.plot=TRUE,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,cols.use=gray.colors(10),max.size=10,min.perc=0,norm.exp=NULL,return.plot=FALSE,...) {
            
            features.use=set.ifnull(features.use, object@var.genes[1:20])
            features.use=features.use[features.use %in% rownames(object@data)]
            clust.use = set.ifnull(clust.use, 1)
            print(paste0("Examining genes across samples in ", clust.use))
            object = subsetData(object, cells.use=which.cells(object, clust.use))
            object = set.all.ident(object, ident.use)
            thresh.use=set.ifnull(thresh.use,object@is.expr)
            #Matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Matrix of average transcript levels
            ExpMat = PercMat;
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            idents.in.object = levels(object@ident)
            for (i in idents.in.object){
              cells.in.cluster = names(object@ident)[which(object@ident== i)]
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) =  idents.in.object
            colnames(PercMat) =  idents.in.object
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
            if (!is.null(norm.exp)){
              if (norm.exp < 0 | norm.exp > 1){
                print("Warning: norm.exp should be a value between (0,1). Skipping normalization")
                next
              } else{
                quant.vals = apply(ExpMat,1, function(x) quantile(x, norm.exp))
                ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
              }
            }
            
            
            
            
            if (do.plot){
              
              ExpVal = melt(ExpMat)
              PercVal = melt(PercMat)
              colnames(ExpVal) = c("gene","cluster","nTrans")
              ExpVal$percExp = PercVal$value*100
              
              if (!do.transpose){
                ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
                ExpVal$cluster = factor(ExpVal$cluster, levels= idents.in.object)
                p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                  scale_color_gradient(low ="blue",   high = "red", limits=c(min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
                p = p + xlab("Sample") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                  theme(axis.text.y=element_text(size=12, face="italic")) +  ggtitle(paste0("Expression in Cluster", clust.use)) + 
                  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold')) 
                if (return.plot){ 
                  return(p)
                } else{
                  print(p)
                }
                
              } else {
                ExpVal$gene = factor(ExpVal$gene, levels=features.use)
                ExpVal$cluster = factor(ExpVal$cluster, levels= rev(idents.in.object))
                p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                  scale_color_gradient(low ="blue",   high = "red", limits=c( min(ExpVal$nTrans), max(ExpVal$nTrans) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
                p = p + ylab("Sample") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                  theme(axis.text.y=element_text(size=12, face="italic")) + ggtitle(paste0("Expression in Cluster", clust.use)) + 
                  theme(plot.title = element_text(hjust = 0.5, size=12, face='bold')) 
                if (return.plot){
                  return(p)
                } else{
                  print(p)
                }
                
              }
              
            }else {
              to.return=list()
              to.return$ExpMat = ExpMat;
              to.return$PercMat = PercMat;
              return(to.return)
            }
          }
)


# Distribution of cluster distances by ident
setGeneric("cluster.dist.by.ident", function(object,ident1.use=NULL, ident2.use="orig", ident2.counts=NULL, ident1.names=NULL,  ident2.names=NULL, cells.use=NULL, ylab.use="%", xlab.use="Sample", legend.title="cluster", seed.use=42) standardGeneric("cluster.dist.by.ident"))
setMethod("cluster.dist.by.ident", "scR", 
          function(object,ident1.use=NULL, ident2.use="orig",ident2.counts=NULL, ident1.names=NULL, ident2.names=NULL, cells.use=NULL, ylab.use="%", xlab.use="Sample", legend.title="cluster", seed.use=42) {
            require(randomcoloR)
            cells.use = set.ifnull(cells.use, colnames(object@data))
            ident1.use = set.ifnull(ident1.use,"m")
            if (!(ident1.use %in% colnames(object@data.info)) | !(ident2.use %in% colnames(object@data.info))){
              stop("One of ident1.use or ident2.use is invalid")
            }
            
            
            ident2.counts = set.ifnull(ident2.counts, table(object@data.info[cells.use,c(ident2.use)]))
            
            if (length(ident2.counts) != length(table(object@data.info[cells.use,ident2.use])) ){
              stop("The number of entries in ident2.counts must match the number of unique labels in ident2")
            }
            
            if (is.null(names(ident2.counts))){
              names(ident2.counts) = names(table(object@data.info[cells.use,ident2.use]))
            }
            
            
            print("Using the following base counts for ident.2")
            print(ident2.counts)
            
            # Counts
            object@data.info[,ident1.use] = as.factor(object@data.info[,ident1.use])
            object@data.info[,ident2.use] = as.factor(object@data.info[,ident2.use])
            df = table(object@data.info[cells.use,c(ident1.use, ident2.use)])
            df = scale(df, center=FALSE, scale=ident2.counts)
            df=df*100# 100 to convert to percentages
            df.melt = melt(df)
            df.melt[,ident1.use] = factor(df.melt[,ident1.use])
            df.melt[,ident2.use] = factor(df.melt[,ident2.use])
            df.melt[,ident1.use] = factor(df.melt[,ident1.use], levels=levels(object@data.info[,ident1.use]))
            df.melt[,ident2.use] = factor(df.melt[,ident2.use], levels=levels(object@data.info[,ident2.use]))
            if ((ylab.use) == "%") ylab.use = paste0("Percent", capitalize(ident2.use))
            ylab.use = gsub(" ","", ylab.use)
            colnames(df.melt)[3] = ylab.use
            require(RColorBrewer)
            nColors = length(unique(df.melt[,ident1.use]))
            if (nColors <=8 ){
              col.use= brewer.pal(nColors,"Set2")
            } else {
              require(randomcoloR)
              set.seed(seed.use)
              col.use = unname(distinctColorPalette(nColors))
            }
            p = ggplot() + geom_bar(aes_string(x=ident2.use,y=ylab.use,fill=ident1.use), data=df.melt, stat="identity") + xlab(xlab.use) +
              ylab(ylab.use) + scale_fill_manual(values = col.use)  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            #rainbow(length(unique(df.melt[,ident1.use])))
            #+  scale_fill_brewer(palette="Set2", guide=guide_legend(title = legend.title))  
            print(p)
          }
)

# Distribution of cluster distances by ident
setGeneric("cluster.freq.by.ident", function(object,ident1.use=NULL, clust.use=NULL, ident2.use="orig", ident2.counts=NULL, ident1.names=NULL,  ident2.names=NULL, cells.use=NULL, ylab.use="%", xlab.use="Sample", legend.title="cluster") standardGeneric("cluster.freq.by.ident"))
setMethod("cluster.freq.by.ident", "scR", 
          function(object,ident1.use=NULL,clust.use=NULL, ident2.use="orig",ident2.counts=NULL, ident1.names=NULL, ident2.names=NULL, cells.use=NULL, ylab.use="%", xlab.use="Sample", legend.title="cluster") {
            
            cells.use = set.ifnull(cells.use, colnames(object@data))
            ident1.use = set.ifnull(ident1.use,"m")
            if (!(ident1.use %in% colnames(object@data.info)) | !(ident2.use %in% colnames(object@data.info))){
              stop("One of ident1.use or ident2.use is invalid")
            }
            
            if (!(clust.use %in% unique(object@data.info[,ident1.use]) )) stop("Error: Must specify valid cluster")
            
            
            ident2.counts = set.ifnull(ident2.counts, table(object@data.info[cells.use,c(ident2.use)]))
            
            if (length(ident2.counts) != length(table(object@data.info[cells.use,ident2.use])) ){
              stop("The number of entries in ident2.counts must match the number of unique labels in ident2")
            }
            
            if (is.null(names(ident2.counts))){
              names(ident2.counts) = names(table(object@data.info[cells.use,ident2.use]))
            }
            
            
            print("Using the following base counts for ident.2")
            print(ident2.counts)
            
            # Counts
            object@data.info[,ident1.use] = as.factor(object@data.info[,ident1.use])
            object@data.info[,ident2.use] = as.factor(object@data.info[,ident2.use])
            df = table(object@data.info[cells.use,c(ident1.use, ident2.use)])
            df = as.matrix(scale(df, center=FALSE, scale=ident2.counts))
            #df = df / sum(df)
            
            df=df*100# 100 to convert to percentages
            df = df[clust.use,]
            df = data.frame(df)
            colnames(df) = paste0("Percentage_Clust",clust.use)
            ylab.use = colnames(df)
            df$ident = factor(rownames(df), levels=rownames(df))
            p = ggplot() + geom_bar(aes_string(y=ylab.use, x="ident"), data=df, stat="identity", fill="green",col="black") +
              ylab(ylab.use)
            #rainbow(length(unique(df.melt[,ident1.use])))
            #+  scale_fill_brewer(palette="Set2", guide=guide_legend(title = legend.title))  
            print(p)
          }
)


setGeneric("plot_dendro_withdotplot", function(object, genes.use = NULL, use.raw=FALSE, max.val.exp=5, max.val.perc=1,top.mar=1, left.mar = 1.8,right.mar=1.6, bottom.mar = -1,layout.use=NULL,return.plot=TRUE,...) standardGeneric("plot_dendro_withdotplot"))
setMethod("plot_dendro_withdotplot","scR",
          function(object, genes.use=NULL, use.raw=FALSE, max.val.exp=5, max.val.perc=1,top.mar=1, left.mar = 4.5,right.mar=1.5, bottom.mar = -1,layout.use=NULL,return.plot=TRUE,...) {
            require(ggdendro)
            # Get dendrogram
            tree = object@cluster.tree[[1]]
            data = as.hclust(tree)
            data.tree=as.dendrogram(data)
            ddata <- dendro_data(data.tree, type = "rectangle")
            pdendro <- ggplot(segment(ddata)) + 
              geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
              theme_dendro() + theme(plot.margin = unit(c(top.mar,right.mar,bottom.mar,left.mar), "cm")) # (top, right, bottom, left )
            
            # get node labels according to order
            ident.order=tree$tip.label
            nodes.1=ident.order[getLeftDecendants(tree,tree$Nnode+2)]
            nodes.2=ident.order[getRightDecendants(tree,tree$Nnode+2)]
            ident.order = c(nodes.1,nodes.2)
            
            if (use.raw){
              data.use = object@count.data
            } else {
              data.use = object@data
            }
            
            if (is.null(genes.use)){
              # find 
              genes.use = sample(rownames(data.use),5)
            }
            genes.use1 = genes.use[genes.use %in% rownames(data.use)]
            if (length(genes.use1) < length(genes.use)){
              print(paste0(length(setdiff(genes.use, genes.use1)), " genes not in data. :")) 
              print(setdiff(genes.use, genes.use1))
            }
            genes.use=genes.use1
            print(length(genes.use))
            p <- Perc.pos.by.ident(object,features.use = genes.use,ident.use = ident.order, max.val.exp = max.val.exp, max.val.perc = max.val.perc,return.plot = TRUE, use.raw=use.raw,...)
            print(length(genes.use))
            multiplot(pdendro,p, cols=1, layout=layout.use)
            if (return.plot) return(p)
          }
)


# setGeneric("plot_dendro_withgenesinsampledotplot", function(object, gene.use = NULL, max.val.exp=5, max.val.perc=1,top.mar=1, left.mar = 1.8,right.mar=1.6, bottom.mar = -1,...) standardGeneric("plot_dendro_withgeneinsampledotplot"))
# setMethod("plot_dendro_withgeneinsampledotplot","scR",
#           function(object, gene.use=NULL,max.val.exp=5, max.val.perc=1,top.mar=1, left.mar = 4.5,right.mar=1.5, bottom.mar = -1,...) {
#             require(ggdendro)
#             # Get dendrogram
#             tree = object@cluster.tree[[1]]
#             data = as.hclust(tree)
#             data.tree=as.dendrogram(data)
#             ddata <- dendro_data(data.tree, type = "rectangle")
#             pdendro <- ggplot(segment(ddata)) + 
#               geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#               theme_dendro() + theme(plot.margin = unit(c(top.mar,right.mar,bottom.mar,left.mar), "cm")) # (top, right, bottom, left )
#             
#             # get node labels according to order
#             ident.order=tree$tip.label
#             nodes.1=ident.order[getLeftDecendants(tree,tree$Nnode+2)]
#             nodes.2=ident.order[getRightDecendants(tree,tree$Nnode+2)]
#             ident.order = c(nodes.1,nodes.2)
#             
#             if (is.null(genes.use)){
#               # find 
#               genes.use = sample(rownames(object@data),1)
#             }
#             genes.use1 = genes.use[genes.use %in% rownames(object@data)]
#             if (length(genes.use1) < length(genes.use)){
#               print(paste0(length(setdiff(genes.use, genes.use1)), " genes not in data. :")) 
#               print(setdiff(genes.use, genes.use1))
#             }
#             genes.use=genes.use1
#             p <- Perc.pos.by.ident(object,features.use = genes.use,ident.use = ident.order, max.val.exp = max.val.exp, max.val.perc = max.val.perc,return.plot = TRUE,...)
#             print(length(genes.use))
#             multiplot(pdendro,p, cols=1)
#           }
# )

setGeneric("plot_dendro_with_sample_dist", function(object, id.use = "orig", max.val.perc=100,max.size=10,row.scale=FALSE,top.mar=1, left.mar = 4.5,right.mar=1.5, bottom.mar = -1,do.dendro=TRUE,do.transpose=FALSE,ident.order=NULL,x.lab.rot=0,to.return=TRUE,alpha.use=1,...) standardGeneric("plot_dendro_with_sample_dist"))
setMethod("plot_dendro_with_sample_dist","scR",
          function(object, id.use = "orig",max.val.perc=100,max.size=10,row.scale=FALSE,top.mar=1, left.mar = 4.5,right.mar=1.5, bottom.mar = -1,do.dendro=TRUE,do.transpose=FALSE,ident.order=NULL,x.lab.rot=0,to.return=TRUE,alpha.use=1,...) {
            
            if (do.dendro){
              require(ggdendro)
              # Get dendrogram
              tree = object@cluster.tree[[1]]
              data = as.hclust(tree)
              data.tree=as.dendrogram(data)
              ddata <- dendro_data(data.tree, type = "rectangle")
              pdendro <- ggplot(segment(ddata)) + 
                geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
                theme_dendro() + theme(plot.margin = unit(c(top.mar,right.mar,bottom.mar,left.mar), "cm")) # (top, right, bottom, left )
              
              # get node labels according to order
              ident.order=tree$tip.label
              nodes.1=ident.order[getLeftDecendants(tree,tree$Nnode+2)]
              nodes.2=ident.order[getRightDecendants(tree,tree$Nnode+2)]
              ident.order = c(nodes.1,nodes.2)
              print(ident.order)
            } else{
              if (is.null(ident.order)){
                ident.order = levels(object@ident)
              }
              ident.order = rev(ident.order)
              print(ident.order)
            }
            
            
            # Get sample distribution
            ExpMat = table(object@data.info[,id.use],object@ident)
            ExpMat = ExpMat[,ident.order]
            print(colnames(ExpMat))
            if(row.scale){
              ExpMat = t(scale(t(ExpMat), center=FALSE, scale=rowSums(ExpMat)))*100
            } else {
              ExpMat = scale(ExpMat, center=FALSE, scale=colSums(ExpMat))*100
            }
            ExpVal = melt(ExpMat)
            colnames(ExpVal) = c("sample","cluster","perc")
            ExpVal$sample = factor(ExpVal$sample, levels=rev(levels(object@data.info[,id.use])))
            ExpVal$cluster = factor(ExpVal$cluster, levels= ident.order)
            p=ggplot(ExpVal, aes(y = factor(sample),  x = factor(cluster))) + geom_point(aes(colour = perc,  size =perc),alpha=alpha.use) + 
              scale_color_gradient(low ="white",   high = "darkblue", limits=c(min(ExpVal$perc), max(ExpVal$perc) ))+scale_size(range = c(1, max.size))+   theme_bw() +nogrid
            p = p + xlab("Cluster") + ylab("Sample") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
              theme(axis.text.y=element_text(size=12, face="italic")) + theme(axis.text.x = element_text(angle = x.lab.rot))
            
            if (do.dendro){
              multiplot(pdendro,p, cols=1)
            } else {
              
              if (do.transpose){
                ExpVal$sample = factor(ExpVal$sample, levels=levels(object@data.info[,id.use]))
                p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(sample))) + geom_point(aes(colour = perc,  size =perc), alpha=alpha.use) +
                  scale_color_gradient(low ="white",   high = "darkblue", limits=c(min(ExpVal$perc), max(ExpVal$perc) )) + 
                  scale_size(range = c(1, max.size))+   theme_bw() +nogrid
                p = p + xlab("Sample") + ylab("Cluster") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                  theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1))
              }

              print(p)
            }
            
            if (to.return){
              return(p)
            }
              
          }
)

# Reset Par
#
# Reset the graphing space to
# mfrow = c(1, 1)
#
# @param ... Extra parameters for par
#
ResetPar <- function(...) {
  par(mfrow = c(1, 1), ...)
  
}
