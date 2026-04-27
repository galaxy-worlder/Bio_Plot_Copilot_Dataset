DOTplot<-function (object, assay = NULL, features, cols = c("lightgrey", 
                                                   "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, 
          dot.scale = 6, idents = NULL, group.by = NULL, split.by = NULL, 
          cluster.idents = FALSE, scale = TRUE, scale.by = "radius", 
          scale.min = NA, scale.max = NA) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", 
              call. = FALSE, immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", 
                                         color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", 
                     no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) +theme_bw()+
    theme(panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          #panel.border =element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 20, hjust = 0.5),
          title = element_text(size = 15),
          legend.text =element_text(size=10),  # Font size of legend labels.
          legend.key.size=unit(0.2, "inches")
    )+theme(
      axis.line = element_line(colour = 'black', size = 0.5), 
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12,color = 'black'),
      axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
    )

  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups,labeller = label_wrap_gen(10), scales = "free_x", 
                              space = "free_x", switch = "y") + 
      theme(panel.spacing = unit(x = 0.25, units = "lines"), 
            #geom_text(aes(label = feature.groups), x = Inf, y = Inf, hjust = 1.5, vjust = 1.5),
            strip.background  = element_rect(fill = "#5BBCD6",colour = "black", size = 0),
            strip.text =element_text(size = 12,color = 'black')#,strip.background = element_blank()
            )
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

FEATUREplot<-function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("lightgrey", "#ff0000", "#00ff00")
} else {
  c("lightgrey", "blue")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
reduction = NULL, split.by = NULL, keep.scale = "feature", 
shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, 
label = FALSE, label.size = 4, label.color = "black", 
repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE, 
sort.cell = NULL, interactive = FALSE, combine = TRUE, raster = NULL, 
raster.dpi = c(512, 512)) {
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ", 
            "parameter instead for equivalent functionality.", 
            call. = FALSE, immediate. = TRUE)
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(object = object, feature = features[1], 
                        dims = dims, reduction = reduction, slot = slot))
  }
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature","all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
                                                                                           size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                             default.colors[1], "' for double-negatives", 
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, "ident", 
                                              features), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", 
         paste(features, collapse = ", "), " in slot ", 
         slot, call. = FALSE)
  }else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", 
         call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
  ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- SetQuantile(cutoff = min.cutoff[index - 
                                                                                    3], data.feature)
                                       max.use <- SetQuantile(cutoff = max.cutoff[index - 
                                                                                    3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  }else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells, 
                                                            drop = TRUE], object[[split.by, drop = TRUE]][cells, 
                                                                                                          drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, 
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, 
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], 
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, 
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, 
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ", 
             paste(no.expression, collapse = ", "), 
             call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")], 
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, 
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      }
      else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", 
                                   feature, shape.by)]
      plot <- SingleDimPlot(data = data.single, dims = dims, 
                            col.by = feature, order = order, pt.size = pt.size, 
                            cols = cols.use, shape.by = shape.by, label = FALSE, 
                            raster = raster, raster.dpi = raster.dpi) + scale_x_continuous(limits = xlims) + 
        scale_y_continuous(limits = ylims) + theme_bw()+
        theme(panel.grid.major =element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              panel.border = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              title = element_text(size = 15),
              legend.text =element_text(size=10),  # Font size of legend labels.
              #legend.title = element_blank(), 
              legend.key.size=unit(0.2, "inches")
        )+theme(
          axis.line = element_line(colour = 'black', size = 0.5), 
          axis.title = element_text(size = 15, color = 'black'),
          axis.text = element_text(size = 12)+NoLegend()
        )+labs(x = 'UMAP1', y = "UMAP2")+CenterTitle()
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident", 
                              repel = repel, size = label.size, color = label.color)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, 
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident), 
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(), 
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(), 
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", 
                    unique.feature.exp, ") of ", feature, 
                    ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                                                                       guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && 
          !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, 
                                                              limits = c(min.feature.value, max.feature.value)))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots, 
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")), 
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2], 
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol, 
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  }else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() + 
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""), 
                                                                   limits = ylims) + no.right))+
          theme(plot.title = element_text(family="serif",size = 20, hjust = 0.5),
                 axis.text = element_text(family="serif",size = 12, color = 'black'), 
                 axis.title = element_text(family="serif",size = 20, color = 'black'),
                 legend.text=element_text(family="serif",size = 12, color = 'black'),
                 legend.title=element_text(family="serif",size = 15, color = 'black'))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), 
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% 
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" && 
        !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, 
                                                              limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}


VINPLOT<-function (object, features, cols = NULL, pt.size = NULL, idents = NULL, 
          sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL, 
          adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE, 
          ncol = NULL, slot = "data", split.plot = FALSE, stack = FALSE, 
          combine = TRUE, fill.by = "feature", flip = FALSE, 
          raster = NULL) {
  if (!is.null(x = split.by) & getOption(x = "Seurat.warn.vlnplot.split", 
                                         default = TRUE)) {
    message("The default behaviour of split.by has changed.\n", 
            "Separate violin plots are now plotted side-by-side.\n", 
            "To restore the old behaviour of a single split violin,\n", 
            "set split.plot = TRUE.\n      \nThis message will be shown once per session.")
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot(object = object, type = ifelse(test = split.plot, 
                                                yes = "splitViolin", no = "violin"), features = features, 
                 idents = idents, ncol = ncol, sort = sort, assay = assay, 
                 y.max = y.max, same.y.lims = same.y.lims, adjust = adjust, 
                 pt.size = pt.size, cols = cols, group.by = group.by, 
                 split.by = split.by, log = log, slot = slot, stack = stack, 
                 combine = combine, fill.by = fill.by, flip = flip, raster = raster))
}

DimPlot<-function (plot, id, clusters = NULL, labels = NULL, split.by = NULL, 
                   repel = TRUE, box = FALSE, geom = "GeomPoint", position = "median", 
                   ...) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), 
                    use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, 
            " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, 
                                                                         id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, 
                                                                            id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", 
         paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == "GeomSpatial") {
    xrange.save <- layer_scales(plot = plot)$x$range$range
    yrange.save <- layer_scales(plot = plot)$y$range$range
    data[, xynames["y"]] = max(data[, xynames["y"]]) - 
      data[, xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - 
        pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + 
        sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- data[data[, id] == group, , drop = FALSE]
    data.medians <- if (!is.null(x = split.by)) {
      do.call(what = "rbind", args = lapply(X = unique(x = data.use[, 
                                                                    split.by]), FUN = function(split) {
                                                                      medians <- apply(X = data.use[data.use[, split.by] == 
                                                                                                      split, xynames, drop = FALSE], MARGIN = 2, 
                                                                                       FUN = median, na.rm = TRUE)
                                                                      medians <- as.data.frame(x = t(x = medians))
                                                                      medians[, split.by] <- split
                                                                      return(medians)
                                                                    }))
    }
    else {
      as.data.frame(x = t(x = apply(X = data.use[, xynames, 
                                                 drop = FALSE], MARGIN = 2, FUN = median, na.rm = TRUE)))
    }
    data.medians[, id] <- group
    data.medians$color <- data.use$color[1]
    return(data.medians)
  })
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) == 
                           as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1, 
                                                                               2)]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, 
                                                                        id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels), 
         ") must be equal to the number of clusters being labeled (", 
         length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel, 
                       no = geom_label)
    plot <- plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames["x"], 
                                                                    y = xynames["y"], label = id, fill = id), show.legend = FALSE, 
                            ...) + scale_fill_manual(values = labels.loc$color[order(labels.loc[, 
                                                                                                id])])
  }
  else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, 
                       no = geom_text)
    plot <- plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames["x"], 
                                                                    y = xynames["y"], label = id), show.legend = FALSE, 
                            ...)
  }
  if (geom == "GeomSpatial") {
    plot <- suppressMessages(expr = plot + coord_fixed(xlim = xrange.save, 
                                                       ylim = yrange.save))
  }
  return(plot)
}







