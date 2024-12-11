library(ggpubr)
library(ggplot2)
library(scatterpie)
library(RColorBrewer)
library(grDevices)
library(Seurat)
library(viridis)
library(tibble)
library(data.table)
set.seed(1234)
library(dplyr)
library(stringr)

# make sure that rownames(ot) is the barcodes of single-cell reference data
# ref is the seurat Object
# st can be a spatial seurat object or a spatial expression matrix
# provided meta must have cell names as rownames
# default balance is TRUE
# `stGeneExp` function is only for non seuratObject data
# `mapCluster` used for mapping single-cell cluster to corresponding spatial data
# `stClusterPie` used for ploting spatial scatter pie plot
# `stClusterExp` used for ploting single cluster proportion

theme_pie <- function(){
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.position = 'right',
        legend.title=element_text(size=12),
        legend.text = element_text(size = 10))
}
theme_gene <- function(){
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = .5, face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank() #panel.border = element_rect(colour = "black", fill=NA, size=1)
        )
}

theme_cluster <- function(){
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = .5, face = "bold", size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank())
}

mapCluster <- function(ot, ref = NULL, cluster = NULL, meta = NULL, balance = F, min_cut = NULL){
  if(is(ref,"Seurat")){
    meta <- ref@meta.data
  } else {
    if(!is.null(meta)){
      meta <- meta
    } else {
      print('Please provide meta file of scRNA data')
    }
  }
  if(!is.null(cluster)){
    ref_cluster <- meta[cluster]
  } else {
    print('Please provide reference cluster name')
  }
  ref_cluster$cell <- rownames(ref_cluster)
  if(isTRUE(unique(rownames(ot) %in% ref_cluster$cell))){
    ot[] <- apply(ot, 2, as.numeric)
  } else {
    print('Please make sure that rownames of ot is the barcodes of scRNA reference data or rownames of meta file is cell name!')
  }
  if(isTRUE(balance)){
    ot_map <- lapply(unique(ref_cluster[,cluster]),function(x){
      cell = ref_cluster[ref_cluster[,cluster] == x,]$cell
      dt = t(as.data.frame(apply(ot[cell,],2,sum)))
      rownames(dt) <- x
      dt <- dt/length(cell)
      return(dt)
    }) %>%
      do.call('rbind',.) %>% as.data.frame() %>% t()
  } else {
    ot_map <- lapply(unique(ref_cluster[,cluster]),function(x){
      cell = ref_cluster[ref_cluster[,cluster] == x,]$cell
      dt = t(as.data.frame(apply(ot[cell,],2,sum)))
      rownames(dt) <- x
      dt <- dt*(length(cell)/nrow(ref_cluster))
      return(dt)
    }) %>%
      do.call('rbind',.) %>% as.data.frame() %>% t()
  }

  if(is.null(min_cut)){

    return(ot_map)
  } else {
    ot_map <- apply(ot_map,1,function(x){
      x[x < max(x)*min_cut] = 0
      return(x)
    }) %>% as.data.frame() %>% t()
    return(ot_map)
  }
}


stClusterPie <- function(ot_map, st = NULL, slice = NULL, coord = NULL,
                         img_alpha = 0.8, pie_alpha = 0.9, pie_scale = 0.38,
                         color = colorRampPalette(brewer.pal(9,"Set1"))(ncol(cord_ot)-4)){
  if(is(st,"Seurat")){
    if (is.null(slice) && !is.null(names(st@images))){
      slice <- names(st@images)[1]
      warning(sprintf("Using slice %s for ploting", slice))
      img <- st@images[[slice]]@image
    } else {
      img <- st@images[[slice]]@image
    }
    img_grob <- grid::rasterGrob(
      matrix(rgb(img[,,1],img[,,2],img[,,3],alpha = img_alpha),nrow=dim(img)[1]),
      width = unit(1,"npc"))
    cord <- st@images[[slice]]@coordinates[,c('imagerow','imagecol')]
    cord$imagerow_scaled <- cord$imagerow*st@images[[slice]]@scale.factors$lowres
    cord$imagecol_scaled <- cord$imagecol*st@images[[slice]]@scale.factors$lowres
    cord_ot <- cbind(cord,ot_map)

    spatial_scatterPie_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grob = img_grob,
        xmin = 0,
        xmax = ncol(img),
        ymin = 0,
        ymax = -nrow(img)) +
      scatterpie::geom_scatterpie(
        data = cord_ot,
        ggplot2::aes(x = imagecol_scaled,
                     y = imagerow_scaled),
        cols = colnames(cord_ot)[5:ncol(cord_ot)],
        color = NA,
        alpha = pie_alpha,
        pie_scale = pie_scale) +
      ggplot2::ylim(nrow(img),0) +
      ggplot2::xlim(0, ncol(img)) +
      ggplot2::coord_fixed(ratio = 1,xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      scale_fill_manual('cluster',values = color) +
      theme_pie()
    return(spatial_scatterPie_plot)

  } else {
    if(!is.null(coord)){
      colnames(coord) <- c('row_ind','col_ind')
    } else {
      print('Please provide coordinates of st data!')
    }
    cord_ot <- cbind(coord,ot_map)

    spatial_scatterPie_plot <- ggplot() +
      scatterpie::geom_scatterpie(
        aes(x= row_ind, y= col_ind),
        data = cord_ot,
        cols= colnames(cord_ot)[3:ncol(cord_ot)],
        color = NA,
        pie_scale = pie_scale) +
      coord_fixed() + theme(aspect.ratio = 1)+
      #scale_fill_manual('cluster', values = colorRampPalette(brewer.pal(9,"Set1"))(ncol(cord_ot)-2)) +
      scale_fill_manual('cluster', values = color) +
      theme_pie()
    return(spatial_scatterPie_plot)
  }
}


stClusterExp <- function(ot_map, st = NULL, slice = NULL, coord = NULL, cluster, cut = NULL,
                         img_alpha = 0.8, point_size = 0.8, scale.size =  c(2, 6.5)){
  dt.use = t(apply(ot_map, 1, function(x) x/sum(x, na.rm = TRUE)))
  dt.use[is.nan(dt.use)] <- 0
  if(!is.null(cut)){
    dt.use[,cluster][dt.use[,cluster] < max(dt.use[,cluster])*cut] <- NA
  }
  if(is(st,"Seurat")){
    if (is.null(slice) && !is.null(names(st@images))){
      slice <- names(st@images)[1]
      warning(sprintf("Using slice %s for ploting", slice))
      img <- st@images[[slice]]@image
    } else {
      img <- st@images[[slice]]@image
    }
    img_grob <- rasterGrob(
      matrix(rgb(img[,,1],img[,,2],img[,,3],alpha = img_alpha),nrow=dim(img)[1]),
      width = unit(1,"npc"))
    cord <- st@images[[slice]]@coordinates[,c('imagerow','imagecol')]
    cord$imagerow_scaled <- cord$imagerow*st@images[[slice]]@scale.factors$lowres
    cord$imagecol_scaled <- cord$imagecol*st@images[[slice]]@scale.factors$lowres

    cord_ot <- cbind(cord,dt.use)

    spatial_cluster_plot <- ggplot2::ggplot(
      data = cord_ot,aes(x = imagecol_scaled, y = imagerow_scaled, color = get(cluster))) +
      ggplot2::annotation_custom(grob = img_grob, xmin = 0,xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
      ggplot2::geom_point(size = point_size,na.rm = TRUE, colour="grey100") +
      ggplot2::scale_color_gradientn(na.value = "transparent",colors = rev(c('#d7191c','#fdae61','#ffffbf'))) +
      ylim(nrow(img),0) +
      xlim(0, ncol(img)) +
      coord_fixed(ratio = 1,xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
      labs(color = "Proportion") +
      ggtitle(cluster) +
      theme_pie()
    return(spatial_cluster_plot)

  } else {
    if(!is.null(coord)){
      colnames(coord) <- c('row_ind','col_ind')
    } else {
      print('Please provide coordinates of st data!')
    }
    dt.use = t(apply(ot_map, 1, function(x) x/sum(x, na.rm = TRUE)))
    dt.use[is.nan(dt.use)] <- 0
    dt.use <- cbind(coord,dt.use[,cluster,drop=F])
    max.val <- quantile(dt.use[, cluster],0.95, na.rm = TRUE)
    if (max.val > 0) {
      dt.use[, cluster] <- ifelse(dt.use[, cluster] > max.val, max.val, dt.use[, cluster])
    } else {
      max.val <- max(dt.use[, cluster])
    }
    #dt.use[dt.use = 0] <- NA
    spatial_cluster_plot <- ggplot(
      dt.use, aes(row_ind, col_ind, color = get(cluster), size = get(cluster))) +
      geom_point() +
      scale_size(range = scale.size) +
      ggtitle(cluster) +
      scale_color_viridis(alpha = 0.8,breaks = c(min(dt.use[, cluster]),max.val*1.15),
                          labels = c("low","high"),option = "inferno",
                          limits = c(min(dt.use[,cluster])/1.1, max.val*1.2)) +
      theme_gene()+
      coord_fixed() + theme(aspect.ratio = 1)
    return(spatial_cluster_plot)
  }
}

stGeneNorm <- function(st, assay_name = NULL, scale.factor = 10000, norm.method = "LogNormalize"){
  if(is(st,"Seurat") && !is.null(st@assays)){
    if(is.null(assay_name)){
      assay <- DefaultAssay(st)
      warning(sprintf("Using default assay %s", assay))
    } else {
      assay <- assay_name
    }
    stRNA <- NormalizeData(st, normalization.method = norm.method, scale.factor = scale.factor)
    stRNA <- ScaleData(stRNA, features = rownames(stRNA),do.center = F)
  }
  else {
    stRNA <- CreateSeuratObject(counts = st)
    stRNA <- NormalizeData(stRNA, normalization.method = norm.method, scale.factor = scale.factor)
    stRNA <- ScaleData(stRNA, features = rownames(stRNA),do.center = F)
    assay <- DefaultAssay(stRNA)
  }
  dt <- as.data.frame(stRNA@assays$RNA@scale.data)
  return(dt)
}

stGeneExp <- function(exp, gene, coord = NULL, ncol = NULL, size_range = c(2, 6.5), color_alpha = 0.85){
  colnames(coord) <- c("row_ind","col_ind")
  cord_gene <- cbind(coord,t(exp[gene,]))

  if(isTRUE(length(gene) == 1)){
    fig <- ggplot(cord_gene, aes(row_ind, col_ind, color = get(gene), size = get(gene))) +
      geom_point() +
      scale_size(range = size_range) +
      ggtitle(gene) +
      scale_color_viridis(alpha = color_alpha,
                          breaks = c(min(cord_gene[,gene]), max(cord_gene[,gene])),
                          labels = c("low","high")) +
      theme_gene()
    return(fig)
  }
  else {
    if(!is.null(ncol)){
      ncol = ncol
    } else {
      if(length(gene) > 2){
        ncol = 3
      } else {
        ncol = 2
      }
    }
    fig.list <- lapply(colnames(cord_gene)[3:ncol(cord_gene)], function(gene){
      dt.use <- cord_gene
      dt.use <- dt.use[,c("row_ind","col_ind",gene)]
      ggplot(dt.use, aes(row_ind, col_ind, color = get(gene), size = get(gene))) +
        geom_point() +
        scale_size(range = size_range) +
        ggtitle(gene) +
        scale_color_viridis(alpha = color_alpha,
                            breaks = c(min(dt.use[,gene]), max(dt.use[,gene])),
                            labels = c("low","high")) +
        theme_gene()
      })
    fig <- ggpubr::ggarrange(plotlist = fig.list, ncol = ncol)
    return(fig)
  }
}




#' Stacked bar plot showing the proportion of cells across certrain cell groups
#'
#' @param object seurat object
#' @param x Name of one metadata column to show on the x-axis
#' @param fill Name of one metadata column to compare the proportion of cells
#' @param facet Name of one metadata column defining faceting groups
#' @param colors.use defining the color of stacked bar plot; either a char vector defining a color for each cell group or a palette name from brewer.pal
#' @param n.colors Number of colors when setting colors.use to be a palette name from brewer.pal
#' @param n.row Number of rows in facet_grid()
#' @param title.name Name of the main title
#' @param legend.title Name of legend
#' @param xlabel Name of x label
#' @param ylabel Name of y label
#' @param width bar width
#' @param show.legend Whether show the legend
#' @param x.lab.rot Whether rorate the xtick labels
#' @param text.size font size
#' @param flip Whether flip the cartesian coordinates so that horizontal becomes vertical
#' @return ggplot2 object
#' @export
#'
#' @import ggplot2
#' @importFrom plyr ddply as.quoted
computeComposition <- function(object = NULL, y.group = NULL, x.group = NULL, x.levels = NULL, y.levels = NULL,  cutoff.prop = 0.1, type = "heatmap",
                               color.use = NULL, color.heatmap = "RdPu", n.colors = 8, annotation.cols = FALSE, color.use.cols = NULL, group.cols = NULL,
                               dot.size = c(1, 6), angle.x = 45,
                               xlabel = NULL, ylabel = NULL,title.name = NULL, legend.title = NULL,
                               font.size = 10, font.size.title = 12, cluster.rows = TRUE, cluster.cols = FALSE,clustering_distance_rows = "euclidean",
                               x.lab.rot = 45,row.show = NULL, col.show = NULL, remove.isolate = TRUE) {

  # df <- plyr::ddply(object@meta.data, plyr::as.quoted(c(x,fill)), nrow)
  if (!is.null(object)) {
    if (is.character(x.group) | is.null(x.group)) {
      if (is.null(x.group)) {
        x.cell.state <- Seurat::Idents(object)
        x.cell.state.level <- levels(x.cell.state)
      } else if (x.group %in% colnames(object@meta.data) == TRUE) {
        x.cell.state <- object@meta.data[,x.group]
        x.cell.state.level <- levels(x.cell.state)
      } else if (x.group %in% colnames(object@meta.data) == FALSE) {
        stop("'x.group' is not the column of `object@meta.data`! \n")
      }
    } else {
      x.cell.state <- x.group
      x.cell.state.level <- levels(x.cell.state)
    }



    if (is.character(y.group)) {
      if (y.group %in% c("nicheOut", "nicheIn") == TRUE) {
        pred <- t(object@reductions[[y.group]]@cell.embeddings)
        if (y.group == "nicheOut") {
          title.name <- "outgoing"
          ylabel <- "Outgoing communication niches"
        } else if (y.group == "nicheIn"){
          title.name <- "incoming"
          ylabel <- "Incoming communication niches"
        }
      } else if (y.group %in% colnames(object@meta.data) == TRUE) {
        y.cell.state <- object@meta.data[,y.group]
        y.cell.state.level <- levels(y.cell.state)
        pred <- matrix(0, nrow = length(y.cell.state.level), ncol = length(y.cell.state))
        for (i in 1:length(y.cell.state.level)) {
          pred[i, y.cell.state == y.cell.state.level[i]] <- 1
        }
        rownames(pred) <- cell.state.level
      } else {
        stop("'y.group' is not the column of `object@meta.data`! \n")
      }
    } else if (inherits(y.group, what = c("data.frame","matrix", "Matrix", "dgCMatrix"))) {
      pred <- y.group
      cat("The rownames of the input `y.group` are ", toString(rownames(y.group)), " which will be used as group names. \n")
    } else {
      stop("Please check your input `y.group`! \n")
    }

    labels <- x.cell.state
    labels.level <- x.cell.state.level
    prop <- matrix(0, nrow = length(labels.level), ncol = nrow(pred))
    for (i in 1:length(labels.level)) {
      cell.use <- which(labels == labels.level[i])
      prop[i, ] <- apply(pred[, cell.use],1,function(x) thresholdedMean(x, trim = 0.1, na.rm = TRUE))
    }
    colnames(prop) <- rownames(pred)
    rownames(prop) <- labels.level
    prop <- t(prop)

  } else if (inherits(x.group, what = c("data.frame","matrix", "Matrix", "dgCMatrix"))) {
    pred <- x.group # this is the prediction score of each cell type for each spot for 10X visium
    cat("The rownames of the input `x.group` are ", toString(rownames(x.group)), " which will be used as group names. \n")
    if (!is.null(object)) {
      if (y.group %in% colnames(object@meta.data) == TRUE) {
        y.cell.state <- object@meta.data[,y.group]
        y.cell.state.level <- levels(y.cell.state)
      } else {
        stop("When `x.group` is a matrix/dataframe, `y.group` can only be one of the columns of 'object@meta.data'! \n")
      }
    } else {
      y.cell.state <- y.group
      y.cell.state.level <- levels(y.cell.state)
    }

    labels <- y.cell.state
    labels.level <- y.cell.state.level
    prop <- matrix(0, nrow = length(labels.level), ncol = nrow(pred))
    for (i in 1:length(labels.level)) {
      cell.use <- which(labels == labels.level[i])
      prop[i, ] <- apply(pred[, cell.use],1,function(x) thresholdedMean(x, trim = 0.1, na.rm = TRUE))
    }
    colnames(prop) <- rownames(pred)
    rownames(prop) <- labels.level

  } else {
    stop("Please check your input `y.group`! \n")
  }

  mat <- prop
  mat <- sweep(mat, 2L, apply(mat, 2, sum), '/', check.margin = FALSE) # 1 is row, 2 is column
  if (remove.isolate) {
    idx <- which(apply(mat, 2, max) < cutoff.prop)
    if (length(idx) > 0) {
      mat <- mat[, -idx]
    }
  }

  if (!is.null(x.levels)) {
    x.levels <- x.levels[x.levels %in% colnames(mat)]
    mat <- mat[, order(factor(colnames(mat), levels = x.levels)), drop = FALSE]
  }
  if (!is.null(y.levels)) {
    y.levels <- y.levels[y.levels %in% rownames(mat)]
    mat <- mat[order(factor(rownames(mat), levels = y.levels)), , drop = FALSE]
  }

  if (is.null(color.use)) {
    color.use <- CellChat::scPalette(nrow(mat))
  }
  names(color.use) <- rownames(mat)
  if (is.null(color.use.cols)) {
    color.use.cols <- CellChat::scPalette(ncol(mat))
  }
  names(color.use.cols) <- colnames(mat)

  if (type == "dot") {
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("x","y","proportion")
    gg <- ggplot(df, aes(x, y, size = proportion, color = proportion)) +
      geom_point() +
      scale_size_continuous(range = dot.size) +
      theme_linedraw() +
      labs(x = xlabel,
           y = ylabel,
           size = "Proportion") +
      scale_color_viridis_c(option = "D") +
      theme(legend.key.height = grid::unit(0.15, "in"))+
      guides(color = guide_colourbar(barwidth = 0.5, title = "Proportion"))
    # gg <- gg + guides(color = guide_colorbar(barwidth = legend.width, title = "Scaled expression"),size = guide_legend(title = 'Percent expressed'))
    gg <- gg + theme(text = element_text(size = 10),
                     axis.text.x = element_text(angle = x.lab.rot, hjust=1),
                     axis.text.y = element_text(angle = 0, hjust=1),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank()) +
      theme(axis.line.x = element_line(linewidth = 0.25), axis.line.y = element_line(linewidth = 0.25)) +
      theme(panel.grid.major = element_line(colour="grey90", size = (0.1)))
    return(gg)

  } else if (type == "heatmap") {
    if (!is.null(row.show)) {
      mat <- mat[row.show, ]
      color.use <- color.use[row.show]
    }
    if (!is.null(col.show)) {
      mat <- mat[ ,col.show]
      color.use.cols <- color.use.cols[col.show]
    }
    color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)

    df<- data.frame(group = rownames(mat)); rownames(df) <- rownames(mat)
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))
    if (annotation.cols == TRUE) {
      df<- data.frame(group = group.cols); rownames(df) <- colnames(mat)
      col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use.cols), which = "column",
                                          show_legend = FALSE, show_annotation_name = FALSE,
                                          simple_anno_size = grid::unit(0.2, "cm"))
    } else {
      col_annotation <- NULL
    }

    # ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
    # mat[mat == 0] <- NA
    color.heatmap.use = c("white", color.heatmap.use)
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Proportion",
                  left_annotation = row_annotation, bottom_annotation = col_annotation,
                  cluster_rows = cluster.rows,cluster_columns = cluster.cols, clustering_distance_rows = clustering_distance_rows,
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  # width = unit(width, "cm"), height = unit(height, "cm"),
                  row_title = ylabel,row_title_gp = gpar(fontsize = font.size.title),
                  column_title = paste0(" ",title.name, " "),column_title_gp = gpar(fontsize = font.size.title),column_names_rot = x.lab.rot,
                  heatmap_legend_param = list(title = "Proportion", title_gp = gpar(fontsize = font.size, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA,
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
    )

    return(ht1)
  }
}

#' Compute the average expression per cell group when the percent of expressing cells per cell group larger than a threshold
#' @param x a numeric vector
#' @param trim the percent of expressing cells per cell group to be considered as zero
#' @param na.rm whether remove na
#' @return
#' @importFrom Matrix nnzero
#' @export
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) {
    return(0)
  } else {
    return(mean(x, na.rm = na.rm))
  }
}
