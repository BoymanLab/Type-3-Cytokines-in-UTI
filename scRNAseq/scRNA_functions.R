####Introduction####
#Here, we define functions for scRNA-Seq  downstream analysis. 
#Many are inspired or adapted from CATALYST (https://bioconductor.org/packages/release/bioc/html/CATALYST.html)

####Utility Functions####
#scaling function, will be used for heatmaps and feature plots
# scale_exprs <- function(x){
#   margin <- 1
#   q <- 0.01
#   qs <- c(rowQuantiles, colQuantiles)[[margin]]
#   qs <- qs(as.matrix(x), probs = c(0.01, 1-0.01))
#   qs <- matrix(qs, ncol = 2)
#   x <- switch(margin, #margin defines scaling: 1 -> column-wise scaling (without the double transposition)
#               "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
#               "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
#   
#   x[x < 0 | is.na(x)] <- 0 #seems sus, no?
#   x[x > 1] <- 1 #seems sus, no? 
#   return(as.data.frame(t(x)))
# }

scale_exprs <- function(x, margin = 1, q = 0.01) {
  if (!is(x, "matrix")) x <- as.matrix(x)
  qs <- c(rowQuantiles, colQuantiles)[[margin]]
  qs <- qs(x, probs = c(q, 1-q))
  qs <- matrix(qs, ncol = 2)
  x <- switch(margin,
              "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
              "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1
  return(x)
}

#gives barplot annotation for heatmaps
anno_counts <- function(x, perc) {
  ns <- table(x)
  fq <- round(ns/sum(ns)*100, 2)
  if (perc) {
    #txt <- sprintf("%s%%(%s)", fq, names(fq)) #puts the cluster_id in brackets
    txt <- sprintf("%s%%", fq)
    foo <- row_anno_text(txt, 
                         just = "center", 
                         gp = gpar(fontsize = 8),  
                         location = unit(0.5, "npc"))
  } else foo <- NULL
  rowAnnotation(
    "n_cells" = row_anno_barplot(
      x = as.matrix(ns), width = unit(2, "cm"),
      gp = gpar(fill = "grey", col = "white"),
      border = FALSE, axis = TRUE, bar_width = 0.8),
    "foo" = foo)
}

equalize_cells <- function(f){ #f is a data frame with rows for cells and columns for markers and meta data (at least a split_by column, originally e.g. condition)
  minn <- min(table(f$split_by))
  
  namemin <- names(table(f$split_by))[which.min(table(f$split_by))]
  namemax <- names(table(f$split_by))[which.max(table(f$split_by))]
  
  set.seed(42)
  keep <- sample(f[f$split_by == namemax,] %>% rownames(), minn, replace=FALSE)
  keep <- c(keep, f[f$split_by == namemin,] %>% rownames())
  f <- f[keep,]
  return(f)
}

#Takes DGE results from FindMarkers() and converts them to ggpubr compatible format.
#group1 and group2 are the two groups that are tested in the DGE i.e. ident.1 and ident.2 in FindMarkers().
#if the two groups are part of the ordinal scale on the x-axis, no group_by is required. 
#if this is not the case, i.e. if there is a split_by variable and the two groups are separated by dodge, group_by is the column on the x-axis.
dge2ggpubr <- function(dgeres, group1, group2, group_by = NULL){
  formatted <- dgeres %>%
    mutate("group1" = all_of(group1)) %>%
    mutate("group2" = all_of(group2)) %>%
    rownames_to_column("key") %>%
    mutate("p" = signif(p_val, 2)) %>%
    mutate("p.adj" = signif(p_val_adj, 2)) %>%
    dplyr::select("key", "group1", "group2", "p", "p.adj") %>%
    mutate(".y." = "value")
  
  if(!is.null(group_by)){
    formatted$group_by <- all_of(group_by) 
  }
  return(formatted)
}

#Renames the features names of the adt assay of the seaurat object
# RenameGenesSeurat(obj = SeuratObj, newnames = ADTs$name)
RenameADTSeurat <- function(obj = ls.Seurat[[i]], newnames = ADTs$name) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  integrated.adt <- obj@assays$integrated.adt
  if (nrow(integrated.adt) == length(newnames)) {
    if (length(integrated.adt@counts)) integrated.adt@counts@Dimnames[[1]]            <- newnames
    if (length(integrated.adt@data)) integrated.adt@data@Dimnames[[1]]                <- newnames
    if (length(integrated.adt@var.features)) integrated.adt@var.features              <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$integrated.adt <- integrated.adt
  return(obj)
}


#This could be a way: build a wrapper function around this so that it takes gasp as input and returns an rstatix data frame and insert it into switch function
zi.wilcox <- function(gasp = gasp, comparisons = comparisons){
  require(ZIR)
  if(length(group_vars(gasp)) == 1){  #if split_by == NULL and comparisons are made between groups
    test <- data.frame("key" = rep(levels(gasp$key), each = length(comparisons)), 
                       "group1" = rep(as.data.frame(comparisons) %>% t() %>% .[,1], times = length(levels(gasp$key))), 
                       "group2" = rep(as.data.frame(comparisons) %>% t() %>% .[,2], times = length(levels(gasp$key))),
                       "p" = NA,
                       ".y." = "pos.fract") %>%
      dplyr::mutate("groups" = paste(group1, group2))
    for(feature in unique(gasp$key)){
      int <- gasp[gasp$key == feature,]
      for(comp in comparisons){
        test[test$groups == paste(comp %>% unlist, collapse = " ") & test$key == feature, "p"] <- ziw(int[int$group_by == comp[1], "pos.fract"] %>% unlist(),
                                                                                                      int[int$group_by == comp[2], "pos.fract"] %>% unlist())$p.value        
      }
    }
    
  }else if(length(group_vars(gasp)) == 2){#if split_by != NULL and comparisons are made between conditions/timepoints
    for(group in unique(gasp$group_by)){
      for(feature in unique(gasp$key)){
        int <- gasp[gasp$group_by == group & gasp$key == feature,]
        test[test$group_by == group & test$key == feature, "p"] <- ziw(int[int$split_by == levels(int$split_by)[1], "pos.fract"] %>% unlist(),
                                                                      int[int$split_by == levels(int$split_by)[2], "pos.fract"] %>% unlist())$p.value
      }
    }    
  }

  return(test %>% na.omit())
}


####Dimensionality Reduction Functions####
#subsetting the dfs is rather slow, so I combined operations that would require running these steps multiple times.
#this combines: UMAP showing cluster id, and, if features are !NULL, feature plots + heatmap.
clusterplusfeature <- function(x, #a Seurat object
                               protein.features = NULL, #Protein features
                               rna.features = NULL, #RNA features
                               cluster_column = "seurat_clusters", 
                               hml = TRUE, #Heat-map logical: should a heatmap of the clusters be plotted? default is TRUE
                               assay.rna = "integratedRNA", 
                               assay.adt = "integrated.adt",
                               dim = "wnn.umap", 
                               ncol = round(sqrt(length(protein.features)+length(rna.features)))){
  
  #Extract dimensionality reduction embeddings and cluster ids from Seurat object
  u <- x@reductions[[all_of(dim)]]@cell.embeddings %>% as.data.frame()
  colnames(u) <- c("dim1", "dim2")
  u$cluster_id <- x@meta.data[,all_of(cluster_column)] %>% as.factor()
  
  #If Protein or RNA features are given as input, extract the corresponding rows and add them to the DimRed coordinates
  features <- c(protein.features, rna.features)
  if(!is.null(features)){
    if(!is.null(protein.features) & !is.null(rna.features)){
      pf <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      rf <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      f <- cbind(pf, rf) %>% as.data.frame()
    }else if(!is.null(protein.features) & is.null(rna.features)){
      f <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }else if(!is.null(rna.features) & is.null(protein.features)){
      f <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }
    
    #f2 <- scale_exprs(f, margin = 2, q = 0.01) #scale the expression levels of the markers
    u <- cbind(u, f)
  }
  
  clustlabs <- u[c("cluster_id", "dim1", "dim2")]
  clustlabs <- clustlabs %>%
    group_by(cluster_id) %>%
    summarise_all("median")
  
  clusterplot <- ggplot(u, aes(x = dim1, y = dim2))+
    geom_point(aes(col = cluster_id), size = min(4000/nrow(u), 3), shape = 16)+
    theme_bw()+
    scale_color_manual(values = brewer.paired(n = length(unique(u$cluster_id))))+
    geom_text_repel(data = clustlabs, aes(x = dim1, y = dim2, label = cluster_id), size = 5)+
    guides(col = guide_legend(override.aes = list(size = 3)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(clusterplot)
  
  if(!is.null(features)){
    ga <- u %>% gather(key = "key", value = "value", all_of(features))
    ga$key <- factor(ga$key, levels = all_of(features))
    featureplot <- ggplot(ga, aes(x = dim1, y = dim2))+
      geom_point(aes(col = value), size = min(4000/nrow(u), 3), shape = 16)+
      facet_wrap("key", ncol = ncol)+
      theme_bw()+
      scale_color_gradientn(colors = hcl.colors(10, "Viridis"))+
      labs(col = "Scaled expression")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank())
    
    print(featureplot)
    
    if(hml){
      hmdf <- u %>% dplyr::select(-c(dim1, dim2))
      hmdf <- hmdf %>%
        group_by(cluster_id) %>%
        summarise_all("median") %>%
        column_to_rownames("cluster_id")
      
      cluster_colvec <- structure(brewer.paired(n = nrow(hmdf)), names = rownames(hmdf))
      ra <- rowAnnotation(`Cluster ID` = rownames(hmdf),
                          col = list(`Cluster ID` = cluster_colvec), 
                          show_legend = FALSE)
      
      set.seed(1234)
      hm <- Heatmap(as.matrix(hmdf), 
                    col = rev(brewer.pal(11, "RdYlBu")),
                    right_annotation = anno_counts(x@meta.data[,all_of(cluster_column)], perc = TRUE),
                    left_annotation = ra)
      print(hm)
    }
  }
  return(clusterplot)
}
# This one adds a raster function to the featureplot, but doesn't work on Mac
clusterplusfeature <- function(x, #a Seurat object
                               protein.features = NULL, #Protein features
                               rna.features = NULL, #RNA features
                               cluster_column = "seurat_clusters", 
                               hml = TRUE, #Heat-map logical: should a heatmap of the clusters be plotted? default is TRUE
                               assay.rna = "integratedRNA", 
                               assay.adt = "integrated.adt",
                               dim = "wnn.umap", 
                               ncol = round(sqrt(length(protein.features)+length(rna.features)))){
  
  #Extract dimensionality reduction embeddings and cluster ids from Seurat object
  u <- x@reductions[[dim]]@cell.embeddings %>% as.data.frame()
  colnames(u) <- c("dim1", "dim2")
  u$cluster_id <- x@meta.data[,cluster_column] %>% as.factor()
  UTI@assays[["integrated.adt"]]@data@Dimnames[[1]]
  #If Protein or RNA features are given as input, extract the corresponding rows and add them to the DimRed coordinates
  features <- c(protein.features, rna.features)
  if(!is.null(features)){
    if(!is.null(protein.features) & !is.null(rna.features)){
      pf <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      rf <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      f <- cbind(pf, rf) %>% as.data.frame()
    }else if(!is.null(protein.features) & is.null(rna.features)){
      f <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }else if(!is.null(rna.features) & is.null(protein.features)){
      f <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }
    
    #f2 <- scale_exprs(f, margin = 2, q = 0.01) #scale the expression levels of the markers
    u <- cbind(u, f)
  }
  
  clustlabs <- u[c("cluster_id", "dim1", "dim2")]
  clustlabs <- clustlabs %>%
    group_by(cluster_id) %>%
    summarise_all("median")
  
  clusterplot <- ggplot(u, aes(x = dim1, y = dim2))+
    geom_point(aes(col = cluster_id), size = min(4000/nrow(u), 3), shape = 16)+
    theme_bw()+
    scale_color_manual(values = brewer.paired(n = length(unique(u$cluster_id))))+
    geom_text_repel(data = clustlabs, aes(x = dim1, y = dim2, label = cluster_id), size = 5)+
    guides(col = guide_legend(override.aes = list(size = 3)))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  print(clusterplot)
  if(!is.null(features)){
    ga <- u %>% gather(key = "key", value = "value", features)
    ga$key <- factor(ga$key, levels = features)
    featureplot <- ggplot(ga, aes(x = dim1, y = dim2))+
      geom_point_rast(aes(col = value), size = min(4000/nrow(u), 3), raster.dpi = 600, shape = 16)+
      facet_wrap("key", ncol = ncol)+
      theme_bw()+
      scale_color_gradientn(colors = hcl.colors(10, "Viridis"))+
      labs(col = "Scaled expression")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank())
    
    print(featureplot)
    
    if(hml){
      hmdf <- u %>% dplyr::select(-c(dim1, dim2))
      hmdf <- hmdf %>%
        group_by(cluster_id) %>%
        summarise_all("median") %>%
        column_to_rownames("cluster_id")
      
      cluster_colvec <- structure(brewer.paired(n = nrow(hmdf)), names = rownames(hmdf))
      ra <- rowAnnotation(`Cluster ID` = rownames(hmdf),
                          col = list(`Cluster ID` = cluster_colvec), 
                          show_legend = FALSE)
      
      set.seed(1234)
      hm <- Heatmap(as.matrix(hmdf), 
                    col = rev(brewer.pal(11, "RdYlBu")),
                    right_annotation = anno_counts(x@meta.data[,cluster_column], perc = TRUE),
                    left_annotation = ra)
      print(hm)
    }
  }
  return(clusterplot)
}


plot_DR <- function(x, #a Seurat object
                               color_by = "seurat_clusters",
                               hml = FALSE, #Heat-map logical: should a heatmap of the clusters be plotted? default is TRUE
                               assay.rna = "integratedRNA", 
                               assay.adt = "integrated.adt",
                               dim = "wnn.umap", 
                               ncol = round(sqrt(length(color_by)))){
  
  #Extract dimensionality reduction embeddings and cluster ids from Seurat object
  u <- x@reductions[[all_of(dim)]]@cell.embeddings %>% as.data.frame()
  colnames(u) <- c("dim1", "dim2")
  
  #If Protein or RNA features are given as input, extract the corresponding rows and add them to the DimRed coordinates
  if(all(color_by %in% c(rownames(x@assays[[assay.rna]]), rownames(x@assays[[assay.adt]])))){
    if(all(color_by %in% rownames(x@assays[[assay.adt]]))){
      f <- x@assays[[assay.adt]][] %>% as.data.frame() %>% .[all_of(color_by),] %>% t() %>% as.data.frame()
    }else if(all(color_by %in% rownames(x@assays[[assay.rna]]))){
      f <- x@assays[[assay.rna]][] %>% as.data.frame() %>% .[all_of(color_by),] %>% t() %>% as.data.frame()
    }else{
      stop("Currently, only markers from one assay can be plotted.")
      #use grep to find out which markers belong to which assay, extract and combine to one df, see clusterplusfeature
    }
    
    f2 <- scale_exprs(f, margin = 2, q = 0.01) #scale the expression levels of the markers
    u <- cbind(u, f2)
    u <- u %>% gather(key = "key", value = "value", all_of(color_by))
    u$key <- factor(u$key, levels = color_by)
  }else if(color_by %in% colnames(x@meta.data)){
    u$cluster_id <- x@meta.data[,all_of(color_by)] %>% as.factor()
    
    clustlabs <- u[c("cluster_id", "dim1", "dim2")]
    clustlabs <- clustlabs %>%
      group_by(cluster_id) %>%
      summarise_all("median")
  }else{
    stop("Some markers are not available or an inexistent cluster column has been indicated.")
  }
  
  p1 <- ggplot(u, aes(x = dim1, y = dim2))+
    theme_cowplot()
    

  if(all(color_by %in% c(rownames(x@assays[[assay.rna]]), rownames(x@assays[[assay.adt]])))){
    p1 <- p1 + 
      geom_point(aes(col = value), size = min(4000/nrow(u), 3))+
      facet_wrap("key", ncol = ncol)+
      theme(strip.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())+
      scale_color_gradientn(colors = hcl.colors(10, "Viridis"))+
      labs(col = "Scaled expression")
    
  }else if(color_by %in% colnames(x@meta.data)){
    p1 <- p1 + 
      geom_point(aes(col = cluster_id), size = min(4000/nrow(u), 3))+
      geom_text_repel(data = clustlabs, aes(x = dim1, y = dim2, label = cluster_id), size = 5)+
      scale_color_manual(values = brewer.paired(n = length(unique(u$cluster_id))))+
      guides(col = guide_legend(override.aes = list(size = 3)))+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
  }
  
  return(p1)
}


plot_density <- function(x, dim = "UMAP", color_by = "condition", xlim = NULL, ylim = NULL, bins = 10, ncol= NULL){
  
  u <- x@reductions[[all_of(dim)]]@cell.embeddings %>% as.data.frame()
  colnames(u) <- c("dim1", "dim2")
  u <- cbind(u, x@meta.data[all_of(color_by)])
  

  if(is.null(xlim)){
    xlim <- c(min(u$dim1), max(u$dim1)) * 1.15 #fixed the references
  }
  if(is.null(ylim)){
    ylim <- c(min(u$dim2), max(u$dim2)) * 1.15
  }
  
  ggplot(u, aes(x = dim1, y = dim2))+
    geom_density_2d(aes(col = get(color_by)), na.rm = TRUE, bins = bins)+
    theme_bw()+
    facet_wrap(all_of(color_by), ncol = ncol)+
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())+
    scale_color_manual(values = colvec)+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    labs(x = paste(all_of(dim), "1", sep = "_"), y = paste(all_of(dim), "2", sep = "_"), col = color_by)
}

plot_density_fill <- function(x, dim = "UMAP", color_by = "condition", xlim = NULL, ylim = NULL, bins = 10, ncol = NULL){
  
  u <- x@reductions[[all_of(dim)]]@cell.embeddings %>% as.data.frame()
  colnames(u) <- c("dim1", "dim2")
  u <- cbind(u, x@meta.data[all_of(color_by)])
  
  
  if(is.null(xlim)){
    xlim <- c(min(u$dim1), max(u$dim1)) * 1.15
  }
  if(is.null(ylim)){
    ylim <- c(min(u$dim2), max(u$dim2)) * 1.15
  }
  
  ggplot(u, aes(x = dim1, y = dim2))+
    geom_density_2d_filled(na.rm = TRUE, bins = bins)+
    scale_fill_brewer()+
    theme_bw()+
    facet_wrap(all_of(color_by))+
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank())+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    labs(x = paste(all_of(dim), "1", sep = "_"), y = paste(all_of(dim), "2", sep = "_"), col = color_by)
}
#Plot just heatmap of clusters and features
clusthm <- function(x, protein.features = NULL, rna.features = NULL, 
                    cluster_column = "seurat_clusters", 
                    assay.rna = "integratedRNA", assay.adt = "integrated.adt"){

  #If Protein or RNA features are given as input, extract the corresponding rows and add them to the DimRed coordinates
  features <- c(protein.features, rna.features)
  if(!is.null(features)){
    if(!is.null(protein.features) & !is.null(rna.features)){
      pf <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      rf <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t()  %>% scale_exprs(margin = 2, q = 0.01)
      f <- cbind(pf, rf) %>% as.data.frame()
    }else if(!is.null(protein.features) & is.null(rna.features)){
      f <- x@assays[[assay.adt]]@data[protein.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }else if(!is.null(rna.features) & is.null(protein.features)){
      f <- x@assays[[assay.rna]]@data[rna.features,] %>% as.matrix() %>% t() %>% scale_exprs(margin = 2, q = 0.01) %>% as.data.frame()
    }
  }
  f2 <- scale_exprs(f, margin = 2, q = 0.01) %>% as.data.frame() #scale the expression levels of the markers
  
  f2$cluster_id <- x@meta.data[,cluster_column]

  f2 <- f2 %>%
    group_by(cluster_id) %>%
    summarise_all("median") %>%
    column_to_rownames("cluster_id")
  
  cluster_colvec <- structure(brewer.paired(n = nrow(f2)), names = rownames(f2))
  
  ra <- rowAnnotation(`Cluster ID` = rownames(f2),
                      col = list(`Cluster ID` = cluster_colvec), 
                      show_legend = FALSE)
  
  hm <- Heatmap(as.matrix(f2), 
                col = rev(brewer.pal(11, "RdYlBu")),
                right_annotation = anno_counts(x@meta.data[,cluster_column], perc = TRUE),
                left_annotation = ra)
  return(hm)
}

#plot mean expression levels of features by patient


####Functions visualizing selected features####
# DotPlot of RNA expression 
# Here I extract a dataframe from the Seurat object and select only the genes of interest. Then, summarize to get the scaled expression and the fraction of cells expressing it for a dotplot. I make my own dotplot, because the DotPlot function is not exactly what I need
PlotDot <- function(x, #a Seurat object
                    cluster = "primary",
                    split_by = NULL,
                    Assay_counts = "RNA",
                    Assay_scale = "integrated",
                    RNA_markers = rna_state){ #Give a vector containing the genes of interest
  counts.df <- cbind(FetchData(x, vars = if(is.null(split_by)){c(cluster)}
                               else{c(split_by, cluster)}),GetAssayData(x[[Assay_counts]], slot = "counts") %>%
                       as.matrix %>% t() %>% as.data.frame()) %>%
    select(if(is.null(split_by)){cluster
    }else{c(split_by, cluster)},all_of(RNA_markers)) %>% pivot_longer(cols= -if(is.null(split_by)){cluster
    }else{c(split_by, cluster)},
    names_to="gene",
    values_to="measurement") %>%
    group_by(across(.cols = everything()[-ncol(.)])) %>% 
    summarise(meancount = mean(measurement), n= n(), counts = sum(measurement > 0)) %>%
    mutate("Fraction" = counts/n*100)
  
  scale.df <- cbind(FetchData(x, vars = if(is.null(split_by)){c(cluster)}
                              else{c(split_by, cluster)}),GetAssayData(x[[Assay_scale]], slot = "scale.data") %>%
                      as.matrix %>% t() %>% as.data.frame()) %>%
    select(if(is.null(split_by)){cluster
    }else{c(split_by, cluster)},all_of(RNA_markers)) %>% pivot_longer(cols= -if(is.null(split_by)){cluster
    }else{c(split_by, cluster)},
    names_to="gene",
    values_to="measurement") %>%
    group_by(across(.cols = everything()[-ncol(.)])) %>% 
    summarise(ExpressionMean = mean(measurement))
  
  appended.df <- cbind(counts.df, scale.df["ExpressionMean"]) %>% select(-c("meancount","n", "counts"))
  appended.df$gene  <- factor(appended.df$gene, levels = c(RNA_markers))
  if(is.null(split_by)){
    p1 <- ggplot(appended.df, aes(x = get(cluster), y = gene, color = ExpressionMean , size = Fraction)) +
      geom_point() +
      scale_color_gradientn(colors = viridis(9)) +
      scale_size(range = c(1, 30), limits = c(0.01, 100)) +
      theme(panel.spacing = unit(0, "lines"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_line(color = "grey"),  # add this line to set only horizontal grid lines
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 10, face = "bold"),
            strip.text = element_text(size = 14))
    print(p1)
  }else{
    p1 <- ggplot(appended.df, aes(x = get(split_by), y = gene, color = ExpressionMean , size = Fraction)) +
      geom_point() +
      scale_color_gradientn(colors = viridis(9)) +
      facet_wrap(cluster) +
      scale_size(range = c(1, 30), limits = c(0.01, 100)) +
      theme(panel.spacing = unit(0, "lines"),
            panel.background = element_blank(),
            strip.background = element_blank(),
            panel.grid.major = element_line(color = "grey"),  # add this line to set only horizontal grid lines
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 10, face = "bold"),
            strip.text = element_text(size = 14))
    print(p1)
  }
}

plot_pbfeatures <- function(x, pf = NULL, rf = NULL, 
                            assay.rna = "RNA", assay.adt = "Protein", 
                            group.by = "timepoint", bulk.by = "patient_id", paired = TRUE){
  features <- c(pf, rf)
  p <- x@assays[[assay.adt]][] %>% as.data.frame() %>% .[all_of(pf),] %>% t() %>% as.data.frame()
  r <- x@assays[[assay.rna]][] %>% as.data.frame() %>% .[all_of(rf),] %>% t() %>% as.data.frame()
  f <- cbind(p, r)
  f$gp <- x@meta.data[[all_of(group.by)]]
  f$bb <- x@meta.data[[all_of(bulk.by)]]


  f <- f %>% tidyr::gather(key = "key", value = "value", all_of(features))
  
  f <- f %>%
    group_by(gp, bb, key) %>%
    summarise_all("mean")

  f$key <- factor(f$key, levels = features)

  p1 <- ggplot(f, aes(x = gp, y = value))+
    geom_violin(aes(col = gp, fill = gp), alpha = 0.5)+
    geom_point(aes(col = gp))+
    ylim(0, NA)+
    facet_wrap("key", scales = "free_y")+
    theme_bw()+
    labs(fill = "Time point", col = "Time point", y = "Expression")+
    theme(axis.title.x = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_blank())

  if(paired){
    p1 <- p1 + geom_line(aes(group = bb), col = "grey")
  }
  
  p1
  
}

# x <- tregs
# protein.stats <- Tregs.treatment.Protein.markers
# rna.stats <- Tregs.treatment.RNA.markers
#plot expression levels by cell, split by timepoint
plot_features <- function(x, pf = NULL, rf = NULL, 
                          assay.rna = "RNA", assay.adt = "Protein", 
                          evencells = FALSE, 
                          protein.stats = NULL, rna.stats = NULL){
  features <- c(pf, rf)
  p <- x@assays[[assay.adt]][] %>% as.data.frame %>% .[all_of(pf),] %>% t() %>% as.data.frame()
  r <- x@assays[[assay.rna]][] %>% as.data.frame %>% .[all_of(rf),] %>% t() %>% as.data.frame()
  f <- cbind(p, r)
  
  f$sample_id <- x@meta.data$hash.ID
  f <- f %>%
    separate(sample_id, into = c("patient_id", "timepoint"))
  
  ###filter cells such that each time point has the same abundance to make the plots more comparable
  ###if we ever have more than two time points, write sample() as a for loop with for (name in f$timepoint != namemin)....
  if(evencells){
    if(!is.null(split_by)){
      f <- equalize_cells(f)
    }else{
      stop("Please set split_by (e.g. timepoint or condition) for which to equalize cells.")
    }
    
  }

  f <- f %>% gather(key = "key", value = "value", all_of(features))

  
  f$key <- factor(f$key, levels = features)
  
  p1 <- ggplot(f, aes(x = timepoint, y = value))+
    geom_violin(aes(col = timepoint, fill = timepoint), alpha = 0.5)+
    stat_summary(aes(group=1), fun=mean, colour="black", geom="line")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1)+
    geom_jitter(aes(col = timepoint), size = 0.5, height = 0)+
    facet_wrap("key", scales = "free_y")+
    theme_bw()+
    scale_color_manual(values = colvec)+
    scale_fill_manual(values = colvec)+
    labs(fill = "Time point", col = "Time point", y = "Expression")+
    theme(axis.title.x = element_blank(),
          strip.background = element_blank())
  
  stats <- rbind(protein.stats, rna.stats)
  if(!is.null(stats)){
    ypos <- f %>% 
      dplyr::select(key, value) %>% 
      group_by(key) %>%
      summarise(y.position = max(value))
    stat.test <- tibble(".y." = "value", "group1" = "V2", "group2" = "V3", 
                        "p" = stats[rownames(stats) %in% features, "p_val"], 
                        "p.adj" = signif(stats[rownames(stats) %in% features, "p_val_adj"], 2),
                        "key" = factor(stats[rownames(stats) %in% features,] %>% rownames(), levels = features)) %>% 
      merge(ypos)
   p1 <- p1 + stat_pvalue_manual(stat.test, label = "p.adj")+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
   
   # df_fraction <- stats[rownames(stats) %in% features, c("pct.1", "pct.2")] %>% rownames_to_column("marker") %>% tidyr::gather(key = "timepoint", value = "fraction", pct.1, pct.2)
   # df_fraction$timepoint <- ifelse(df_fraction$timepoint == "pct.1", "V3", "V2")
   # 
   # ggplot(df_fraction, aes(x = marker, y = fraction))+
   #   geom_bar(stat = "identity", position = "dodge", aes(fill = timepoint))+
   #   scale_fill_manual(values = colvec)+
   #   theme_bw()+
   #   theme(axis.text.x = element_text(angle = 90),
   #         axis.title.x = element_blank(),
   #         plot.title = element_text(hjust = 0.5))+
   #   labs(y = "Positive Fraction", fill = "Timepoint", title = "Positive Fraction by marker and timepoint")
   
  }
  print(p1)
    
}


plot_features2 <- function(x, pf = NULL, rf = NULL, 
                          group_by = "Treg_phenotpe", split_by = "timepoint",
                          assay.rna = "RNA", assay.adt = "Protein", 
                          evencells = FALSE, dgeres = NULL,
                          comparisons = NULL, method = "wilcox.test", paired = FALSE,
                          write.excel = TRUE, excel.path = "feature_boxplots.xlsx"){
  
  if(!is.null(pf)){
    pf <- x@assays[[assay.adt]][] %>% as.data.frame %>% .[all_of(pf),] %>% t() %>% as.data.frame()
  }
  if(!is.null(rf)){
    rf <- x@assays[[assay.rna]][] %>% as.data.frame %>% .[all_of(rf),] %>% t() %>% as.data.frame()
  }
  f <- merge(pf %>% as.data.frame(), 
             rf %>% as.data.frame(), 
             all = TRUE, by = 0) %>%
    column_to_rownames("Row.names")
  features <- colnames(f)
  
  f$group_by <- x@meta.data[,all_of(group_by)]
  f$cell_id <- colnames(x)
  if(is.factor(x@meta.data[,all_of(group_by)])){
    f$group_by <- factor(f$group_by, levels = levels(x@meta.data[,all_of(group_by)]))
  }
  if(!is.null(split_by)){
    f$split_by <- x@meta.data[,all_of(split_by)]
  }
  
  ###filter cells such that each time point has the same abundance to make the plots more comparable
  ###if we ever have more than two time points, write sample() as a for loop with for (name in f$timepoint != namemin)....
  if(evencells){
    if(!is.null(split_by)){
      f <- equalize_cells(f)
    }else{
      stop("Please set split_by (e.g. timepoint or condition) for which to equalize cells.")
    }
    
  }
  
  f <- f %>% tidyr::gather(key = "key", value = "value", features)
  f$key <- factor(f$key, levels = features)
  
  p1 <- ggplot(f, aes(x = group_by, y = value))+
    geom_boxplot(aes(col = split_by, fill = split_by), alpha = 0.5, outlier.shape = NA)+
    # stat_summary(aes(group=1), fun=mean, colour="black", geom="line")+
    # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1)+
    facet_wrap("key", scales = "free_y")+
    theme_bw()+
    # scale_color_manual(values = colvec)+
    # scale_fill_manual(values = colvec)+
    labs(fill = all_of(split_by), col = all_of(split_by), y = "Expression")+
    theme(axis.title.x = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90))
  
  if(!is.null(split_by)){
    p1 <- p1 + geom_point(aes(fill = split_by), size = 0.5, position = position_jitterdodge(jitter.height = 0, jitter.width = 0.2, dodge.width = 0.75))
  }else{
    p1 <- p1 + geom_point(size = 0.5, position = position_jitter(height = 0, width = 0.2))
  }
  
  if(!is.null(dgeres) & any(features %in% dgeres$key)){
    dgeres <- dgeres[dgeres$key %in% features,]
    if(all(c(dgeres$group1, dgeres$group2) %in% f$split_by)){
      dgeres <- dgeres %>% add_xy_position(x = "group_by", dodge = 0.75, data = f, formula = value ~ split_by)
    }else{ # i think this could be done with some rstatix/ggpubr functions more elegantly, but this works for now.
      dgeres$y.position <- NA
      for(k in unique(dgeres$key)){
        m <- f[f$key == k & f$group_by %in% c(dgeres$group1, dgeres$group2), "value"] %>% max(na.rm = TRUE)
        getypos <- function(m, l){
          cc <- c(m*1.2)
          for(k in seq(2, l)){
            cc <- c(cc, m*(1+k*0.2))
          }
          return(cc)
        }
        ypos <- getypos(m, length(dgeres[dgeres$key == k,"y.position"]))
        dgeres[dgeres$key == k,"y.position"] <- ypos
      }
    }
    
    p1 <- p1 + stat_pvalue_manual(data = dgeres, label = "p.adj")+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }else if(!is.null(comparisons)){
    if(comparisons == "split_by"){
      stat.test <- f %>%
        group_by(group_by, key) %>%
        wilcox_test(value ~ split_by) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj")%>%
        add_significance("p") %>%
        add_xy_position(x = "group_by", dodge = 0.75)
      
      p1 <- p1 + stat_pvalue_manual(data = stat.test)+
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    }else{
      p1 <- p1 + stat_compare_means(comparisons = comparisons, method = method, paired = paired, step.increase = 0.25)
    }
  }
  
  if(write.excel){
    write.xlsx(x = f, file = excel.path)
  }
  
  return(p1)
  
}

get_y_pos <- function(test, data, formula){
  outcome <- deparse(formula[[2]])
  group <- deparse(formula[[3]])
  test$y.position <- NA
  for(k in unique(test$key)){
    test[test$key == k, "y.position"] <- data[data$key == k, outcome] %>% max()
    test[test$key == k, "y.position"] <- test[test$key == k, "y.position"] + (seq(length(test[test$key == k, "y.position"]))-1) * (test[test$key == k, "y.position"][1]/5) + 0.03
  }
  return(test)
}

plot_features3 <- function(x, protein.features = NULL, rna.features = NULL, 
                           group_by = "Treg_phenotpe", split_by = "timepoint",
                           assay.rna = "RNA", assay.adt = "Protein", 
                           evencells = FALSE, pointsize = 0.5,
                           comparisons = NULL, #three options: Data frame with DGE results (dge2ggpubr output, can combine multiple dges), comparison list, or "split_by"
                           method = "wilcox.test", paired = FALSE,
                           write.excel = TRUE, excel.path = "feature_boxplots.xlsx"){
  
  if(!is.null(protein.features)){
    protein.features <- x@assays[[assay.adt]][all_of(protein.features),] %>% as.data.frame() %>% t() %>% as.data.frame()
  }
  if(!is.null(rna.features)){
    rna.features <- x@assays[[assay.rna]][all_of(rna.features),] %>% as.data.frame() %>% t() %>% as.data.frame()
  }
  f <- merge(protein.features %>% as.data.frame(), 
             rna.features %>% as.data.frame(), 
             all = TRUE, by = 0) %>%
    column_to_rownames("Row.names")
  features <- colnames(f)
  
  f$group_by <- x@meta.data[,all_of(group_by)]
  f$cell_id <- colnames(x)
  if(is.factor(x@meta.data[,all_of(group_by)])){
    f$group_by <- factor(f$group_by, levels = levels(x@meta.data[,all_of(group_by)]))
  }
  if(!is.null(split_by)){
    f$split_by <- x@meta.data[,all_of(split_by)]
  }
  
  ###filter cells such that each time point has the same abundance to make the plots more comparable
  ###if we ever have more than two time points, write sample() as a for loop with for (name in f$timepoint != namemin)....
  if(evencells){
    if(!is.null(split_by)){
      f <- equalize_cells(f)
    }else{
      stop("Please set split_by (e.g. timepoint or condition) for which to equalize cells.")
    }
    
  }
  
  f <- f %>% tidyr::gather(key = "key", value = "value", features)
  f$key <- factor(f$key, levels = features)
  
  p1 <- ggplot(f, aes(x = group_by, y = value))+
    geom_boxplot(aes(col = split_by, fill = split_by), alpha = 0.5, outlier.shape = NA)+
    # stat_summary(aes(group=1), fun=mean, colour="black", geom="line")+
    # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.1)+
    facet_wrap("key", scales = "free_y")+
    theme_bw()+
    # scale_color_manual(values = colvec)+
    # scale_fill_manual(values = colvec)+
    labs(fill = all_of(split_by), col = all_of(split_by), y = "Expression")+
    theme(axis.title.x = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90))
  
  ##Add points
  if(!is.null(split_by)){
    p1 <- p1 + geom_point(aes(fill = split_by), size = pointsize, position = position_jitterdodge(jitter.height = 0, jitter.width = 0.2, dodge.width = 0.75))
  }else{
    p1 <- p1 + geom_point(size = pointsize, position = position_jitter(height = 0, width = 0.2))
  }
  
  ##Add stats
  if(!is.null(comparisons)){ #check if statistics should be performed at all
    if(is.data.frame(comparisons)){ #check if dge results are provided to annotate statistics
      if(any(features %in% comparisons$key)){
        if(all(c(comparisons$group1, comparisons$group2) %in% f$split_by)){ #comparison is data frame and groups are split_by
          #Filter dge-results for plottet features; add factor levels; arrange by factors (this is essential for y-position of pval)
          comp <- comparisons %>% dplyr::filter(key %in% features)
          comp$key <- factor(comp$key, levels = features)
          comp$group_by <- factor(comp$group_by, levels = levels(f$group_by))
          comp <- comp %>% arrange(key, group_by)
          
          #Add xy-position
          stat.test <- comp %>%
            add_significance("p.adj") %>%
            add_x_position(x = "group_by", group = "split_by", dodge = 0.75) %>%
            add_y_position(data = f %>% group_by(key, group_by), formula = value ~ split_by, step.increase = 0)
  
          p1 <- p1 + stat_pvalue_manual(stat.test, label = "p.adj")+
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
        }else if(all(c(comparisons$group1, comparisons$group2) %in% f$group_by)){
          comp <- comparisons %>% dplyr::filter(key %in% features)
          comp$key <- factor(comp$key, levels = features)
          comp <- comp %>% arrange(key)
  
          stat.test <- comp %>%
            add_significance("p.adj") %>%
            add_x_position(x = "group_by")
          stat.test <- get_y_pos(test = stat.test, data = f %>% group_by(key), formula = value ~ group_by)
  
          p1 <- p1 + stat_pvalue_manual(stat.test, label = "p.adj")+
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
        }
      }
    }else if(is.list(comparisons)){ #check if list is provided to calculate stats using ggpubr (for values on ordinal x-scale, i.e. no split_by-variables)
      p1 <- p1 + stat_compare_means(comparisons = comparisons, method = method, paired = paired, step.increase = 0.25)
    }else if(comparisons == "split_by"){ #option to calculate stats using ggpubr for split_by variables
      stat.test <- f %>%
        group_by(group_by, key) %>%
        wilcox_test(value ~ split_by) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj")%>%
        add_significance("p") %>%
        add_xy_position(x = "group_by", dodge = 0.75)
      
      p1 <- p1 + stat_pvalue_manual(data = stat.test)+
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
      
    }
  }
  
  if(write.excel){
    write.xlsx(x = f, file = excel.path)
  }
  
  return(p1)
  
}

### Attention when integrating this into plot_features: probably cannot be run together with even_cells == TRUE
#x <- tregs


fractions2 <- function(x, pf = NULL, rf = NULL, 
                           group_by = "cluster_id", bulk_by = "patient_id", split_by = "condition",
                           assay.rna = "RNA", assay.adt = "Protein", comparisons = NULL, method = "t.test", paired = FALSE, write.excel = FALSE, 
                       excel.path = "feature_fraction_boxplots.xlsx"){
  
  if(!is.null(pf)){
    pf <- x@assays[[assay.adt]][] %>% as.data.frame %>% .[all_of(pf),] %>% t() %>% as.data.frame()
  }
  if(!is.null(rf)){
    rf <- x@assays[[assay.rna]][] %>% as.data.frame %>% .[all_of(rf),] %>% t() %>% as.data.frame()
  }
  f <- merge(pf %>% as.data.frame(), 
             rf %>% as.data.frame(), 
             all = TRUE, by = 0) %>%
    column_to_rownames("Row.names")
  features <- colnames(f)
  f[f > 0] <- 1
  
  f$group_by <- x@meta.data[,all_of(group_by)]
  f$bulk_by <- x@meta.data[,all_of(bulk_by)]
  if(!is.null(split_by)){
    f$split_by <- x@meta.data[,all_of(split_by)]
  }
  
  ga <- f %>% 
    rownames_to_column("cell_id") %>% 
    gather(key = "key", value = "value", all_of(features)) %>%
    dplyr::select(-cell_id)
  
  ga <- ga %>% 
    table() %>% 
    as.data.frame()

  gasp <- ga %>% 
    spread(key = "value", value = "Freq") %>%
    dplyr::rename("neg" = `0`, "pos" = `1`) %>%
    mutate("pos.fract" = pos/(neg + pos))
  gasp[is.na(gasp$pos.fract), "pos.fract"] <- 0
  
  p1 <- ggplot(gasp, aes(x = group_by, y = pos.fract))+
    geom_boxplot(aes(col = split_by))+
    geom_point(aes(fill = split_by), position = position_dodge(width = 0.75), size = 0.5)+
    facet_wrap("key")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          strip.background = element_blank())
  
  if(!is.null(split_by)){
    p1 <- p1 + geom_point(aes(fill = split_by), size = 0.5, position =position_dodge(width = 0.75))
  }else{
    p1 <- p1 + geom_point(size = 0.5)
  }

  
  if(!is.null(comparisons)){
    if(comparisons == "split_by"){
      gasp <- gasp %>%
        group_by(group_by, key)
      stat.test <- switch(method,
               "t.test" = t_test(pos.fract ~ split_by, data = gasp),
               "wilcox.test" = wilcox_test(pos.fract ~ split_by, data = gasp)) %>% #we now have the ZIwilcox function, just needs some rewriting to be compatible with ggpubr
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj")%>%
        add_significance("p") %>%
        add_xy_position(x = "group_by", dodge = 0.75)
      
      p1 <- p1 + stat_pvalue_manual(data = stat.test)+
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    }else if(is.list(comparisons)){
      gasp <- gasp %>%
        group_by(key)
      stat.test <- switch(method,
                          "t.test" = t_test(pos.fract ~ group_by, data = gasp),
                          "wilcox.test" = wilcox_test(pos.fract ~ group_by, data = gasp),
                          "ziwilcox.test" = zi.wilcox(gasp = gasp, comparisons = comparisons)) %>% #we now have the ZIwilcox function, just needs some rewriting to be compatible with ggpubr
        adjust_pvalue(method = "fdr") %>%
        add_significance("p.adj") %>%
        add_significance("p") %>%
        add_xy_position(x = "group_by", dodge = 0.75, data = gasp, formula = pos.fract ~ group_by)
      
      p1 <- p1 + stat_pvalue_manual(data = stat.test)+
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))    }
  }
  
  if(write.excel){
    write.xlsx(x = f, file = excel.path)
  }
  
  return(p1)
}

#### DGE functions ####
#I was thinking whether edgeR could be a better option for scRNA DGE. Robinson wrote this nature methods on the topic: https://doi.org/10.1038/nmeth.4612
#TLDR: I think using wilcoxon is fine.

plot_volcano <- function(dgeres = dgeres, ylims = c(0,NA), xlims = "free", labels = NULL, suffix = "-RNA"){
  dgeres <- dgeres %>%
    na.omit() %>%
    rownames_to_column("gene_nomen") %>%
    dplyr::mutate("avgExprs" = (pct.1 + pct.2)/2)

  if(!is.null(suffix)){
    dgeres <- dgeres %>% dplyr::mutate("gene_nomen" = str_replace_all(gene_nomen, all_of(suffix), ""))
  }
  dgeres$delabel <- NA
  if(!is.null(labels)){
    dgeres[dgeres[,"gene_nomen"] %in% all_of(labels), "delabel"] <- dgeres[dgeres[,"gene_nomen"] %in% all_of(labels), "gene_nomen"]
  }else{
    dgeres[dgeres$p_val_adj < 0.05 & abs(dgeres$avg_log2FC) > 1.0, "delabel"] <- dgeres[dgeres$p_val_adj < 0.05 & abs(dgeres$avg_log2FC) > 1.0, "gene_nomen"]
    dgeres[duplicated(dgeres$delabel), "delabel"] <- NA #this should not be necessary if the same genes are used to summarize and to plot    
  }
  dgeres$col <- "NO"
  dgeres$col[dgeres$p_val_adj < 0.05 & dgeres$avg_log2FC > 1.0] <- "UP"
  dgeres$col[dgeres$p_val_adj < 0.05 & dgeres$avg_log2FC < -1.0] <- "DOWN"
  
  p1 <- ggplot(dgeres, aes(x = avg_log2FC, y = -log10(p_val_adj)))+
    geom_point(aes(fill = col, col = col, size = avgExprs, alpha = 0.8), shape = 16)+
    scale_color_manual(breaks = c("DOWN", "NO", "UP"), values = c("blueviolet", "grey", "seagreen"))+
    labs(x = "Log2 FC", y = "-log10 FDR", col = "Regulation", size = "Positive Cell Fration")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(limits = ylims)
  
  if(!is.null(labels)){
    p1 <- p1 + geom_text_repel(aes(label = delabel), size = 8,  fontface = "bold", min.segment.length = unit(0.1, 'lines'), max.overlaps = 10000, max.time = 5, max.iter = 100000)
  }else{
    p1 <- p1 + geom_text_repel(aes(label = delabel), size = 8,  fontface = "bold", min.segment.length = unit(0.1, 'lines'))
  }
  
  if(xlims == "symmetrical"){
    p1 <- p1 + scale_x_continuous(limits = c(-max(abs(dgeres$avg_log2FC)), max(abs(dgeres$avg_log2FC))))
  }else if(xlims != "free"){
    p1 <- p1 + scale_x_continuous(limits = xlims)
  }
  
  return(p1)
}

####GSEA functions####
run_fgsea <- function(dge, gene_set, suffix = NULL, sampleSize = 1000){
  
  fgsea_sets<- gene_set %>% split(x = .$gene_symbol, f = .$gs_name)
  
  Markers <- dge %>% 
    na.omit() %>%
    mutate(ranking_metric = -log(p_val)*sign(avg_log2FC)) %>%
    arrange(desc(ranking_metric)) %>% 
    rownames_to_column(var = "feature")
  
  if(!is.null(suffix)){
    Markers$feature <- Markers$feature %>% str_replace_all(pattern = suffix, "")
  } 
  
  preranked_list <- Markers %>% 
    dplyr::select(feature, ranking_metric) %>% 
    deframe()
  
  set.seed(42)
  fgseaResults <- fgsea(fgsea_sets, stats = preranked_list, sampleSize = sampleSize) %>%
    as_tibble() %>%
    arrange(desc(abs(NES))) %>% 
    mutate("LeadingEdge" = NA)
  
  for(q in 1:nrow(fgseaResults)){
    fgseaResults$LeadingEdge[q] <- fgseaResults$leadingEdge[q] %>% unlist() %>% .[1:6] %>% na.omit() %>% paste(collapse = ", ") %>% unlist()
  }
  fgseaResults <- fgseaResults %>% dplyr::select(-leadingEdge)
  
  return(fgseaResults)
}

gsea_barplot <- function(fgsea_res, top_n = NULL, p_max = 0.05, pathways = NULL){
  df <- fgsea_res %>% arrange(padj)
  if(!is.null(p_max)){
    df <- df %>% dplyr::filter(padj < p_max)
  }
  if(!is.null(pathways)){
    df <- df %>% dplyr::filter(pathway %in% pathways)
    df$pathway <- factor(df$pathway, levels = rev(pathways))}else {
    df <- df[1:top_n,]}
    
  
  ggplot(df, aes(x = pathway, y = -log10(padj)))+
    geom_bar(stat = "identity", fill = "black", col = "black", width = 0.7)+
    theme_classic()+
    coord_flip()+
    geom_text(aes(label = LeadingEdge), hjust = -0.05, col = "black")+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.7)))+
    labs(y = "-log10 FDR")+
    theme(axis.title.y = element_blank())+
    scale_x_discrete(limits = rev)  
}

plot_fgsea <- function(dge, gene_set, top_n = NULL, pathways = NULL,
                            gseaParam=1,
                            ticksSize=0.2, suffix = NULL, sampleSize = 1000) {
  
  fgsea_sets<- gene_set %>% split(x = .$gene_symbol, f = .$gs_name)
  
  Markers <- dge %>% 
    na.omit() %>%
    mutate(ranking_metric = -log(p_val)*sign(avg_log2FC)) %>%
    arrange(desc(ranking_metric)) %>% 
    rownames_to_column(var = "feature")
  
  if(!is.null(suffix)){
    Markers$feature <- Markers$feature %>% str_replace_all(pattern = suffix, "")
  } 
  
  preranked_list <- Markers %>% 
    dplyr::select(feature, ranking_metric) %>% 
    deframe()
  
  set.seed(42)
  fgseaResults <- fgsea(fgsea_sets, stats = preranked_list, sampleSize = sampleSize) %>%
    as_tibble() %>%
    arrange(desc(abs(NES)))

  rnk <- rank(-preranked_list)
  ord <- base::order(rnk)
  
  statsAdj <- preranked_list[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  if(!is.null(top_n)){
    if(top_n<nrow(fgseaResults)){
      genesets <- fgseaResults[1:top_n, "pathway"] %>% unlist() %>% unname()
    }else{
      genesets <- fgseaResults[, "pathway"] %>% unlist() %>% unname()
    }
  }else if(!is.null(pathways)){
    genesets <- pathways
  }else{
    genesets <- fgseaResults[fgseaResults$padj <= 0.05, "pathway"] %>% unlist() %>% unname()
  }
  plotlist <- list()
  for(pw in genesets){
    pathway <- unname(as.vector(na.omit(match(fgsea_sets[[pw]], names(statsAdj)))))
    pathway <- unique(pathway)
    pathway <- sort(pathway)
    
    
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)
    
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
    
    diff <- (max(tops) - min(bottoms)) / 8
    
    x=y=NULL
    plotlist[[pw]] <- ggplot(toPlot, aes(x=x, y=y)) +
      #geom_point(color=linecolor, size=0.1) +
      geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
      geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
      geom_hline(yintercept=0, colour="black") +
      geom_line(color="green", size = 0.8) + theme_classic() +
      geom_segment(data=data.frame(x=pathway),
                   mapping=aes(x=x, y=-diff/2,
                               xend=x, yend=diff/2),
                   size=ticksSize) +
      
      theme(panel.border=element_blank(),
            panel.grid.minor=element_blank(),
            title= element_text(size = 18, face = "bold"),
            axis.text.x = element_text(size = 14),# I've added some text size changes
            axis.text.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18)) +
      labs(x="Rank", y="Enrichment score", title = pw)+
      annotate(geom = "text",
               size = 8,
               x = length(preranked_list) * 0.7,
               y = max(tops) * 0.8,
               #vjust = 20,
               label = paste0("p = ", signif(fgseaResults[fgseaResults$pathway == pw, "padj"], 2),"\n", "NES = ", round(fgseaResults[fgseaResults$pathway == pw, "NES"], 2)))
  }# I've added the NES score as well
  return(plotlist)
}

run_gsva <- function(x, cluster_column = "Treg_phenotpe", assay = "RNA", genesets = Hallmarks, suffix = "-RNA", verbose = TRUE, ...){
  df <- x@assays[[all_of(assay)]][]  %>% as.data.frame() %>% t() %>% as.data.frame() %>% .[,!colSums(.) == 0]
  df$cluster <- x@meta.data[,all_of(cluster_column)]
  colnames(df) <- str_replace_all(colnames(df), all_of(suffix), "")
  df <- df %>% as.data.table()
  if(verbose){message("aggregating data")}
  df <- setDT(df)[, lapply(.SD, mean), by = .(cluster)] %>% as.data.frame()

  df <- df %>% column_to_rownames("cluster")

  gs <- genesets %>% split(x = .$gene_symbol, f = .$gs_name)
  
  if(verbose){message("running GSVA")}
  gsva_res <- gsva(t(df), gs, ...)
  
}


####Clustering functions####
flowsom <- function(x, pf = NULL, rf = NULL, assay.adt = "integrated.adt", assay.rna = "integratedRNA", maxK = 20){
  df.adt <- x@assays[[assay.adt]][] %>% as.data.frame() %>% .[all_of(pf),] %>% as.data.frame()
  df.rna <- x@assays[[assay.rna]][] %>% as.data.frame() %>% .[all_of(rf),] %>% as.data.frame()
  
  df <- rbind(df.adt, df.rna)
  
  fsom <- ReadInput(flowFrame(t(df)))
  
  xdim <- 10
  ydim <- 10
  set.seed(1234)
  som <- BuildSOM(fsom = fsom, xdim = xdim, ydim = ydim)
  maxK <- maxK
  mc <- suppressWarnings(suppressMessages(
    ConsensusClusterPlus(t(som$map$codes), 
                         maxK = maxK, reps = 100, 
                         distance = "euclidean", plot = "pdf", seed = 1234)))
  
  codes <- data.frame(seq_len(xdim*ydim), map(mc[-1], "consensusClass"))
  codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
  colnames(codes) <- c(sprintf("som%s", xdim*ydim), sprintf("meta%s", seq_len(maxK)[-1]))
  
  flowsom <- merge(data.frame(cluster.id = som$map$mapping[,1], order = 1:length(som$map$mapping[,1])), 
                   data.frame(seq_len(100), map(mc[-1], "consensusClass"))[c(1,maxK)], by.x = "cluster.id", by.y = "seq_len.100.")
  flowsom <- flowsom[order(flowsom$order),]
  x@meta.data$flowsom <- flowsom[,3]
  
  clustree(codes, "meta") %>% print()
  return(x)
}

phenograph <- function(x, pf = NULL, rf = NULL, assay.adt = "integrated.adt", assay.rna = "integratedRNA", k = 30){
  df.adt <- x@assays[[assay.adt]][] %>% as.data.frame() %>% .[all_of(pf),] %>% as.data.frame()
  df.rna <- x@assays[[assay.rna]][] %>% as.data.frame() %>% .[all_of(rf),] %>% as.data.frame()
  
  df <- rbind(df.adt, df.rna)
  
  set.seed(42)
  phenoresults <- Rphenograph(t(df), k = k)
  x@meta.data$pheno <- as.factor(membership(phenoresults[[2]]))
  
  return(x)
}

merge_clusters <- function(x, cluster_column, merging_table, merging_name = "merging1", levels = NULL){
  x@meta.data$order <- 1:nrow(x@meta.data)
  coldata <- x@meta.data[c("order", cluster_column)]
  colnames(merging_table) <- c(cluster_column, merging_name)
  coldata <- merge(coldata, merging_table)
  coldata <- coldata[order(coldata$order),]
  
  x@meta.data[,merging_name] <- coldata[,merging_name]
  if(is.null(levels)){
    x@meta.data[,merging_name] <- factor(x@meta.data[,merging_name])
  }else if(!is.null(levels)){
    x@meta.data[,merging_name] <- factor(x@meta.data[,merging_name], levels = levels)
  }
  
  Idents(object = x) <- x@meta.data[merging_name]
  
  return(x)
}

####Archive####
fractions <- function(x, pf = NULL, rf = NULL, 
                      assay.rna = "RNA", assay.adt = "Protein"){
  features <- c(pf, rf)
  p <- x@assays[[assay.adt]][] %>% as.data.frame %>% .[all_of(pf),] %>% t() %>% as.data.frame()
  r <- x@assays[[assay.rna]][] %>% as.data.frame %>% .[all_of(rf),] %>% t() %>% as.data.frame()
  f <- cbind(p, r)
  
  f$sample_id <- x@meta.data$hash.ID
  f <- f %>%
    separate(sample_id, into = c("patient_id", "timepoint"))
  
  
  
  ga <- f %>% rownames_to_column("cell_id") %>% gather(key = "key", value = "value", all_of(features))
  
  
  ga$key <- factor(ga$key, levels = features)
  ga[ga$value > 0, "value"] <- 1
  tbl <- table(ga$value, ga$patient_id, ga$timepoint, ga$key) %>% as.data.frame()
  tbl <- tbl %>% spread(key = "Var1", value = "Freq")
  colnames(tbl) <- c("patient_id", "timepoint", "marker", "negative", "positive")
  tbl$fraction <- tbl$positive / (tbl$negative + tbl$positive)
  tbl$mp <- paste(tbl$marker, tbl$patient_id)
  
  p1 <- ggplot(tbl, aes(x = marker, y = fraction))+
    #geom_boxplot(aes(col = timepoint))+
    geom_bar(stat = "summary", fun = "mean", position = "dodge", aes(fill = timepoint), alpha = 0.5)+
    geom_point(aes(col = timepoint), position = position_dodge(width = 0.75))+
    #geom_line(aes(group = interaction(marker, patient_id)), col = "grey", position = position_dodge(width = 0.75))+
    #lemon::geom_pointline(aes(col = timepoint, group = interaction(marker, patient_id)), position = position_dodge(width = 0.75))+
    scale_color_manual(values = colvec)+
    scale_fill_manual(values = colvec)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    labs(y = "Positive Fraction", col = "Timepoint", fill = "Timepoint", title = "Positive Fraction by marker and timepoint")
  print(p1)
  
  compare_means(formula = fraction~timepoint, data = tbl %>% dplyr::filter(marker == "IL2RA-Protein"), method = "t.test", paired = TRUE)
  
  p2 <- ggplot(tbl, aes(x = timepoint, y = fraction))+
    #geom_boxplot(aes(col = timepoint))+
    geom_bar(stat = "summary", fun = "mean", position = "dodge", aes(fill = timepoint), alpha = 0.5)+
    geom_point(aes(col = timepoint), position = position_dodge(width = 0.75))+
    geom_line(aes(group = patient_id), col = "grey")+
    scale_color_manual(values = colvec)+
    scale_fill_manual(values = colvec)+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())+
    labs(y = "Positive Fraction", col = "Timepoint", fill = "Timepoint", title = "Positive Fraction by marker and timepoint")+
    facet_wrap("marker", nrow = 1)+
    stat_compare_means(comparisons = list(c("V2", "V3")), paired = TRUE, method = "t.test", tip.length = 0.01, p.adjust.method = "holm")
  print(p2)
}
