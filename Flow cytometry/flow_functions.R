
#### Utility functions ####
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

set_path <- function(...){
  s <- file.path(...)
  s <- str_replace_all(s, "//", "/")
  return(s)
}

#reproducible downsampling
#There are other functions to down-sample data sets, I have found the one used in the CytoNorm pipeline not to be reproducible, thus we use this for-loop.

downsample <- function(fs, n = 10000){
  dsFilt <- sampleFilter(size = n, filterId="dsFilter")
  result <- flowCore::filter(fs, dsFilt)
  
  for(i in 1:length(result)){
    a <- result@.Data[[i]]@subSet
    l <- a %>% length()
    if(n>l){
      next
    }
    set.seed(42)
    ids <- sample(1:l, n, replace=FALSE) %>% sort()
    x <- rep(FALSE, l)
    x[ids] <- TRUE
    result@.Data[[i]]@subSet <- x
  }
  
  fs.ds <- Subset(fs, result) 
  return(fs.ds)
}


#### Clustering algorithms ####
phenograph <- function(x = sce, markers = type_markers(x), assay = "exprs", k = 30, seed = 1234, scale = TRUE){
  set.seed(seed)
  df <- x %>% assay("exprs") %>% .[markers,]
  
  if(scale){
    df <- scale_exprs(df, margin = 1, q = 0.01)
  }
  
  phenoresults <- Rphenograph(t(df), k = k)
  
  colData(x)$pheno <- as.factor(membership(phenoresults[[2]]))
  x$cluster_id <- as.factor(membership(phenoresults[[2]]))
  metadata(x)$cluster_codes <- data.frame(
    custom = factor(levels(as.factor(membership(phenoresults[[2]]))), levels = levels(as.factor(membership(phenoresults[[2]])))))
  
  return(x)
}


AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x){
    cbind(get(x),source = x)
  }))
}

flowsom <- function(x = sce, markers = type_markers(x), assay = "exprs", maxK = 20, seed = 1234, xdim = 10, ydim = 10, scale = TRUE, plot.tree = TRUE){
  df <- x %>% assay(all_of(assay)) %>% .[all_of(markers),]
  
  if(scale){
    df <- scale_exprs(df, margin = 1, q = 0.01)
  }
  
  fsom <- ReadInput(flowFrame(t(df)))
  
  xdim <- xdim
  ydim <- ydim
  set.seed(seed = seed)
  som <- BuildSOM(fsom = fsom, xdim = xdim, ydim = ydim)
  maxK <- maxK
  mc <- suppressWarnings(suppressMessages(
    ConsensusClusterPlus(t(som$map$codes), 
                         maxK = maxK, reps = 100, 
                         distance = "euclidean", plot = "pdf", seed = seed)))
  
  codes <- data.frame(seq_len(xdim*ydim), map(mc[-1], "consensusClass"))
  codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
  colnames(codes) <- c(sprintf("som%s", xdim*ydim), sprintf("meta%s", seq_len(maxK)[-1]))
  x@metadata$cluster_codes <- codes
  x@metadata$SOM_codes <- som$data
  
  flowsom <- merge(data.frame("cluster.id" = som$map$mapping[,1], "order" = 1:length(som$map$mapping[,1])), 
                   data.frame(seq_len(xdim*ydim), map(mc[-1], "consensusClass"))[c(1,maxK)], by.x = "cluster.id", by.y = "seq_len.xdim...ydim.")
  flowsom <- flowsom[order(flowsom$order),]
  
  colData(x)$cluster_id <- as.factor(flowsom[,1])
  colData(x)$flowsom <- factor(flowsom[,3], levels = c(1:maxK))
  
  if(plot.tree){
    clustree(codes, "meta") %>% print()
  }
  
  return(x)
}


#### Differential Analysis ####
plot_DA <- function(x = sce, cluster_column = "primary", subject_id = "patient_id", condition = "condition", comparisons = NULL, method = "t.test", paired = FALSE){
  dtp <- colData(x)[c(cluster_column, condition, subject_id)] %>%
    table() %>% 
    as.data.frame() %>%
    dplyr::rename("counts" = Freq) %>%
    mutate("sample_id" = paste(subject_id, condition, sep = "_"))
  
  
  totdtp <- colData(x)[all_of(c(subject_id, condition))] %>%
    table() %>% 
    data.frame() %>% 
    dplyr::rename("totcounts" = Freq)%>%
    mutate("sample_id" = paste(all_of(subject_id), all_of(condition), sep = "_"))
  
  dtp <- merge(dtp, totdtp) %>% mutate("Freq" = counts/totcounts*100)
  
  p1 <- ggplot(dtp, aes(x = get(condition), y = Freq, col = get(condition)))+
    geom_violin(aes(fill = get(condition)), alpha = 0.5, draw_quantiles = 0.5)+
    #geom_jitter(width = 0.1)+
    geom_point()+
    facet_wrap(cluster_column, scales = "free_y")+
    theme_bw()+
    labs(col = all_of(condition), fill = all_of(condition), title = paste0("Cluster Abundance by ", all_of(condition)), y = "Frequency (%)")+
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank())
  
  if(paired == TRUE){
    p1 <- p1 + geom_line(aes(group = get(subject_id)), col = "grey")
    
  }
  
  if(!is.null(comparisons)){
    p1 <- p1 + 
      stat_compare_means(comparisons = comparisons, method = method, paired = paired)+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    
  }
  return(p1)
}


calculate_DA <- function(x = sce, cluster_column = "primary", sample_id = "sample_id", condition = "condition", timepoint = "timepoint"){
  
  dtp <- colData(x)[all_of(c(cluster_column, condition, timepoint, sample_id))] %>%
    table() %>% 
    as.data.frame() %>%
    merge(ei(x)) %>% 
    select(-"n_cells") %>% 
    dplyr::rename("counts" = Freq)
  
  totdtp <- colData(x)[all_of(c(sample_id, condition, timepoint))] %>%
    table() %>% 
    data.frame() %>% 
    merge(ei(x)) %>% 
    select(-"n_cells") %>% 
    dplyr::rename("totcounts" = Freq)
  
  df <- merge(dtp, totdtp) %>% dplyr::mutate("Freq" = counts/totcounts)
  
  return(df)
}

# plot Differential State analysis as boxplot facet wrapped by marker
plot_DS <- function(x = sce, assay = "exprs", markers = state_markers(x), cluster_column = "primary", subject_id = "patient_id", condition = "condition", timepoint = "timepoint"){
  
  mfidf <- cbind(assay(x, all_of(assay)) %>% t() %>% as.data.frame() %>% .[all_of(markers)], 
                 colData(x) %>% as.data.frame() %>% .[all_of(c(cluster_column, subject_id, condition))])
  mfidfga <- mfidf %>% 
    #dplyr::rename("cluster_column" = all_of(cluster_column), "condition" = all_of(condition), "subject_id" = all_of(subject_id)) %>%
    dplyr::group_by(across(all_of(c(condition, subject_id, cluster_column)))) %>% 
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    gather(key = "marker", value = "mfi", all_of(markers))
  
  p1 <- ggplot(mfidfga, aes(x = get(cluster_column), y = mfi, col = get(condition), fill = get(condition)))+
    geom_boxplot(alpha = 0.5)+
    geom_point(position = position_dodge(width = 0.75))+
    #geom_line(aes(group = patient_id), col = "grey")+
    facet_wrap("marker", scales = "free_y")+
    theme_bw()+
    labs(col = all_of(condition), fill = all_of(condition), y = "MFI")+
    theme(axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  return(p1)
}

calculate_DS <- function(x = sce, assay = "exprs", markers = state_markers(x), cluster_column = "primary", subject_id = "patient_id", condition = "condition", timepoint= "timepoint"){
  
  mfidf <- cbind(assay(x, all_of(assay)) %>% t() %>% as.data.frame() %>% .[all_of(markers)], 
                 colData(x) %>% as.data.frame() %>% .[all_of(c(cluster_column, subject_id, condition, timepoint))])
  mfidfga <- mfidf %>% 
    #dplyr::rename("cluster_column" = all_of(cluster_column), "condition" = all_of(condition), "subject_id" = all_of(subject_id)) %>%
    dplyr::group_by(across(all_of(c(condition, subject_id, cluster_column, timepoint)))) %>% 
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    gather(key = "marker", value = "mfi", all_of(markers))

  return(mfidfga)
}

#### Dimensionality reduction & Visualization ####
run_DR <- function(x, 
                   dr = c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap"), 
                   cells = NULL, markers = type_markers(x), assay = "exprs", scale = TRUE) {
  
  # check validity of input arguments
  dr <- match.arg(dr, c("UMAP", "TSNE", "PCA", "MDS", "DiffusionMap"))
  
  if (is.null(cells)) {
    # use all cells
    cs <- TRUE 
  } else {
    if (is.null(x$sample_id))
      stop("colData column sample_id not found,\n ", 
           " but is required to downsample cells.")
    stopifnot(
      is.numeric(cells), length(cells) == 1,
      as.integer(cells) == cells, cells > 0)
    # split cell indices by sample
    cs <- split(seq_len(ncol(x)), x$sample_id)
    # sample at most 'n' cells per sample
    cs <- unlist(lapply(cs, function(u)
      sample(u, min(cells, length(u)))))
  }
  
  # run dimension reduction
  df <- x %>% assay(all_of(assay)) %>% .[all_of(markers), cs]
  
  if(scale){
    df <- scale_exprs(df, margin = 1, q = 0.01)
  }
  
  fun <- get(paste0("calculate", dr))
  embeddings <- fun(df)
  
  # return SCE when no cell subsetting has been done
  if (is.null(cells)){
    reducedDim[x, all_of(dr)] <- embeddings
    return(x)
  } else{  
    # else, match coordinates to subsetted cells from SCE
    m <- matrix(NA, nrow = ncol(x), ncol = ncol(embeddings))
    m[cs, ] <- embeddings
    reducedDim(x, all_of(dr)) <- m
    return(x)
  }
}

# Plot clusters in UMAP space
plot_clusters <- function(x = sce, color_by = "pheno", dimred = "UMAP", add_labels = TRUE){
  df <- cbind(colData(x) %>% as.data.frame() %>% dplyr::select(all_of(color_by)), 
              reducedDim(x, all_of(dimred)) %>% as.data.frame()) %>% na.omit()
  
  if(add_labels){
    clustlabs <- df[c(all_of(color_by), "V1", "V2")]
    clustlabs <- clustlabs %>%
      group_by(across(all_of(color_by))) %>%
      summarise_all("median")
  }

  
  p1 <- ggplot(df, aes(x = V1, y = V2))+
    geom_point(aes(col = get(color_by)), size = min(4000/nrow(df), 3))+
    theme_minimal()+
    scale_color_manual(values = brewer.paired(n = length(unique(df[,all_of(color_by)]))))+
    guides(col = guide_legend(override.aes = list(size = 3)))+
    labs(col = color_by, x = paste(dimred, "1", sep = "_"), y = paste(dimred, "2", sep = "_"))
  
  if(add_labels){
    p1 <- p1 + geom_text_repel(data = clustlabs, aes(x = V1, y = V2, label = get(color_by)), size = 5)
  }
  return(p1)
}

# Plot features in UMAP space
plot_features <- function(x = sce, features = type_markers(x), dimred = "UMAP", assay = "exprs"){
  df <- cbind(assay(x, all_of(assay)) %>% as.data.frame() %>% .[all_of(features),] %>% scale_exprs(margin = 1) %>% t() %>% as.data.frame(), 
              reducedDim(x, all_of(dimred)) %>% as.data.frame()) %>% na.omit()
  
  
  ga <- df %>% gather(key = "key", value = "value", all_of(features))
  ga$key <- factor(ga$key, levels = all_of(features))
  featureplot <- ggplot(ga, aes(x = V1, y = V2))+
    geom_point(aes(col = value), size = min(4000/nrow(df), 3))+
    facet_wrap("key")+
    theme_minimal()+
    scale_color_gradientn(colors = hcl.colors(10, "Viridis"))+
    labs(col = "Scaled expression")
  return(featureplot)
}

#combining plot_clusters and plot_features, should be the most flexible of these functions
plot_DR <- function(x = sce, color_by = "pheno", dimred = "UMAP", add_labels = TRUE, assay = "exprs", ncol = NULL, nrow = NULL){
  if(all(all_of(color_by) %in% rownames(x))){
    a <- assay(x, all_of(assay)) %>% as.data.frame() %>% .[all_of(color_by),] %>% scale_exprs(margin = 1) %>% t() %>% as.data.frame()
    add_labels <- FALSE
  }else if(length(all_of(color_by)) == 1 & all_of(color_by) %in% colnames(colData(x))){
    a <- colData(x) %>% as.data.frame() %>% dplyr::select(all_of(color_by)) %>% dplyr::rename("color_by" = all_of(color_by))
  }else(
    stop("color_by: Please specify any number of features or one column of colData")
  )
  
  b <- reducedDim(x, all_of(dimred)) %>% as.data.frame()
  
  df <- cbind(a, b) %>% na.omit()
  
  if(add_labels){
    clustlabs <- df[c("color_by", "V1", "V2")]
    clustlabs <- clustlabs %>%
      group_by(color_by) %>%
      summarise_all("median")
  }
  
  #if > 1 features are used for color_by, gather features
  if(length(all_of(color_by)) > 1){
    df <- df %>% gather(key = "key", value = "color_by", all_of(color_by))
    df$key <- factor(df$key, levels = all_of(color_by))
  }
  
  
  p1 <- ggplot(df, aes(x = V1, y = V2))+
    ## rasterize points to avoid overplotting, on mac this sometimes gives an error.
    geom_point_rast(aes(col = color_by), size = min(4000/nrow(df), 3), raster.dpi = 600)+
    theme_minimal()+
    labs(x = paste(dimred, "1", sep = "_"), y = paste(dimred, "2", sep = "_"))
  
  if(add_labels){
    p1 <- p1 + geom_text_repel(data = clustlabs, aes(x = V1, y = V2, label = color_by), size = 6)
  }
  
  if(length(all_of(color_by)) > 1){
    p1 <- p1 + facet_wrap("key", ncol = ncol, nrow = nrow)
  }
  
  if(all(all_of(color_by) %in% rownames(x))){
    p1 <- p1 + scale_color_gradientn(colors = hcl.colors(10, "Viridis")) + labs(col = "scaled expression")
  }else{
    p1 <- p1 + 
      #scale_color_manual(values = custom_colors)+
      scale_color_manual(values = brewer.paired(n = length(unique(df[,"color_by"]))))+
      guides(col = guide_legend(override.aes = list(size = 3)))
  }
  
  return(p1)
}


calc_dist <- function(){
  a <- df[2,1] - df[1,1]
  b <- df[2,2] - df[1,2]
  c <- sqrt(a^2+b^2)
  return(c)
}

calc_density <- function(df){
  dist_mtrx <- matrix(nrow = nrow(df), ncol = nrow(df))
  for(i in 1:(nrow(df)-1)){
    for(j in (i+1):nrow(df)){
      a <- df[j,1]-df[i,1]
      b <- df[j,2]-df[i,2]
      c <- sqrt(a^2+b^2)
      dist_mtrx[i,j] <- c
      dist_mtrx[j,i] <- c
    }
  }
  dist_mtrx[is.na(dist_mtrx)] <- 0
  dist_mtrx[dist_mtrx<1] <- 1
  dist_mtrx[dist_mtrx>1] <- 0
  df$nn <- rowSums(dist_mtrx)
  return(df)
}

plot_density <- function(x = sce, dimred = "UMAP", color_by = "condition", xlim = NULL, ylim = NULL, bins = 10, colvec = NULL, contures = TRUE, ncol = NULL){
  
  df <- cbind(reducedDim(x, all_of(dimred)), colData(x)[all_of(color_by)]) %>% as.data.frame() %>% na.omit()
  
  if(is.null(xlim)){
    xlim <- c(min(df$V1), max(df$V1)) * 1.15
  }
  if(is.null(ylim)){
    ylim <- c(min(df$V2), max(df$V2)) * 1.15
  }
  
  p1 <- ggplot(df, aes(x = V1, y = V2))+
    theme_bw()+
    facet_wrap(all_of(color_by), ncol = ncol)+
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    # scale_color_manual(values = brewer.paired(n = length(unique(df[,all_of(color_by)]))))+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    labs(x = paste(dimred, "1", sep = "_"), y = paste(dimred, "2", sep = "_"), col = color_by)
  
  if(contures){
    p1 <- p1 + geom_density_2d(aes(col = get(color_by)), na.rm = TRUE, bins = bins, linewidth = 0.2)
    if(!is.null(colvec)){
      p1 <- p1 + scale_color_manual(values = colvec)
    }else{
      p1 <- p1 + scale_color_brewer(palette = "Set1")
    }
  }else{ #calculate density: define radius, count cells in this radius for each cell and set color relative to total cell number
    df <- calc_density(df)
    p1 <- p1 + geom_point(data = df, aes(col = nn)) + scale_color_gradient(low = "blue", high = "red")
    #p1 + geom_hex(bins = 20) + scale_fill_gradient(low = "blue", high = "red")
  }

  return(p1)
}



#### Normalization ####
CytoNorm.train <- function(files,
                           labels,
                           channels,
                           transformList,
                           outputDir = "./tmp",
                           FlowSOM.params = list(nCells = 1000000,
                                                 xdim = 15,
                                                 ydim = 15,
                                                 nClus = 10,
                                                 scale = FALSE),
                           normMethod.train = QuantileNorm.train,
                           normParams = list(nQ = 101),
                           seed = NULL,
                           clean = TRUE,
                           plot = FALSE,
                           verbose = FALSE,
                           ...){
  
  if (length(labels) != length(files)) {
    stop("Input parameters 'labels' and 'files'",
         " should have the same length")
  }
  
  # Create output directory
  dirCreated = FALSE
  if (!dir.exists(outputDir)) {
    dirCreated = dir.create(outputDir)
  }
  
  if(!file.exists(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))){
    
    nCells <- FlowSOM.params[["nCells"]]
    if(is.null(FlowSOM.params[["channels"]])){
      FlowSOM.channels <- channels
    } else {
      FlowSOM.channels <- FlowSOM.params[["channels"]]
    }
    FlowSOM.params <- FlowSOM.params[grep("nCells|channels",
                                          names(FlowSOM.params),
                                          invert = TRUE)]
    fsom <- prepareFlowSOM(files = files,
                           nCells = nCells,
                           FlowSOM.params = FlowSOM.params,
                           transformList = transformList,
                           colsToUse = FlowSOM.channels,
                           seed = seed,
                           ...)
    
    saveRDS(fsom, file.path(outputDir, "CytoNorm_FlowSOM.RDS"))
    
    if (plot) {
      FlowSOM::FlowSOMmary(fsom,
                           plotFile = file.path(outputDir, "CytoNorm_FlowSOM.pdf"))
    }
  } else {
    fsom <- readRDS(file.path(outputDir, "CytoNorm_FlowSOM.RDS"))
    warning("Reusing previously saved FlowSOM result.")
  }
  
  # Split files by clusters
  for(i in seq_along(files)) {
    if(is(files, "flowSet")) {
      file <- sampleNames(files)[i]
      ff <- files[[i]]
    } else {
      file <- files[i]
      ff <- flowCore::read.FCS(file, ...) 
    }
    if(verbose) message("Splitting ", file)
    if (!is.null(transformList)) {
      ff <- flowCore::transform(ff,transformList)
    }
    # Map the file to the FlowSOM clustering
    fsom_file <- FlowSOM::NewData(fsom, ff)
    # Get the metacluster label for every cell
    cellClusterIDs <- FlowSOM::GetMetaclusters(fsom_file) #fsom$metaclustering[GetClusters(fsom_file)]
    for (cluster in unique(fsom$metaclustering)) {
      if (sum(cellClusterIDs == cluster) > 0) {
        suppressWarnings(
          flowCore::write.FCS(
            ff[cellClusterIDs == cluster,],
            file=file.path(outputDir,
                           paste0(gsub("[:/]","_",file),
                                  "_fsom",cluster,".fcs"))))
      }
    }
  }
  
  # file names
  if(is(files, "flowSet")) {
    file_names <- sampleNames(files)
  } else {
    file_names <- files
  }
  
  # Learn quantiles for each cluster
  clusterRes <- list()
  for (cluster in unique(fsom$metaclustering)) {
    if(verbose) message("Processing cluster ",cluster)
    if (plot) {
      grDevices::pdf(file.path(outputDir,
                               paste0("CytoNorm_norm_Cluster",
                                      cluster, ".pdf")),
                     height = 3*(2*length(files)+2),
                     width = 3*(length(channels)+1))
    }
    
    normParams_tmp <- c(normParams,
                        list(files = file.path(outputDir,
                                               paste0(gsub("[:/]", "_", file_names),
                                                      "_fsom", cluster, ".fcs")),
                             labels = as.character(labels),
                             channels = channels,
                             transformList = NULL,
                             verbose = verbose,
                             plot = plot))
    normParams_tmp <- normParams_tmp[unique(names(normParams_tmp))]
    if(is.list(normParams[["goal"]])){
      normParams_tmp[["goal"]] <- normParams[["goal"]][[cluster]]
    }
    clusterRes[[cluster]] <- do.call(normMethod.train,
                                     normParams_tmp)
    
    if (plot) { grDevices::dev.off() }
  }
  
  if(clean){
    for(cluster in unique(fsom$metaclustering)){
      tmp_files <- file.path(outputDir,
                             paste0(gsub("[:/]", "_", file_names),
                                    "_fsom", cluster, ".fcs"))
      
      file.remove(tmp_files[file.exists(tmp_files)])
    }
    if(dirCreated & !plot){
      unlink(outputDir, recursive=TRUE)
    }
  }
  named.list(fsom, clusterRes)
}




CytoNorm.normalize <- function(model,
                               files,
                               labels,
                               transformList,
                               transformList.reverse,
                               outputDir = ".",
                               prefix = "Norm_",
                               clean = TRUE,
                               verbose = FALSE,
                               normMethod.normalize = QuantileNorm.normalize,
                               write = TRUE,
                               ...){
  if(is.null(model$fsom) |
     is.null(model$clusterRes)){
    stop("The 'model' paramter should be the result of using the
             trainQuantiles function.")
  }
  
  if(length(labels) != length(files)){
    stop("Input parameters 'labels' and 'files' should have the same length")
  }
  
  # Create output directory
  if(!dir.exists(outputDir)){
    dir.create(outputDir)
  }
  
  fsom <- model$fsom
  clusterRes <- model$clusterRes
  
  # file names
  if(is(files, "flowSet")) {
    file_names <- sampleNames(files)
  } else {
    file_names <- files
  }
  
  # Split files by clusters
  cellClusterIDs <- list()
  meta <- list()
  cluster_files <- list()
  for(i in seq_along(files)) {
    if(is(files, "flowSet")) {
      file <- file_names[i]
      ff <- files[[i]]
    } else {
      file <- files[i]
      ff <- flowCore::read.FCS(file, ...)
    }
    if(verbose) message("Splitting ",file)
    if(!is.null(transformList)){
      ff <- flowCore::transform(ff, transformList)
      # meta[[file]] <- list()
      # meta[[file]][["description_original"]] <- ff@description
      # meta[[file]][["parameters_original"]] <- ff@parameters
    }
    
    fsom_file <- FlowSOM::NewData(fsom,ff)
    
    cellClusterIDs[[file]] <- FlowSOM::GetMetaclusters(fsom_file)
    
    for(cluster in unique(fsom$metaclustering)){
      if (sum(cellClusterIDs[[file]] == cluster) > 0) {
        f <- file.path(outputDir,
                       paste0(gsub("[:/]","_",file),
                              "_fsom", cluster, ".fcs"))
        suppressWarnings(
          flowCore::write.FCS(ff[cellClusterIDs[[file]] == cluster],
                              file = f)
        )
      }
    }
  }
  
  # Apply normalization on each cluster
  for(cluster in unique(fsom$metaclustering)){
    if(verbose) message("Processing cluster ",cluster)
    files_tmp <- file.path(outputDir,
                           paste0(gsub("[:/]",
                                       "_",
                                       file_names),
                                  "_fsom",
                                  cluster,
                                  ".fcs"))
    labels_tmp <- labels[file.exists(files_tmp)]
    files_tmp <- files_tmp[file.exists(files_tmp)]
    normMethod.normalize(model = clusterRes[[cluster]],
                         files = files_tmp,
                         labels = labels_tmp,
                         outputDir = file.path(outputDir),
                         prefix = "Norm_",
                         transformList = NULL,
                         transformList.reverse = NULL,
                         removeOriginal = TRUE,
                         verbose = verbose)
  }
  
  # Combine clusters into one final fcs file
  res <- flowSet(
    lapply(
      seq_along(files),
      function(i) {
        
        if(is(files, "flowSet")) {
          file <- file_names[i]
          ff <- files[[i]]
        } else {
          file <- files[i]
          ff <- flowCore::read.FCS(file, ...)
        }
        if(verbose) message("Rebuilding ",file)
        for(cluster in unique(fsom$metaclustering)){
          file_name <- file.path(outputDir,
                                 paste0("Norm_",gsub("[:/]","_",file),
                                        "_fsom",cluster,".fcs"))
          if (file.exists(file_name)) {
            ff_subset <- flowCore::read.FCS(file_name, ...)
            flowCore::exprs(ff)[cellClusterIDs[[file]] == cluster,] <- flowCore::exprs(ff_subset)
          }
        }
        if(!is.null(transformList.reverse)){
          ff <- flowCore::transform(ff, transformList.reverse)
          # ff@description <- meta[[file]][["description_original"]]
          # ff@parameters <- meta[[file]][["parameters_original"]]
        }
        
        
        # Adapt to real min and max because this gets strange values otherwise
        ff@parameters@data[,"minRange"] <- apply(ff@exprs, 2, min)
        ff@parameters@data[,"maxRange"] <- apply(ff@exprs, 2, max)
        ff@parameters@data[,"range"] <- ff@parameters@data[,"maxRange"] -
          ff@parameters@data[,"minRange"]
        
        if(clean){
          file.remove(file.path(outputDir,
                                paste0("Norm_",gsub("[:/]","_",file),
                                       "_fsom",unique(fsom$metaclustering),".fcs")))
        }
        
        if(write) {
          suppressWarnings(
            flowCore::write.FCS(
              ff,
              file= file.path(outputDir,paste0(prefix,gsub(".*/","",file)))
            )
          )
        }
        return(ff)
      }
    )
  )
  
  # dremove empty output directory
  if(length(list.files(outputDir)) == 0){
    unlink(outputDir)
  }
  
  # normalized flowSet
  return(res)
  
}



save_flow_set <- function(fs = fs, md = md, folder_col = "batch", file_name_col = "file_name", out.dir = path.files){
  if(!exists(file.path(out.dir))){dir.create(file.path(out.dir))}
  
  for(a in md[,all_of(folder_col)] %>% unique()){
    if(!exists(file.path(out.dir, a))){dir.create(file.path(out.dir, a))}
      
    ffs <- md[md[all_of(folder_col)] == a, all_of(file_name_col)]
    write.flowSet(fs[ffs,], outdir = file.path(out.dir, a))
  }
}


roberts_flowset_saver <- function(fs = fs, md = md, 
                                  folder1 = "organ", folder2 = "batch", 
                                  file_name_col = "file_name", out.dir = path.files){
  
  for(a in unique(md[,all_of(folder1)])){
    if(!exists(set_path(out.dir, a))){dir.create(set_path(out.dir, a))}
    
    for(b in unique(md[,all_of(folder2)])){
      if(!exists(set_path(out.dir, a, b))){dir.create(set_path(out.dir, a, b))}
      
      ffs <- md[md[all_of(folder1)] == a & md[all_of(folder2)] == b, all_of(file_name_col)]
      write.flowSet(fs[ffs,], outdir = set_path(out.dir, a, b))
    }
  }
}






cycombine_cluster <- function(df, features, xdim = 10, ydim = 10, maxK = 20, seed = 1234, scale = TRUE, plot.tree = TRUE){
  asdf <- t(df[,all_of(features)])
  
  if(scale){
    asdf <- scale_exprs(asdf, margin = 1, q = 0.01)
  }
  
  fsom <- ReadInput(flowFrame(t(asdf)))
  
  xdim <- xdim
  ydim <- ydim
  set.seed(seed = seed)
  som <- BuildSOM(fsom = fsom, xdim = xdim, ydim = ydim)
  maxK <- maxK
  mc <- suppressWarnings(suppressMessages(
    ConsensusClusterPlus(t(som$map$codes), 
                         maxK = maxK, reps = 100, 
                         distance = "euclidean", plot = "pdf", seed = seed)))
  
  codes <- data.frame(seq_len(xdim*ydim), map(mc[-1], "consensusClass"))
  codes <- mutate_all(codes, function(u) factor(u, levels = sort(unique(u))))
  colnames(codes) <- c(sprintf("som%s", xdim*ydim), sprintf("meta%s", seq_len(maxK)[-1]))
  
  flowsom <- merge(data.frame("cluster.id" = som$map$mapping[,1], "order" = 1:length(som$map$mapping[,1])), 
                   data.frame(seq_len(xdim*ydim), map(mc[-1], "consensusClass"))[c(1,maxK)], by.x = "cluster.id", by.y = "seq_len.xdim...ydim.")
  flowsom <- flowsom[order(flowsom$order),]
  
  df$labels <- factor(flowsom[,3], levels = c(1:maxK))
  
  if(plot.tree){
    clustree(codes, "meta") %>% print()
  }  
  return(df)
}


cycombine_run_DR <- function(df, cells = NULL, features, scale = TRUE) {
  
  # check validity of input arguments
  dr <- "UMAP"
  
  if (is.null(cells)) {
    # use all cells
    cs <- TRUE 
  } else {
    if (is.null(df$sample))
      stop("colData column sample_id not found,\n ", 
           " but is required to downsample cells.")
    stopifnot(
      is.numeric(cells), length(cells) == 1,
      as.integer(cells) == cells, cells > 0)
    # split cell indices by sample
    cs <- split(seq_len(nrow(df)), df$sample)
    # sample at most 'n' cells per sample
    cs <- unlist(lapply(cs, function(u)
      sample(u, min(cells, length(u)))))
  }
  
  # run dimension reduction
  asdf <- df %>% .[cs, all_of(features)] %>% t()
  
  if(scale){
    asdf <- scale_exprs(asdf, margin = 1, q = 0.01)
  }
  
  fun <- get(paste0("calculate", dr))
  embeddings <- fun(asdf)
  
  # return SCE when no cell subsetting has been done
  if (is.null(cells)){
    df$UMAP1 <- embeddings[,1]
    df$UMAP2 <- embeddings[,2]
    return(df)
  }else{  
    # else, match coordinates to subsetted cells from SCE
    m <- matrix(NA, nrow = nrow(df), ncol = ncol(embeddings))
    m[cs, ] <- embeddings
    df$UMAP1 <- m[,1]
    df$UMAP2 <- m[,2]
    return(df)
  }
}



cycombine_plot_clusters <- function(df, color_by = "labels", add_labels = TRUE){
  asdf <- df %>% dplyr::select("UMAP1", "UMAP2", all_of(color_by))
  
  if(add_labels){
    clustlabs <- asdf[c(all_of(color_by), "UMAP1", "UMAP2")]
    clustlabs <- clustlabs %>%
      group_by(across(all_of(color_by))) %>%
      summarise_all("median")
  }
  
  
  p1 <- ggplot(df, aes(x = UMAP1, y = UMAP2))+
    geom_point(aes(col = get(color_by)), size = min(4000/nrow(df), 3))+
    theme_minimal()+
    #scale_color_manual(values = brewer.paired(n = length(unique(df[,all_of(color_by)]))))+
    guides(col = guide_legend(override.aes = list(size = 3)))+
    labs(col = color_by, x = "UMAP1", y = "UMAP2")
  
  if(add_labels){
    p1 <- p1 + geom_text_repel(data = clustlabs, aes(x = UMAP1, y = UMAP2, label = get(color_by)), size = 5)
  }
  return(p1)
}


# Plot features in UMAP space
cycombine_plot_features <- function(df, features){
  df[,all_of(features)] <- df %>% dplyr::select(all_of(features)) %>% t() %>% scale_exprs(margin = 1) %>% t() %>% as.data.frame()
  
  ga <- df %>% gather(key = "key", value = "value", all_of(features))
  ga$key <- factor(ga$key, levels = all_of(features))
  featureplot <- ggplot(ga, aes(x = UMAP1, y = UMAP2))+
    geom_point(aes(col = value), size = min(4000/nrow(df), 3))+
    facet_wrap("key")+
    theme_minimal()+
    scale_color_gradientn(colors = hcl.colors(10, "Viridis"))+
    labs(col = "Scaled expression")
  featureplot
  return(featureplot)
}


correct_data <- function(df,
                         label,
                         covar = NULL,
                         anchor = NULL,
                         markers = NULL,
                         parametric = TRUE,
                         truncate = TRUE){
  message("Batch correcting data..")
  # Check for batch column
  cyCombine:::check_colname(colnames(df), "batch", "df")
  if (is.null(markers)){
    # Get markers
    markers <- df %>%
      cyCombine::get_markers()
  }
  
  # Add ID column to retain data order
  if(!"id" %in% colnames(df)) df$id <- 1:nrow(df)
  
  # Add label to df
  if(length(label) == 1){
    cyCombine:::check_colname(colnames(df), label, "df")
  }else{
    df$label <- label
    label <- "label"
  }
  
  # Add covar to df, if given
  if(!is.null(covar)){
    if(length(covar) == 1){
      cyCombine:::check_colname(colnames(df), covar, "df")
      df[[covar]] <- as.factor(df[[covar]])
    } else{
      # Covar was given as a vector
      df$covar <- as.factor(covar)
      covar <- "covar"
    }
    # Ensure there is more than 1 factor level
    if(nlevels(df[[covar]]) == 1) covar <- NULL
  }
  
  # Add anchor to df, if given
  if(!is.null(anchor)){
    if(length(anchor) == 1){
      cyCombine:::check_colname(colnames(df), anchor)
      df[[anchor]] <- as.factor(df[[anchor]])
    } else{
      # Anchor was given as a vector
      df$anchor <- as.factor(anchor)
      anchor <- "anchor"
    }
    # Ensure there is more than 1 factor level
    if(nlevels(df[[anchor]]) == 1) anchor <- NULL
  }
  
  corrected_data <- df %>%
    dplyr::group_by(.data[[label]]) %>%
    # Correct (modify) each label group with ComBat
    dplyr::group_modify(.keep = TRUE, function(df, ...){
      # Initiate anchor and covar counter
      num_covar <- 1
      num_anchor <- 1
      # Detect if only one batch is present in the node
      num_batches <- df$batch %>%
        factor() %>%
        nlevels()
      lab <- df[[label]][1] # Current label group
      if(num_batches == 1){
        batch <- df$batch[1]
        message(paste("Label group", lab, "only contains cells from batch", batch))
        df <- df %>% dplyr::select(-label) # Not removed from output, but removed here to prevent bug
        return(df)
      }
      message(paste("Correcting Label group", lab))
      # Calculate number of covars in the node
      if(!is.null(covar)){
        
        # Only use covar, if it does not confound with batch
        if(!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]]))){
          num_covar <- df[[covar]] %>%
            factor() %>%
            nlevels()
          
          # If a node is heavily skewed to a single covar, it should be treated as having only 1 covar.
          # Get number of cells in the condition with most cells
          covar_counts <- df %>%
            dplyr::count(.data[[covar]]) %>%
            dplyr::pull(n)
          
          if(sum(covar_counts) < max(covar_counts) + num_covar*5){
            message("The label group almost exclusively consists of cells from a single covar. Therefore, covar is ignored for this label group")
            num_covar <- 1
          }
        } else{
          message("Covar is confounded with batch. Ignoring covar in this label group")
        }
      }
      # Do a similar check on anchor
      if(!is.null(anchor)){
        if(!cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[anchor]]))){
          num_anchor <- df[[anchor]] %>%
            factor() %>%
            nlevels()
          
          # If a node is heavily skewed to a single anchor, it should be treated as having only 1 covar.
          # Get number of cells in the anchor with most cells
          anchor_counts <- df %>%
            dplyr::count(.data[[anchor]]) %>%
            dplyr::pull(n)
          
          if(sum(anchor_counts) < max(anchor_counts) + num_anchor*5){
            message("The label group almost exclusively consists of cells from a single anchor group. Therefore, anchor is ignored for this label group")
            num_anchor <- 1
          }
        } else{
          message("Anchor is confounded with batch. Ignoring anchor in this label group")
        }
      }
      if(num_covar > 1 & num_anchor > 1){
        # If neither covar nor anchor confounds with batch but they do each other, prioritise covar
        if(cyCombine:::check_confound(df$batch, stats::model.matrix(~df[[covar]] + df[[anchor]]))){
          num_anchor <- 1
          message("Anchor and covar are confounded. Ignoring anchor in this label group")
        }
      }
      # Compute ComBat correction
      ComBat_output <- df %>%
        dplyr::select(dplyr::all_of(markers)) %>%
        t() %>%
        # The as.character is to remove factor levels not present in the SOM node
        purrr::when(num_covar > 1 & num_anchor == 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~ df[[covar]]),
                                  par.prior = parametric),
                    num_covar > 1 & num_anchor > 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~df[[covar]] + df[[anchor]]),
                                  par.prior = parametric),
                    num_covar == 1 & num_anchor > 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  mod = stats::model.matrix(~df[[anchor]]),
                                  par.prior = parametric),
                    num_covar == 1 & num_anchor == 1 ~
                      sva::ComBat(.,
                                  batch = as.character(df$batch),
                                  par.prior = parametric)
        ) %>%
        t() %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(
          dplyr::select(df,
                        -dplyr::all_of(c(markers, label))))
        # Cap values to range of input data
        if(truncate){
          ComBat_output <- ComBat_output %>% dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                      function(x) {
                                        min <- min(df[[dplyr::cur_column()]])
                                        max <- max(df[[dplyr::cur_column()]])
                                        x <- ifelse(x < min, min, x)
                                        x <- ifelse(x > max, max, x)
                                        return(x)
                                      }))
        }

      # Only add covar column, if it is not null
      if(!is.null(covar)) ComBat_output[[covar]] <- df[[covar]]
      # Only add anchor column, if it is not null
      if(!is.null(anchor)) ComBat_output[[anchor]] <- df[[anchor]]
      
      return(ComBat_output)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(id) %>%
    dplyr::select(id, dplyr::everything()) %>%
    dplyr::mutate(batch = as.factor(batch))
  return(corrected_data)
}

Dataplot <- function(df = All,
                     type = "Percentages", #the type can be either: Percentages, Counts or Both
                     # group_by = c(Organ, Analysis, Timepoint),
                     x_plot = "timepoint",
                     y_plot = "Freq",
                     cluster_column = "cluster",
                     color_by,
                     rows = NULL, # Needs to be in the format "vars(Organ)"
                     cols = NULL, # Needs to be in the format "vars(Timepoint)"
                     text = 12,
                     title = 14,
                     ref = "H",
                     LOD= NULL) {
  pattern <- ifelse(type == "Percentages", "%", "Count")
  scale <- ifelse(type == "Percentages", percent, scientific)
  dftyped <- df #%>% dplyr::filter(grepl(pattern, Analysis))
  for(t in unique(dftyped$cluster_column)){
    int <- dftyped %>% dplyr::filter(cluster_column == t) #filter for measurements of the current type
    
    k <- int$cluster_column %>% unique() %>% length()
    
    # With the stat.test table we calculate the significance values, which we can later add to the plot.
    # Watch out that the T-test can only be performed if there's n>1 measurement per condition.
    # stat.test <- int %>%
    #   group_by(!!!rlang::syms(grouping_vars)) %>%
    #   t_test(Measurement ~ Timepoint, ref.group = ref) %>%
    #   #adjust_pvalue(method = "bonferroni") %>%
    #   add_significance("p")
    # stat.test
    # stat.test <- stat.test %>%
    #   add_xy_position(x = x_plot, dodge = 0.75, fun = "max", scales = "free")%>%
    #   dplyr::filter(p != "NaN")
    # 
    
    p1 <-ggplot(int, aes(x = get(x_plot), y = get(y_plot), interaction(get(x_plot), get(color_by))))+
      geom_boxplot(aes(color = get(color_by), fill = get(color_by)), alpha = 0.5, width = 0.75)+
      geom_point(aes(color = get(color_by), fill = get(color_by)), 
                 position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
                 alpha = 0.5, size = 2) +
      scale_y_continuous(labels = scale, expand = expansion(mult = c(.05, .3)))+ # we expand the y-axis by 30% to leave space for significance
      facet_grid(rows = rows, cols = cols, scales="free_y", switch = "both")+ #scales function is to account for groups that don't have a timepoint
      theme_bw()+
      stat_pvalue_manual(stat.test,  label = "p.signif", tip.length = 0.02, hide.ns = TRUE, remove.bracket = TRUE)+
      theme(panel.grid = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = text),
            axis.text.y = element_text(size = text),
            strip.text = element_text(size = title, face = "bold"),
            legend.title = element_text(size = title, face = "bold"),
            legend.text = element_text(size = text))+
      labs(y = ifelse(type == "Percentages", paste0("% of ", t), "Count"))+
      labs(color = color_by, fill = color_by)+
      scale_fill_manual(values = colvec)+
      scale_color_manual(values = colvec)
    if (!is.null(LOD)) {
      p1 <- p1 + geom_hline(yintercept = LOD, linetype = "dashed")
    }
    print(p1) #print plot (has to be specified in a for-loop)
  } 
}
