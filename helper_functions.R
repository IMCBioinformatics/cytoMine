### Created by Hena R. Ramay 
### IMC Bioinformatics Core
### This file includes custom written functions and modified ones from Cytofkit and from the bioconductor CyTOF workflow 



## numeric0
#' Title numeric0
#'
#' @param x 
#'
#' @return Logical value depending on if the scalar is numeric0 or not
#'
#' @examples 
#' v= Null
#' is.numeric0(v)
is.numeric0 <- function(x) {
  return(identical(x, numeric(0)))
}



## generic function for randomizing 0s to between -1 and -0.1
#' Title
#'
#' @param x depends on the type of object to be randomized
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
randomize_zeros<-function(x,...){
  UseMethod("randomize_zeros")
}

##
#' Title  randomize for expression matrix
#'
#' @param x is of type matrix
#' @param ... 
#'
#' @return matrix where 0s are randomized
#' @export
#'
#' @examples
randomize_zeros.matrix<-function(x,...)
{ apply(x,2,function(y) {
    rand <- runif(10000, -1, -0.1);
    idx<-which(y == 0);
    y[idx]<-sample(rand,length(idx),replace = TRUE);
    return(y)})-> exprs
  return(exprs)
}


## randomize for flowFrames
#' Title
#'
#' @param x an object of type flowFrame
#' @param ... 
#'
#' @return an object of type flowFrame
#' @export
#'
#' @examples
randomize_zeros.flowFrame<-function(x,...)
{
  
  apply(x@exprs,2,function(y) {
    rand <- runif(10000, -1, -0.1);
    idx<-which(y == 0);

    y[idx]<-sample(rand,length(idx),replace = TRUE);

    return(y)})-> x@exprs
  return(x)
}


## 
#' Randomize for flowSet
#'
#' @param x object of type flowSet
#' @param ... 
#'
#' @return obejct of type flowset
#' @export
#'
#' @examples
randomize_zeros.flowSet<-function(x,...)
{
  message("CytoMine: Randomizing 0s...")
fsApply(x,function(y){
  randomize_zeros(y)
}) 
}

#################################CODE##########################
### FILE read in options
## Individual files or flowSets can be read in.
## If concatenation option is on the files in the sample_info folder
## will be concatenated and returned
###Read FCS files
#' readFCS
#'
#' @param files, reads in FCS files if concat=TRUE concatenates them
#'
#' @return an object of type flowFrame, flowSet
#' @export
#'
#' @examples
readFCS <- function(files)
{
  if (length(files) < 2)
  { message("CytoMine: Reading in a single FlowFrame....")
    invisible(capture.output(fcs <- (read.FCS(files,
                                              transformation = FALSE,
                                              truncate_max_range = FALSE))))
  }
  else if ((length(files) >= 2) & (mode == "concat"))
  {
    message("CytoMine: Concataneting Files....")
    fcs <- concatFCS(files)
    message("CytoMine: Writing concatenated fcs file to ",outputDir)
    write.FCS(fcs,filename=paste0(outputDir,"/concatenated.fcs"))
    quit(status=0)
  }
  else if ((length(files) >= 2) & (concat == FALSE))
  {
    message("CytoMine: Reading flowset....")
    invisible(capture.output(fcs <- (read.flowSet(files,
                                                  transformation = FALSE,
                                                  truncate_max_range = FALSE))))
  } else {
    stop("CytoMine: No FCS file is found!")
  }
  return(fcs)
}



## Normalize markers to beads.
## Once cal also specify one file in the Flowset to use as baseline
## 
#' Normalize to beads generic function
#'
#' @param x object can be of type flowFrame or flowSet
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
normalize<-function(x,...){
  UseMethod("normalize")
}
## Normalize files for flowFrame
normalize.flowFrame<-function(fcs,beads,remove_beads=TRUE,outpath=NULL,
                              k=300,trim=5,norm_to=NULL)
{ message("CytoMine: Normalizing flowFrame to beads....")
  # The first frame is returned here as that contains the values of the markers without the beads
  temp<-normCytof(x=fcs, y=beads, k=k,out_path=outpath,
                   norm_to=norm_to)
  return(temp@frames[[ls(temp@frames)[1]]])
}

normalize.flowSet<-function(fcs,beads,remove_beads=TRUE,outpath=NULL,
                            k=300,trim=5,norm_to=NULL)
{
  message("CytoMine: Normalizing the FlowSet to beads....")
  fcs_normalize<-fsApply(fcs,function(x){
    normalize(x,beads,remove_beads=TRUE,k=k,outpath=outpath,
              norm_to=norm_to)
  })
}

### Debarcoding is done using barcodescheme file on a single FCS file
#debarcode
#' Title
#'
#' @param flowFrame to debarcode
#' @param bscheme barcode scheme is a dataframe with a specified format. Look at an example data files
#' @param out_path  output directory where to save the debarcoded files
#'
#' @return object re
#' @export
#'
#' @examples
debarcode<-function(fcs,bscheme,out_path)
{
  message("CytoMine:Debarcoding channels...")
  re <- assignPrelim(x=fcs, y=bscheme, verbose=TRUE)
  re <- estCutoffs(x=re, verbose=TRUE)
  re <- applyCutoffs(x=re)

  dedir<-outFCS(x=re,out_path = out_path)
  return(re)
}

## 
#' Density and count plots
#'
#' @param ex expression matris
#' @param ids marker ids
#' @param metadata data from sample_info file
#' @param outputDir output directory for saving plots
#' @param colors colors to use for the different samples
#'
#' @return nothing
#' @export
#'
#' @examples
basicPlots<-function(ex,ids,metadata,outputDir,colors,device)
{
  ggdf <- data.frame(sample_id = ids, ex)
  ggdf <- melt(ggdf, id.var = "sample_id",
               value.name = "expression", variable.name = "antigen")
  mm <- match(ggdf$sample_id, sample_id)
  ggdf$condition <- factor(md$Group[mm])
  
  ggsave(
    filename = paste0(outputDir,"/density_plots.",device),
    ggplot(ggdf, aes(x = expression, color = condition,group = sample_id)) +
      geom_density() +
      facet_wrap(~ antigen, nrow = 4, scales = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text = element_text(size = 7), 
            axis.text = element_text(size = 5)) +
      guides(color = guide_legend(ncol = 1)) +
      scale_color_manual(values = colors),
    width = w,height = 5)
  
  cell_table <- table(ids)
  ggct <- data.frame(sample_id = names(cell_table),
                     cell_counts = as.numeric(cell_table))
  mm <- match(ggct$sample_id,sample_id)
  ggct$condition <- md$Group[mm]
  
  ggsave(
    filename = paste0(outputDir,"/count_plots.",device),ggplot(ggct, aes(x = sample_id, y = cell_counts, fill = condition)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_x_discrete(drop = FALSE),width = w,height = 5)
}


#' Heatmap function
#'
#' @param expr matrix of marker values
#' @param expr01 scaled matrix of values between 0 and 1
#' @param cell_clustering results form Cell clustering
#' @param color_clusters  colors to assign to clusters
#' @param cluster_merging cluster merging strategy if the user is reassigning clusters
#' @param filename output filenames
#'
#' @return
#' @export
#'
#' @examples
plot_clustering_heatmap_wrapper <- function(expr, expr01,
                                            cell_clustering, color_clusters, 
                                            cluster_merging = NULL,filename){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(subset(expr_median,select=-c(cell_clustering)), method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(subset(expr01_median,select=-c(cell_clustering)))
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 10, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend,filename = filename)
  
}


# library(ggridges)
# 
# plot_clustering_distr_wrapper <- function(expr, cell_clustering){
#   # Calculate the median expression
#   cell_clustering <- factor(cell_clustering)
#   expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
#     group_by(cell_clustering) %>% summarize_all(funs(median))
#   print(colnames(expr_median))
#   # Sort the cell clusters with hierarchical clustering
#   d <- dist(expr_median, method = "euclidean")
#   
#   cluster_rows <- hclust(d, method = "average")
#   # Calculate cluster frequencies
#   freq_clust <- table(cell_clustering)
#   freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
#   cell_clustering <- factor(cell_clustering,
#                             labels = paste0(levels(cell_clustering), "  (", freq_clust, "%)"))
#  
#   ### Data organized per cluster
#   ggd <- melt(data.frame(cluster = cell_clustering, expr),
#               id.vars = "cluster", value.name = "expression",
#               variable.name = "antigen")
#   ggd$antigen <- factor(ggd$antigen)
#   print(ggd$antigen)
#   ggd$reference <- "no"
#   ### The reference data
#   ggd_bg <- ggd
#   ggd_bg$cluster <- "reference"
#   ggd_bg$reference <- "yes"
#  
#   ggd_plot <- rbind(ggd, ggd_bg)
#   ggd_plot$cluster <- factor(ggd_plot$cluster,
#                              levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
#   
#   ggplot() +
#     geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
#                                              color = reference, fill = reference), alpha = 0.3) +
#     facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
#     theme_ridges() +
#     theme(axis.text = element_text(size = 7),  
#           strip.text = element_text(size = 7), legend.position = "none")
#   
# }
# 




#' CytofAsinh taken from cytofkit 
#'
#' @param value 
#' @param cofactor 
#'
#' @return transform data to cytofAsinh
#' @export
#'
#' @examples
cytofAsinh <- function(value, cofactor = 5) {
  value <- value-1
  loID <- which(value < 0)
  if(length(loID) > 0)
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  value <- value / cofactor
  value <- asinh(value) 
  return(value)
}




### taken form CytofKit cytof_exprsExtract
## Removed compensation option
## flowframe
#' transformation 
#'
#' @param fcsFile 
#' @param verbose 
#' @param transformMethod 
#' @param scaleTo 
#' @param q 
#' @param l_w 
#' @param l_t 
#' @param l_m 
#' @param l_a 
#' @param a_a 
#' @param a_b 
#' @param a_c 
#' @param marker_id 
#'
#' @return
#' @export
#'
#' @examples
transform<-function (fcsFile, verbose = FALSE, #comp = FALSE, 
                     transformMethod = c("cytofAsinh", "logicle", "arcsinh", "none"), 
                     scaleTo = NULL, 
                     q = 0.05, l_w = 0.1, l_t = 4000, l_m = 4.5, l_a = 0, a_a = 1, 
                     a_b = 1, a_c = 0,marker_id=NULL) 
{
  transformMethod <- match.arg(transformMethod)
  fcs<-fcsFile
  #name <- sub(".fcs$", "", basename(fcsFile))
  # if (verbose) {
  #   fcs <- read.FCS(fcsFile, transformation = FALSE)
  # }
  # else {
  #   fcs <- suppressWarnings(read.FCS(fcsFile, transformation = FALSE))
  # }
  # if (is.matrix(comp) || is.data.frame(comp)) {
  #   fcs <- applyComp(fcs, comp)
  #   cat("    Compensation is applied on", fcsFile, "\\n")
  # }
  # else if (isTRUE(comp)) {
  #   if (!is.null(fcs@description$SPILL)) {
  #     fcs <- applyComp(fcs, fcs@description[["SPILL"]])
  #     cat("    Compensation is applied on ", fcsFile, "\\n")
  #   }
  #   else if (!is.null(fcs@description$SPILLOVER)) {
  #     fcs <- applyComp(fcs, fcs@description[["SPILLOVER"]])
  #     cat("    Compensation is applied on ", fcsFile, "\\n")
  #   }
  #   else if (!is.null(fcs@description$COMP)) {
  #     fcs <- applyComp(fcs, fcs@description[["COMP"]])
  #     cat("    Compensation is applied on ", fcsFile, "\\n")
  #   }
  #   else {
  #     warning("Cannot find compensation matrix in the FCS files!\\n                    Please CHECK the keyword of 'SPILL', 'SPILLOVER', or 'COMP'\\n                    in the FCS file and make sure it stores the compensation matrix.")
  #   }
  # }
  #### find parameters to transform- avoid Time and FSC
  pd <- fcs@parameters@data
  exclude_channels <- grep("Time|Event", colnames(fcs@exprs), 
                           ignore.case = TRUE)
  if(is.null(marker_id)){
  marker_id <- setdiff(seq_along(colnames(fcs@exprs)), exclude_channels)
  } 
  size_channels <- grep("FSC|SSC", colnames(fcs@exprs), ignore.case = TRUE)
  transMarker_id <- setdiff(marker_id, size_channels)
  ### Tranformations
  switch(transformMethod, cytofAsinh = {
    data <- fcs@exprs
    data[, transMarker_id] <- apply(data[, transMarker_id, 
                                         drop = FALSE], 2, cytofAsinh)
    exprs <- data[, marker_id, drop = FALSE]
  }, logicle = {
    #print("inlogicle")
    data <- fcs@exprs
    trans <- flowCore::logicleTransform(w = l_w, t = l_t, 
                                        m = l_m, a = l_a)
    data[, transMarker_id] <- apply(data[, transMarker_id, 
                                         drop = FALSE], 2, trans)
    exprs <- data[, marker_id, drop = FALSE]
  }, arcsinh = {
    data <- fcs@exprs
    trans <- flowCore::arcsinhTransform(a = a_a, b = a_b, 
                                        c = a_c)
    data[, transMarker_id] <- apply(data[, transMarker_id, 
                                         drop = FALSE], 2, trans)
    exprs <- data[, marker_id, drop = FALSE]
  }, none = {
    data <- fcs@exprs
    exprs <- data[, marker_id, drop = FALSE]
  })
  if (length(size_channels) > 0) {
    if (any(size_channels %in% marker_id)) {
      used_size_channel <- size_channels[size_channels %in% 
                                           marker_id]
      used_size_channel_id <- match(used_size_channel, 
                                    marker_id)
      exprs[, used_size_channel_id] <- apply(exprs[, used_size_channel_id, 
                                                   drop = FALSE], 2, function(x) scaleData(x, range = c(0, 
                                                                                                        4.5)))
    }
  }
  if (!is.null(scaleTo)) {
    exprs <- apply(exprs, 2, function(x) scaleData(x, scaleTo))
  }
  
  if(all(is.na(pd$desc)))
  {
    col_names <-pd$name
  }
  else
  {
  #col_names <- paste0(pd$name, "<", pd$desc, ">")
    col_names <- pd$desc
  }
  colnames(exprs) <- col_names[marker_id]
  # row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
  return(exprs)
}



#### ### plot_markers
## generic function
plot_marker_exp<-function(x,...){
  UseMethod("plot_marker_exp")
}


plot_marker_exp.flowFrame<-function(fcs,markers=NULL,prefix=NULL,...)
{
  colnames(fcs@exprs)<-fcs@parameters@data$desc
  
  if(is.null(markers))
  {
    pd <- fcs@parameters@data
    exclude_channels <- grep("Time|Event", colnames(fcs@exprs), 
                             ignore.case = TRUE)
    marker_id <- setdiff(seq_along(colnames(fcs@exprs)), exclude_channels)
    size_channels <- grep("FSC|SSC|Cell_length", colnames(fcs@exprs), ignore.case = TRUE)
    marker_id <- setdiff(marker_id, size_channels)
  }
  else{
    marker_id=which(markers %in% colnames(fcs@exprs))
    print(paste0("markers not found in fcs file: " ,setdiff(markers,colnames(fcs@exprs))))
  }
  
  fcs_selected<-melt(fcs@exprs[, marker_id, drop = FALSE])
  print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  print(colnames(fcs_selected))
  
  p<-ggplot(fcs_selected,aes(x=fcs_selected$value))+
    geom_density() +
    facet_wrap(~ Var2, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5))+
    xlab("expression")
  ggsave(filename = paste0(outputDir,"/plot_expression_markers.",device),plot = p)
}


#####
plot_marker_exp.matrix<-function(exprs,markers=NULL,prefix=NULL,...)
{ print("--------------------------in plot marker matrix-------------")
  #colnames(exprs)<-strsplit2(strsplit2(colnames(exprs),split = "<")[,2],split = ">")
  # print("******")
  # print(markers)
  # print("******")
  # print(colnames(exprs))
  # print("******")
  print(markers)
  print("--------------------------in colnames exprs-------------")
  print(colnames(exprs))
  
  if(is.null(markers))
  {
  marker_id=colnames(exprs)
  #print("in markers")
  }
  else{
    marker_id=which(markers %in% colnames(exprs))
    
    warning(paste0("markers not found in fcs file: " ,setdiff(colnames(exprs)[marker_id],colnames(exprs))))
    # print("in markers else")
  }
  # print(marker_id)
  # print(colnames(exprs))
  
  fcs_selected<-melt(exprs[, marker_id, drop = FALSE])
  
   print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
   print(colnames(fcs_selected))
   #print(fcs_selected)
  
  p<-ggplot(fcs_selected,aes(x=value))+
    geom_density() +
    facet_wrap(~ Var2, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5))+
    xlab("expression")
  ggsave(filename = paste0(outputDir,"/plot_expression_markers_",prefix,".png"),plot = p,device = "png")
}





plot_tsne_marker<-function(exprs,tsne_dims,markers=NULL,prefix=NULL,...)
{
  #print(tsne_dims)
  colnames(tsne_dims)<-c("tSNE1","tSNE2")
  #colnames(exprs)<-strsplit2(strsplit2(colnames(exprs),split = "<")[,2],split = ">")
  #transform values between 0 to 1
  rng <- colQuantiles(as.matrix(exprs), probs = c(0.01, 0.99))
  expr01 <- t((t(as.matrix(exprs)) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1
  exprs<-expr01
  
  if(is.null(markers))
  {
    marker_id=colnames(exprs)
    print("in markers")
  }
  else{
    marker_id=which(markers %in% colnames(exprs))
    
    print(c("markerid",marker_id))
    
    #warning(paste0("markers not found in fcs file: " ,setdiff(colnames(exprs)[marker_id],colnames(exprs))))
    # print("in markers else")
  }
  # print(marker_id)
  # print(colnames(exprs))
  
  
  print("gothere")
  mat<-cbind(as.data.frame(exprs[, marker_id, drop = FALSE]),tsne_dims)
  fcs_selected<-melt(mat,id.vars=colnames(tsne_dims))
  # print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
  # print(colnames(fcs_selected))
  # print(fcs_selected)
  print(fcs_selected)  

  p<-ggplot(fcs_selected,aes(x=tSNE1,y=tSNE2))+
    geom_point(aes(colour =value),size=0.05) + 
    facet_wrap(~ variable)+
    scale_colour_gradientn(colours=rev(rainbow(4)))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5))+
    xlab("expression")
  ggsave(filename = paste0(outputDir,"/plot_tsne_expression_markers_",prefix,".png"),plot = p,device = "png")
}




####################################

#### GROUP PLOTS
group_downsample<-function(dframe,downsample)
{
  ## Data subsampling: create indices by sample
  ind <- split(1:length(dframe$condition), dframe$condition)
  ## How many cells to downsample per-sample
  ncells <- pmin(table(dframe$condition), downsample)
  ## Get subsampled indices
  g_inds <- unlist(lapply(names(ind), function(i){
    s <- sample(ind[[i]], ncells[i], replace = FALSE)
  }))
  return(dframe[g_inds,])
}
