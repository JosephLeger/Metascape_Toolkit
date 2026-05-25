#===============================================================================
# METASCAPE --------------------------------------------------------------------
#===============================================================================

library(rlang)
library(stringr)
library(readxl)
library(forcats)
library(Hmisc)
library(pheatmap)

meta_Read <- function(file){
  
  # Generate a metaObject composed of lists corresponding to metascape results
  #
  # file = path of .zip metascape output
  
  suppressPackageStartupMessages(library(readxl))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(Hmisc))
  
  FINAL_GO <- read.csv(unzip(file, 'Enrichment_GO/_FINAL_GO.csv'), header = T, 
                       check.names = F)
  
  object    <- list()
  n_process <- max(FINAL_GO$GROUP_ID)
  
  # Iterate over every annotated process
  for(id in 1:n_process){
    x           <- FINAL_GO[FINAL_GO$GROUP_ID %in% id,]
    new_process <- list()
    
    # Iterate over every subprocess belonging to current main process
    for(i in 1:nrow(x)){
      new_process[[x$Description[i]]] <- list(
        Description           = x$Description[i],
        GO                    = x$GO[i],
        Category              = x$Category[i],
        Parent_GO             = x$PARENT_GO[i],
        Enrichment            = x$Enrichment[i],
        LogP                  = x$LogP[i],
        LogQ                  = x$`Log(q-value)`[i],
        `Z-score`             = x$`Z-score`[i],
        Hits                  = unlist(strsplit(x$Hits[i], '|', fixed = T)),
        GeneID                = unlist(strsplit(x$GeneID[i], '|', fixed = T)),
        `#TotalGeneInLibrary` = x$`#TotalGeneInLibrary`[i],
        `#GeneInGo`           = x$`#GeneInGO`[i],
        `#GeneInHitList`      = x$`#GeneInHitList`[i],
        `#GeneInGOAndHitList` = x$`#GeneInGOAndHitList`[i],
        `%InGO`               = x$`%InGO`[i],
        `STDV_%InGO`          = x$`STDV %InGO`[i]
      )
      
      # Add current process and subprocesses
      object[[capitalize(x$Description[1])]] <- new_process
    }
  }
  unlink('Enrichment_GO', recursive = T)
  return(object)
}

#-------------------------------------------------------------------------------

meta_Subset <- function(x, feature='LogQ', thr=NA, top=NA){
  
  # Subset a meta_object based on a define threshold or getting top n paths
  #
  # x       = meta object obtained using meta_Read() function
  # feature = name of feature for subset (LogStat will be considered -LogStat)
  # thr     = minimum feature value required (mutually exclusive with 'top')
  # top     = number of top paths to keep (mutually exclusive with 'thr')
  
  # Look for error in specifications
  if(is.na(thr) & is.na(top)){
    warning('Neither \'thr\' or \'top\' specified. Returning original object.')
    return(x)
  }else if(!is.na(thr) & !is.na(top)){
    warning('Both \'thr\' and \'top\' specified. Only \'top\' will be applied.')
    thr=NA
  }
  
  # test
  y <- data.frame()
  for(i in 1:length(x)){
    # Get main corresponding pathways name and feature
    y <- rbind(y, c(names(x)[i], x[[i]][[1]][[feature]]))
  }
  colnames(y) <- c('Path', 'Feature')
  
  # If specified feature is stat (LogP or LogQ), transform it into -LogStat 
  if(feature %in% c('LogP', 'LogQ')){
    y$Feature <- -as.numeric(y$Feature)
  }
  
  # Sort resulting table rows
  y <- y[rev(order(as.numeric(y$Feature))),]
  
  # If thr=NA, apply top to generate x_sub. Otherwise apply thr
  if(is.na(thr) & top < nrow(y)){
    y <- y[c(1:top),]
    x_sub <- x[y$Path]
  }else if(is.na(top)){
    y <- y[y$Feature > as.numeric(thr),]
    x_sub <- x[y$Path]
  }else{
    x_sub <- x
  }
  
  return(x_sub)
}

#-------------------------------------------------------------------------------

meta_GetGenes <- function(meta_object, pathlist, intersect = T){
  
  # Get list of genes belonging to provided pathways
  #
  # meta_object = meta object generated using meta_Read() function
  # patlist     = list containing metascape pathways sorted in categories 
  
  # Initiate final list
  result <- list()
  # For each category of pathways defined
  for(n in names(pathlist)){
    category     <- list(Merged = c())
    all_pathways <- pathlist[[n]]
    # For each pathway in given category
    for(p in all_pathways){
      # If pathway not found in metascape output, print warning
      if(p %!in% names(meta_object)){
        warning(paste('Pathway not found', n))
      }else{
        # Other wise add gene list by pathway and append merged gene list
        category[[p]]        <- meta_object[[p]][[1]][['Hits']]
        category[['Merged']] <- unique(
          c(category[['Merged']], meta_object[[p]][[1]][['Hits']]))
      }
    }
    result[[n]] <- category
  }

  if(intersect){
    # Define intersect by category
    intersect <- list(Merged = c())
    for(i in 1:(length(result)-1)){
      n1      <- names(result)[i]
      for(j in (i+1):length(result)){
        n2    <- names(result)[j]
        genes <- result[[i]][['Merged']][
          result[[i]][['Merged']] %in% result[[j]][['Merged']]]
        intersect[[paste0('[', n1, ']-[', n2, ']')]] <- genes
        intersect[['Merged']] <- unique(c(intersect[['Merged']], genes))
      }
    }
    result[['Intersect']] <- intersect
  }
  
  return(result)
}

#-------------------------------------------------------------------------------

meta_IdentifyGroups <- function(file){
  
  # Return a list with necessary info to identify groups needed for comparisons
  #
  # file = path to Metascape output .zip file as [group1_vs_group2_UP_*.zip]
  
  suppressPackageStartupMessages(library(stringr))
  
  name       <- str_remove_all(file, '.*/')
  name.split <- unlist(str_split(name, '_'))
  group_id   <- list(upin = ifelse('UP' %in% name.split, 
                                   name.split[1], 
                                   name.split[3]),
                     groups = name.split[c(1,3)],
                     dirname = paste0(name.split[1], '_vs_', name.split[3]))
  return(group_id)
}

#-------------------------------------------------------------------------------

meta_DrawEnrichment <- function(x, 
                                stat = 'Q',
                                show.max = length(x),
                                order.by = 'Enrichment',
                                enrich.cutoff = 0,
                                stat.cutoff = 0,
                                count.cutoff = 0,
                                GOterm  = F,
                                removeSpecies = '',
                                title = 'GO Enrichment'){
  
  # Draw Enrichment Factor and statistics of processes contained in a metaObject
  #
  # x              = metaObject generated by meta_Read()
  # stat           = define metric to use, must be 'P' or 'Q'
  # show.max       = define the number of top stat pathways to show 
  # order.by       = define how to order (Enrichment, Count or Stat)
  # enrich.cutoff  = define a minimal enrichment value to consider
  # stat.cutoff    = define a minimal stat value to consider
  # count.cutoff   = define a minimal gene number value to consider
  # GOterm         = whether show GO Term ID in front of process description
  # remove.species = trim process names ex: ' - Mus musculus \\(house mouse\\)'
  # title          = main title of the plot
  
  suppressPackageStartupMessages(library(forcats))
  suppressPackageStartupMessages(library(Hmisc))
  suppressPackageStartupMessages(library(ggplot2))
  
  # Identify metric chosen as stat
  if(tolower(stat) %in% c('p', 'logp')){
    stat        <- 'LogP'
    stat.legend <- '-Log10(P-value)'
  }else if(tolower(stat) %in% c('q', 'logq')){
    stat        <- 'LogQ'
    stat.legend <- '-Log10(q-value)'
  }else{
    warning('stat parameter must be \'P\' or \'Q\'')
  }
  
  # Build a dataframe with required info for metaPlot
  y      <- data.frame()
  for(i in 1:length(x)){
    
    # Define process name depending on selected options
    name <- ifelse(GOterm, 
                   paste0(x[[i]][[1]]$GO, ': ', x[[i]][[1]]$Description),
                   names(x)[i])
    name <- ifelse(removeSpecies %!in% '',
                   unlist(str_split(name, removeSpecies))[1],
                   name)
    
    y <- rbind(y, c(
      name,
      x[[i]][[1]]$`#GeneInGOAndHitList`,
      x[[i]][[1]]$Enrichment,
      x[[i]][[1]]$LogP,
      x[[i]][[1]]$LogQ,
      x[[i]][[1]]$`Z-score`))
  }
  colnames(y) <- c('Process', 'Count', 'Enrichment', 'LogP', 'LogQ', 'Zscore')
  
  # Sort and subset based on defined cutoffs
  #y <- y[rev(order(as.numeric(y$Enrichment))),]
  y <- y[as.numeric(y$Count) >= count.cutoff,]
  y <- y[as.numeric(y$Enrichment) >= enrich.cutoff,]
  y <- y[-as.numeric(y[,stat]) >= stat.cutoff,]
  # Create stat column based on user choice
  y$stat <- y[,colnames(y) %in% stat]
  
  # Keep only remaining top stat pathways if needed
  if(nrow(y) > show.max){
    y <- y[order(as.numeric(y$stat)),]
    y <- y[c(1:show.max),]
  }
  
  # Define the way to order final plot based on order.by parameter
  if(tolower(order.by) %in% c('enrichment')){
    y$order <- as.numeric(y$Enrichment)
  }else if(tolower(order.by) %in% c('count')){
    y$order <- as.numeric(y$Count)
  }else if(tolower(order.by) %in% c('stat')){
    y$order <- -as.numeric(y$stat)
  }
  
  # Draw plot
  y %>% 
    ggplot(aes(as.numeric(Enrichment), 
               fct_reorder(Process, as.numeric(order)))) + 
    geom_segment(aes(xend=0, yend = Process)) +
    geom_point(aes(color=-as.numeric(stat), size = as.numeric(Count))) +
    labs(color = stat.legend, size = 'Gene Number') +
    scale_color_viridis_c(guide=guide_colorbar(reverse=F)) +
    scale_size_continuous(range=c(1, 6)) +
    theme_minimal() + 
    xlab('Enrichment Factor') +
    ylab(NULL) + 
    ggtitle(title)
}

#-------------------------------------------------------------------------------

# Copied from : https://stackoverflow.com/questions/43051525/how-to-draw-
#               pheatmap-plot-to-screen-and-also-save-to-file

save_pheatmap_pdf <- function(x, filename, width=9, height=9) {
  
  # Save pheatmap into PDF file
  #
  # x        = variable containing results of pheatmap()
  # filename = pathway for output 
  
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#-------------------------------------------------------------------------------

meta_DrawHeatmaps <- function(x, 
                              counts, 
                              sampleSheet, 
                              groups = NA,
                              by.group = F,
                              logT = F,
                              cluster_cols = F,
                              cluster_rows = F,
                              scale = 'row',
                              color = colorRampPalette(
                                c('blue', 'white', 'red'))(100),
                              removeSpecies = '', 
                              outdir = NA, 
                              save = F){
  
  # Draw Heatmap of genes expression for each process contained in a metaObject
  #
  # x           = metaObject generated by meta_Read() 
  # counts      = table containing counts for each sample with 'Symbol' column
  # sampleSheet = table containing Sample and Group information
  # groups      = order of groups to consider (default = all groups unsorted)
  # by.group    = whether show mean expression by group rather than by sample
  
  suppressPackageStartupMessages(library(rlang))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(pheatmap))
  
  # Get info from sampleSheet and order samples - - - - - - - - - - - - - - - - 
  if(length(groups) == 1){
    groups  <- names(table(sampleSheet$Group))
  }
  
  # Calculate mean expression if required or just filter count table
  if(by.group){
    mean.counts <- list(Symbol = counts$Symbol)
    for(g in groups){
      samples          <- sampleSheet$Sample[sampleSheet$Group %in% g]
      sub.counts       <- counts[,colnames(counts) %in% samples]
      mean.counts[[g]] <- rowSums(sub.counts)/ncol(sub.counts)
    }
    counts <- data.frame(mean.counts)
    
  }else{
    samples   <- c()
    for(g in groups){
      samples <- c(samples, sampleSheet$Sample[sampleSheet$Group %in% g])
    }  
    counts <- counts[,c('Symbol', samples)]
  }
  
  # Set up for file management - - - - - - - - - - - - - - - - - - - - - - - - - 
  for(i in 1:length(x)){
    dirname <- ifelse(removeSpecies %!in% '', unlist(str_split(
      names(x)[i], removeSpecies))[1], names(x)[i])
    dirname <- str_replace_all(dirname, '[ ><+=():/]', '_') 
    dirname <- paste0(i, '_', dirname)
    
    # Check if dirname is able to contain all files considering path <= 256 char
    if(nchar(paste0(outdir, '/', dirname, '/999_.pdf')) >= 256){
      l_dirname <- 256 - nchar(paste0(outdir, '/./999_.pdf'))
      dirname   <- substr(dirname, 1, l_dirname)
    }
    
    # Draw all Heatmaps - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    for(j in 1:length(x[[i]])){
      title    <-  ifelse(removeSpecies %!in% '',
                          unlist(str_split(names(x[[i]])[j], removeSpecies))[1],
                          names(x[[i]])[j])   
      filename <- str_replace_all(title, '[ ><+=():/]', '_')
      genes    <- x[[i]][[j]]$Hits
      
      # Subset and filter expressed genes
      y        <- counts[counts$Symbol %in% genes,]
      y        <- y[rowSums(y[,colnames(y) %!in% 'Symbol']) != 0,]
      
      # Apply logT if required
      if(logT){
        y.val <- log(y[,colnames(y) %!in% 'Symbol']+1)
        y     <- cbind(Symbol = y$Symbol, y.val)
      }
      # Add Symbol column at the end of this process
      
      
      
      plot <- pheatmap(y[,colnames(y) %!in% 'Symbol'],
                       cluster_cols = cluster_cols, 
                       cluster_rows = (cluster_rows & nrow(y) > 2),
                       scale = scale, 
                       color = color,
                       labels_row = y$Symbol, main = title)
      
      if(save){
        file_path   <- paste0(outdir, '/', dirname, '/', j, '_', 
                              filename, '.pdf')
        # Adapt filename to path lentgh
        if(nchar(file_path) >= 256){
          file_path <- paste0(outdir, '/', dirname, '/', j, '.pdf')
        }
        # Create directory if not already created
        if(!file.exists(paste0(outdir, '/', dirname))){
          dir.create(file.path(outdir, dirname))
        }
        save_pheatmap_pdf(plot, file_path)
      }
    }
  }
}

#-------------------------------------------------------------------------------

meta_DrawPathmap <- function(sub_counts, 
                             genelist, 
                             group.colors=ColorBlind,
                             title = 'Pathmap'){
  
  # Draw a heatmap separating previously defined processes of interest
  #
  # sub_counts = subseted table containing only columnes and rows of interest
  # genelist   = list resulting from meta_GetGenes() function
  # colors     = colors 
  
  # Initialize empty objects
  org_table  <- data.frame()
  gene_group <- c() 
  gaps_id    <- c(0)
  
  # For each process, select genes that are sepcific and append org_table
  for(n in names(genelist)){
    if(n %!in% c('Others', 'Intersect')){
      lines      <- sub_counts[
        sub_counts$Symbol %in% genelist[[n]][['Merged']] &
          sub_counts$Symbol %!in% genelist[['Intersect']][['Merged']],]
      org_table  <- rbind(org_table,lines)
      gene_group <- c(gene_group, rep(n, nrow(lines)))
      gaps_id    <- c(gaps_id, nrow(lines)+gaps_id[length(gaps_id)])
    }
  }
  
  # Add non-categorized genes (Others, Intersect and not annotated bu sig)
  non_cat    <- sub_counts$Symbol[sub_counts$Symbol %!in% org_table$Symbol]
  lines      <- sub_counts[sub_counts$Symbol %in% non_cat,]
  org_table  <- rbind(org_table, lines)
  gene_group <- c(gene_group, rep('Others', nrow(lines)))
  gene_group <- factor(gene_group, level = c(
    names(genelist)[names(genelist) %!in% c('Other', 'Intersect')], 'Others'))
  gene_group <- as.data.frame(gene_group)
  rownames(gene_group) <- rownames(org_table)
  
  # Attribute colors to each process
  cols <- c()
  for(i in 1:length(names(table(gene_group)))){
    cols <- c(cols, group.colors[i])
  }
  names(cols) <- names(table(gene_group))
  ann_colors  <- list(gene_group=cols) 
  
  # Draw plot
  plot <- pheatmap(
    org_table[,-1],
    scale = 'row',
    cluster_cols = F,
    cluster_rows = F,
    gaps_row = gaps_id,
    gaps_col = c(8),
    color = colorRampPalette(c('blue', 'white', 'red'))(100),
    show_rownames = F,
    annotation_row = gene_group,
    annotation_colors = ann_colors,
    main = title)
  
  return(plot)
}

#-------------------------------------------------------------------------------
