#' NEW_collect_cluster_bp
#' 
#' Collect breakpoints from a cluster of cells; this function extracts all breakpoints 
#' from a specified cluster of cells based on the results of the re-clustering step.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Clustering_output A data frame recording the pqArm clustering,
#'  and re-clustering results for each cell. This table tracks the clustering 
#'  history at each step.
#' @param Recluster_label An integer specifying the cluster from the re-clustering 
#'  step for which breakpoints should be extracted.
#'
#' @return A data frame containing breakpoint sites with the following columns:
#' 
#'   - `seqnames`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the breakpoint.
#'   - `end`: End position of the breakpoint.
#'   - `width`: Width of the fixed-bin size in the genomic regions.
#'   - `strand`: Strand information (`+` or `-`).
#'   - `copy.number`: Left copy number value at the breakpoint.
#'   - `cellID`: Identifier of the corresponding cell.
#'
NEW_collect_cluster_bp <- function(input, Clustering_output, Recluster_label)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 5.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  selected_files <- Clustering_output %>%
                    dplyr::filter(.data$Recluster_cluster %in% Recluster_label) %>%
                    dplyr::pull(.data$cellID)
  
  # Use lapply to gather breakpoints for all selected files at once
  breakpoints_list <- lapply(selected_files, function(i) {
    input[[i]]$breakpoints %>% as.data.frame() %>% dplyr::mutate(cellID = i)
  })
  
  # Combine all the results using bind_rows, which is more efficient than rbind in a loop
  breakpoints <- dplyr::bind_rows(breakpoints_list)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(breakpoints)
}