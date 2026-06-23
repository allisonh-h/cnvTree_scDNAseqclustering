#' Based on AneuFinder::clusterHMMs, calcuate distance than hierrachical clustering
#'
#' This function computes the pairwise distance between cells based on their 
#' copy number variations and performs hierarchical clustering.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param selected A character vector specifying the `cellID`s of the cells to be 
#'  included in clustering.
#' @param exclude.regions A `GRanges` object specifying genomic regions to exclude 
#'  from clustering computation. This is useful for filtering out regions with 
#'  artifacts.
#'
#' @return A list containing:
#'   - `ordered_indices`: The ordered indices of cells based on hierarchical clustering.
#' @export
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' Example_data <- changeFormat(file = file_path, core = 4)
#' Clustering_Result <- clusterbyHMM(input = Example_data, 
#'                                   selected = names(Example_data)[1:10])
#' }
#'
clusterbyHMM <- function(input, selected, exclude.regions = NULL)
{
  fucStep <- paste0(" 2.0_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  message("Checking column 'copy.number'  ...")
  
  hmms <- input[selected]
  hmms2use <- numeric()
  for (i1 in 1:length(hmms)) {
    hmm <- hmms[[i1]]
    if (!is.null(hmm$bins$copy.number)) {
      if (is.null(hmm$ID)) {
        stop("Need ID to continue.")
      }
      hmms2use[hmm$ID] <- i1
    }
  }
  hmms <- hmms[hmms2use]
  hc <- NULL

  message("Making consensus template ...")
  
  if (!is.null(hmms[[1]]$bins$copy.number)) {
    constates <- sapply(hmms, function(hmm) {
      hmm$bins$copy.number
    })
  }
  constates[is.na(constates)] <- 0
  vars <- apply(constates, 1, stats::var, na.rm = TRUE)

  message("Clustering ...")
  if (!is.null(exclude.regions)) {
    ind <- GenomicRanges::findOverlaps(hmms[[1]]$bins, exclude.regions)@from
    constates <- constates[-ind, ]
  }

  message("Distance calculating...")
  Dist <- Rfast::Dist(t(constates), method = "euclidean")

  # dist <- parallelDist::parDist(t(constates),
  #                               method = "euclidean",
  #                               threads = 5) # threads

  Dist_as_dist <- stats::as.dist(Dist)
  message("hierarchical clustering...")
  hc <- stats::hclust(Dist_as_dist)

  # message("Reordering ...")
  hmms2use <- hmms2use[hc$order]

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  return(list(IDorder = hmms2use, hclust = hc))
}


#' Construct a phylogenetic tree from copy number variation data
#'
#' This function performs hierarchical clustering on selected cells based on 
#' their copy number variations and derives a phylogenetic tree structure.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param selected A character vector specifying the `cellID`s of the cells to be 
#'  included in the phylogenetic analysis.
#'
#' @return A data frame with two columns:
#'   - `cellID`: The unique identifier of each cell.
#'   - `cluster`: The assigned cluster label (k = 2).
#'
CutTree_final <- function(input, selected)
{
  fucStep <- paste0(" 2.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  # 分群的原始檔，後面要用他作為基底
  message("Clustering and Data processing ...")

  clust <- clusterbyHMM(input = input, selected = selected)
  Clust_cuttree <- data.frame(cluster = stats::cutree(clust[["hclust"]], k = 2),
                              cell = names(clust$IDorder))

  Clust_cuttree <- Clust_cuttree %>%
                   dplyr::mutate(cluster = as.numeric(.data$cluster)) %>%
                   dplyr::arrange(dplyr::desc(.data$cluster))

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Clust_cuttree)
}


#' Divide cells into two groups based on copy number variation in a specific cluster
#'
#' This function separates cells into two groups based on copy number variation 
#' in a specified cluster label from a provided clustering result table.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param Cluster_label An integer specifying the cluster label used to divide the 
#'  cells into two groups.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @return A data frame with two columns:
#'   - `cellID`: The unique identifier of each cell.
#'   - `cluster`: The updated cluster assignment (k = 2).
#'
CutTree <- function(input, Template, Cluster_label)
{
  fucStep <- paste0(" 2.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  # selected.files建立
  # Cluster_label: Clust_cuttree$cluster中的分群數字
  selected.files <- NULL
  selected.files <- subset(Template, .data$cluster %in% c(Cluster_label))$cell

  message("Divided the cluster in k=2 ")
  clust <- clusterbyHMM(input = input, selected = selected.files, exclude.regions = NULL)
  Clust_2 <- data.frame(NewCluster = stats::cutree(clust[["hclust"]], k = 2),
                        cell = names(clust$IDorder))


  Template <- merge(Template, Clust_2, by = "cell", all = TRUE)
  Template$NewCluster <- ifelse(is.na(Template$NewCluster) == T, 
                                0, Template$NewCluster)

  Clust_2 <- which(is.na(Template$NewCluster) == F)

  count = max(Template$cluster, na.rm = TRUE)

  Template$cluster <- dplyr::case_when(
    Template$NewCluster == 1 ~ count + 1,
    Template$NewCluster == 2 ~ count + 2,
    TRUE ~ Template$cluster
  )

  Template <- Template[, !colnames(Template) %in% "NewCluster"]

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Template)
}


#' Count the number of cells in a specific cluster
#'
#' This function calculates the total number of cells that belong to a specified 
#' cluster.
#'
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param Cluster_label An integer specifying the cluster for which the number of 
#'  cells will be counted.
#'
#' @return An integer representing the number of cells in the specified cluster.
#'
Cluster_num <- function(Template, Cluster_label)
{
  fucStep <- paste0(" 2.3_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  # cat("Calculating numbers of cell in Cluster", Cluster_label, "...\n")

  num <- data.frame(table(Template$cluster))
  num <- num[which(num$Var1 == Cluster_label), 2]

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(num)
}


#' Compute the similarity of cells within a cluster
#'
#' This function calculates the overall similarity of cells within a specified 
#' cluster using a precomputed cell-to-cell similarity matrix.
#'
#' @param Template A data frame containing two columns:
#'   - `cellID`: Unique identifier for each cell.
#'   - `cluster`: Cluster assignment for each cell.
#' @param SimCells A square matrix where each element `[i, j]` represents the 
#'  similarity score between `cellID[i]` and `cellID[j]`.
#' @param Cluster_label A character vector containing the `cellID`s of cells that 
#'  belong to the specified cluster.
#'
#' @return A numeric value representing the overall similarity of the cells within 
#'  the specified cluster.
#'
Cluster_sim <- function(Template, SimCells, Cluster_label)
{
  fucStep <- paste0(" 2.4_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  # cat("Calculating cell similarity in Cluster ",  Cluster_label, " ...\n")

  selected <- Template %>% 
              dplyr::filter(.data$cluster %in% Cluster_label) %>% 
              dplyr::pull(.data$cell)
  selected <- which(SimCells$IDorder %in% selected)

  Similarity <- SimCells$similarity[selected, selected]
  Similarity[is.na(Similarity)] <- 0
  Similarity <- mean(Similarity)

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Similarity)
}


#' Computes pairwise similarity between cells
#'
#' This function calculates the similarity values between cells based on their 
#' copy number variations at the bin level.
#'
#' @param binsMatrix A numeric matrix where each row represents a genomic bin
#'   and each column represents a cell. The values indicate copy number variations.
#'
#' @return A square numeric matrix where each element `[i, j]` represents the 
#'  similarity score between `cell[i]` and `cell[j]`.
#'
Cluster_SimTem <- function(binsMatrix) 
{
  fucStep <- paste0(" Cluster_SimTem()_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  message("Making similarity template ... ")

  totalcells <- ncol(binsMatrix)
  num_bins <- nrow(binsMatrix)
  result <- sapply(1:totalcells, function(i) {
    sapply(i:totalcells, function(j) {
      sum(binsMatrix[, i] == binsMatrix[, j]) / num_bins
    })
  })

  similarity <- matrix(0, nrow = totalcells, ncol = totalcells)
  for (i in 1:totalcells) {
    similarity[i, i:totalcells] <- result[[i]]
    similarity[i:totalcells, i] <- result[[i]]
  }

  SIM <- list(IDorder = colnames(binsMatrix),
              similarity = as.matrix(similarity))

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(SIM)
}
