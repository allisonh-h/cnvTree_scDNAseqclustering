#' NEW_pqArm_CN
#'
#' Transform CN matrix to Total Deletion/Loss/Neu/Amp and based on Arm level to
#' smooth the CN; this function processes single-cell copy number data by smoothing
#' copy number variations (CNVs) at the arm level. It utilizes clustering results
#' and cytoband information to assign copy number states (Total Deletion, Loss, 
#' Neutral, or Amplification) for each chromosomal arm.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Cluster_label An integer specifying the cluster for which the copy number 
#'  smoothing is performed.
#' @param Clustering_output A data frame recording the clustering results for each 
#'  cell, including the clustering history at each step.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#' 
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#' @return A matrix with smoothed copy number states (Total Deletion, Loss, 
#'  Neutral, or Amplification) at the arm level across chromosomes.
#'  
NEW_pqArm_CN <- function(input, Cluster_label, Clustering_output, pqArm_file)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  cluster <- getOption("cluster")

  #--------------------- start below ---------------------
  
  Clustering_output <- data.frame(cluster = cluster, cell = names(input))
  selected_files <- subset(Clustering_output, 
                           Clustering_output$cluster %in% c(Cluster_label))
  CN_matrix_temp <- NEW_CN_template(input = input, pqArm_file = pqArm_file)
  CN_matrix <- NEW_CN_seq(input = input, Template = selected_files$cell)
  CN_matrix <- NEW_pqArm_DelNeuAmp(matrix = CN_matrix)

  ### calculate bin position: start & end
  CN_binsLevel <- CN_matrix_temp %>%
                  dplyr::group_by(.data$chr, .data$arm) %>%
                  dplyr::summarise(Freq = dplyr::n(), .groups = "drop") %>%
                  as.data.frame()
  CN_binsLevel$start <- sapply(1:nrow(CN_binsLevel), function(x) {
                               sum(CN_binsLevel$Freq[1:x-1]) + 1 })
  CN_binsLevel$end <- sapply(1:nrow(CN_binsLevel), function(x) {
                             sum(CN_binsLevel$Freq[1:x]) })
  
  ### Binned Smoothing for Single-Cell Copy Number Data
  get_mode <- function(x, bin_idx, col_idx) 
  {
    if (all(is.na(x))) return(NA)
    ux <- sort(unique(x)) # LH: forced to be faithful to its' Code
    counts <- tabulate(match(x, ux))
    
    if (length(counts) > 1) {
      sorted_counts <- sort(counts, decreasing = TRUE)
      if ((sorted_counts[1] == sorted_counts[2]) && cnvTree_msg == TRUE) {
        warning(sprintf("Tie found in Chr Region %d, Cell#›› %d!", bin_idx, col_idx))
      }
    }
    return(ux[which.max(counts)])
  }
  result_list <- lapply(1:nrow(CN_binsLevel), function(i) {
                        pq_CNmatrix <- CN_matrix[CN_binsLevel$start[i]:CN_binsLevel$end[i],
                                                 , drop = FALSE]
                        row_modes <- sapply(1:ncol(pq_CNmatrix), function(j) {
                                            get_mode(pq_CNmatrix[, j], i, j) })
                        return(row_modes) })
  Smooth_pqCN <- do.call(rbind, result_list)
  ### rename Smooth_pqCN
  n_chrs <- length(levels(CN_matrix_temp$chr))
  Smooth_pqCN <- Smooth_pqCN %>%
                 as.data.frame() %>%
                 stats::setNames(c(colnames(CN_matrix))) %>%
                 `rownames<-`(paste0(rep(levels(CN_matrix_temp$chr), each = 2), 
                                     rep(c("p", "q"), times = n_chrs)))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Smooth_pqCN)
}


#' NEW_CN_template
#'
#' Generate a copy number segment template based on chromosomal arms; this function 
#' constructs a template for copy number segmentation using cytoband information 
#' from Giemsa-stained chromosomes.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'   a single cell.
#' @param pqArm_file In-build cytoband template for selection: 
#'  `hg38`, `hg19`, `mm10`, `mm39`. Or a filepath of a table for cytoband information 
#'  seen on Giemsa-stained chromosomes. It should include the following columns:
#'   
#'      - `chrom`: Reference sequence chromosome or scaffold.
#'      - `chromStart`: Start position in genoSeq.
#'      - `chromEnd`: End position in genoSeq.
#'      - `name`: Name of cytogenetic band.
#'      - `gieStain`: Giemsa stain results.
#'
#' @return A data frame containing the defined genomic ranges for p/q arms across 
#'    all chromosomes.
#'    
NEW_CN_template <- function(input, pqArm_file)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  CN_tem <- input[[1]]$bins %>% 
            as.data.frame() %>% 
            dplyr::select(chr = seqnames, start, end, width) # NO chr X,Y

  #============= original pqArm_file.remake() + pqArm_file.pq()=================
  # Add p/q arm information on the template
  # 非p即q
  file_path_Template <- system.file("extdata", 
                                    paste0(pqArm_file, "_arm_map.rds"), 
                                    package = "cnvTree")
  if (file_path_Template == "") {
    stop("The arm map file was not found! Please check the extdata folder.")
  } else {
    pqArm_range_temp <- readRDS(file_path_Template)
  } # chr X,Y
  pqArm_range <- pqArm_range_temp %>% dplyr::filter(.data$arm == "p") # chr X,Y
  #=============================================================================
  
  CN_temp <- vector("list", nrow(pqArm_range)) # a list: optional for chr X,Y

  for (i in 1:nrow(pqArm_range)) {
    CN_temp[[i]] <- CN_tem %>% 
                    dplyr::filter(.data$chr == pqArm_range$chr[i],
                                  .data$start >= pqArm_range$start[i],
                                  .data$end <= pqArm_range$end[i]) %>%
      dplyr::mutate(arm = pqArm_range$arm[i])
  }
  CN_tem_pq <- CN_tem %>% # No chr X,Y
               dplyr::left_join(dplyr::bind_rows(CN_temp), 
                                by = c("chr", "start", "end", "width")) %>%
               tidyr::replace_na(list(arm = "q")) %>%
               dplyr::arrange(.data$chr, .data$start)

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(CN_tem_pq)
}


#' NEW_pqArm_file.pq
#' 
#' Remake pqArm template into new format. Convert UCSC cytoband data to arm-level 
#' chromosomal ranges; this function processes cytoband information from the UCSC 
#' database and reformats it into a structured table containing p/q arm regions 
#' for each chromosome. The output is designed for downstream copy number 
#' variation (CNV) analysis.
#'
#' @param Template A character string specifying the file path to the UCSC cytoband 
#'    data.
#'
#' @return A data frame with labeled chromosomal p/q arm regions, formatted for 
#'    CNV analysis.
#'
NEW_pqArm_file.pq <- function(Template) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.1.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  x <- Template %>%
       dplyr::rename(ChromStart = start, ChromEnd = end) %>%
       dplyr::group_by(.data$arm_category) %>%
       dplyr::summarise(start = min(.data$ChromStart), 
                        end = max(.data$ChromEnd))
  x <- x %>%
       dplyr::mutate(arm = stringr::str_sub(.data$arm_category, -1),
                     chr = stringr::str_sub(.data$arm_category, end = -2)) %>%
       dplyr::select(c("chr", "start", "end", "arm"))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(x)
}


#' NEW_pqArm_DelNeuAmp
#' 
#' Transform CN matrix to Del/Neu/Amp four types. Categorize copy number variations 
#' into three types; this function classifies copy number variations (CNVs) into 
#' three discrete categories: Total Deletion (CN < 0.5), Loss (CN < 2), Neutral 
#' (CN = 2), and Amplification (CN > 2). The input matrix represents copy number 
#' data across genomic regions for multiple cells.
#'
#' @param matrix An integer matrix where columns represent individual cells, and 
#'    rows correspond to fixed-bin size genomic regions across all chromosomes.
#'
#' @return An integer matrix containing only values 0, 1, and 2, representing:
#' 
#'   - `0`: Total Deletion (CN < 0.5)
#'   - `1`: Loss (0.5 <= CN < 2)
#'   - `2`: Neutral (CN = 2)
#'   - `3`: Amplification (CN > 2)
#'
NEW_pqArm_DelNeuAmp <- function(matrix) ## function: 3.1.2 ##
{
  new_matrix <- base::matrix(NA, nrow(matrix), ncol(matrix)) 
  
  # strategy 1
  #new_matrix[matrix <= 0 & matrix < 1] <- 0    # Deletion
  #new_matrix[matrix >= 1 & matrix < 2] <- 1    # Loss
  #new_matrix[matrix == 2] <- 2                 # Neutral
  #new_matrix[matrix > 2 & matrix <= 3]  <- 3   # Gain
  #new_matrix[matrix > 3] <- 4                  # Amplification
  
  # strategy 2
  new_matrix[matrix < 2] <- 1                 # Total Deletion
  new_matrix[matrix == 2] <- 2                # Neutral
  new_matrix[matrix > 2]  <- 3                # Gain
  
  dimnames(new_matrix) <- dimnames(matrix)
  
  return(new_matrix)
}


#' pqArm_file.cen
#' 
#' Remake centromere template into new format (X). Process UCSC cytoband data 
#' for "acen" and "gvar" regions; this function extracts and formats cytoband 
#' information from the UCSC database, specifically for the "acen" (centromeric) 
#' and "gvar" (variable heterochromatic) cytoband types. The output is structured 
#' for defining masking ranges to exclude copy number variations (CNVs) from 
#' downstream analyses.
#'
#' @param FILE A character string specifying the file path to the UCSC cytoband 
#'  data file.
#'
#' @return A data frame containing labeled genomic ranges for "acen" and "gvar" 
#'  cytoband regions, which can be used for masking CNVs in experimental 
#'  analyses.
#'
pqArm_file.cen <- function(FILE)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.1.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  sexchromosome <- getOption("sexchromosome")
  
  #--------------------- start below ---------------------

  if (FILE == "hg38" | FILE == "hg19" | FILE == "mm10" | FILE == "mm39") {
    filename <- paste0(FILE, "_cytoBand.txt.gz")
    FILE <- system.file("extdata", filename, package = "cnvTree")
  }
  x <- utils::read.table(gzfile(FILE), sep = "\t", 
                         col.names = c("chr", "ChromStart", "ChromEnd", "name", 
                                       "gieStain"))
  x <- x %>% dplyr::filter(!grepl("_", .data$chr))
  # set chr levels
  vec <- unique(x$chr)
  nums <- as.numeric(gsub("chr", "", vec)[grepl("\\d", vec)])
  nums <- paste0("chr", nums[order(nums)])
  Levels <- c(nums, vec[!grepl("\\d", vec)])
  x <- x %>%
       dplyr::mutate(cen_category = paste0(.data$chr, .data$gieStain)) %>%
       dplyr::group_by(.data$cen_category) %>%
       dplyr::summarise(Start = min(.data$ChromStart),
                        End  = max(.data$ChromEnd)) %>%
       as.data.frame()
  x <- x %>%
    dplyr::filter(.data$cen_category %in% paste0(rep(Levels,each = 2), c("acen","gvar"))) %>%
    dplyr::mutate(cen = stringr::str_sub(.data$cen_category, -4),
                  chr = stringr::str_sub(.data$cen_category, end = -5)) %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(MaskStart = min(.data$Start),
                     MaskEnd  = max(.data$End)) %>%
    dplyr::arrange(chr = factor(.data$chr, levels = Levels), .data$MaskStart) %>%
    dplyr::select(c("chr", "MaskStart", "MaskEnd"))
  
  if (sexchromosome == FALSE) {
    rows_to_remove <- grepl("chrX|chrY|chrM", x$chr)
    x <- x[!rows_to_remove, ]
  }

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(x)
}


#' pqArm_clustering
#' 
#' Cluster cells based on arm-level copy number patterns; this function performs 
#' clustering on cells using arm-level copy number variations (CNVs). It groups 
#' cells into clusters based on chromosomal arm-level CNV profiles, providing 
#' a hierarchical clustering history at each step.
#'
#' @param matrix An integer matrix where columns represent individual cells, and 
#'  rows correspond to arm-level copy number regions across chromosomes.
#' @param Label An integer specifying the cluster to compute in the arm-level CNV 
#'  analysis.
#'
#' @return A table recorded the clustering result for each cell. This table recorded 
#'  the clustering history in each step.
#'  
pqArm_clustering <- function(matrix, Label)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  sexchromosome <- getOption("sexchromosome")
  
  #--------------------- start below ---------------------
  
  cluster <- sapply(1:ncol(matrix), function(x) {
    paste(matrix[, x], collapse = "_")
  })
  cluster <- tibble::tibble(pqArm_pattern = cluster,
                            cellID = colnames(matrix),
                            cluster = Label)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(cluster)
}


#' NEW_pqArm_clustering_summary
#' 
#' Summarize pqArm clustering step results; this function provides a summary of 
#' the pqArm clustering process, detailing the distribution of arm-level copy 
#' number variation (CNV) patterns across different clusters.
#'
#' @param matrix A data frame recording the clustering results for each cell,
#'  including the clustering history at each step.
#' @param Label An integer specifying the cluster to compute in the arm-level 
#'  CNV analysis.
#'
#' @return A data frame summarizing the pqArm clustering step, containing the 
#'  following columns:
#'   - `pqArm_pattern`: The identified copy number patterns at the arm level.
#'   - `pqArm_pattern_cellnum`: The number of cells associated with each pqArm 
#'                              pattern.
#'   - `cluster`: The assigned cluster for each pattern.
#'   - `pqArm_cluster`: The final cluster grouping based on pqArm patterns.
#'
NEW_pqArm_clustering_summary <- function(matrix, Label)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 3.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  cluster_table <- table(matrix$pqArm_pattern) %>%
                   as.data.frame() %>%
                   dplyr::arrange(dplyr::desc(.data$Freq))
  cluster_table$cluster <- Label
  cluster_table$pqArm_cluster <- seq_len(nrow(cluster_table))
  cluster_table <- cluster_table %>% dplyr::rename(pqArm_pattern = Var1, 
                                                   pqArm_cellnum = Freq)

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(cluster_table)
}


#' NEW_pqArm_file.remake
#'
#' Extract chromosome arm ranges from UCSC cytoband data; this function processes 
#' a UCSC cytoband file to extract the ranges of chromosome long and short arms, 
#' excluding the centromere regions.
#'
#' @param FILE A character string specifying the file path to the UCSC cytoband 
#'    file.
#'
#' @return A table containing the ranges of chromosomal arms, excluding centromeric 
#'   regions. It includes cytobands of type `acen` and `gvar`, along with the 
#'   ranges of the preceding and following cytobands.
#'
NEW_pqArm_file.remake <- function(FILE) # function 3.4
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" NEW_pqArm_file.remake_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  supported_genomes <- c("hg38", "hg19", "mm10", "mm39")
  
  if (FILE %in% supported_genomes) {
    filename <- paste0(FILE, "_cytoBand.txt.gz")
    FILE <- system.file("extdata", filename, package = "cnvTree")
  }
  x <- read.table(FILE, 
                  sep="\t", 
                  col.names = c("chr", "start", "end", "name","gieStain"))
  x <- x[grepl("^chr([0-9]{1,2}|[XY])$", x$chr), ]
  chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  x$chr <- factor(x$chr, levels = chrom_levels)
  x <- x[order(x$chr, x$start), ]
  x <- x %>% mutate(start = .data$start + 1,
                    arm = substring(.data$name, 1, 1),
                    arm_category = paste0(.data$chr, .data$arm))
  Sum_x <- x %>% 
           dplyr::group_by(.data$arm_category) %>%
           dplyr::slice_head(n = 1) %>%
           as.data.frame()
  Sum_x <- x %>% 
           dplyr::group_by(.data$arm_category) %>%
           dplyr::slice_tail(n = 1) %>%
           as.data.frame() %>%
           rbind(Sum_x) %>%
           dplyr::mutate(chr = factor(.data$chr, levels = chrom_levels)) %>% 
           dplyr::arrange(.data$chr, .data$start)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(Sum_x)
}

