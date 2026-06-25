#' scDNA.superimpose
#' 
#' Superimpose defined CNVs onto scDNA-seq copy number results; this function 
#' overlays high-confidence copy number variations (CNVs) onto single-cell DNA 
#' sequencing (scDNA-seq) clustering results to determine whether each cluster 
#' contains the corresponding CNVs.
#'
#' @param Template A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, 
#'   re-clustering, and subclone clustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template 
#'   and its copy number. It includes the following columns:
#'   
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'     
#' @param DefinedCNVs A data frame containing high-confidence CNV regions, with 
#'  the following columns:
#'  
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or 
#'           "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#'
#' @return The function returns the updated `Template` list with an additional table:
#'   - `superimpose`: A data frame where `Template$Subclone_CN` overlaps with 
#'   `DefinedCNVs`, calculating whether each cluster contains or lacks the 
#'   corresponding CNV.
#'
scDNA.superimpose <- function(Template, DefinedCNVs)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  Groups <- unique(Template$Subclone_CN$Subclone)
  Template$Subclone_CN$CNV_state <- Template$Subclone_CN$CN
  # CN only seperate in 3 types: del/neu/amp
  Template$Subclone_CN$CN = dplyr::case_when(
    Template$Subclone_CN$CN <  2 ~ "del",
    Template$Subclone_CN$CN == 2 ~ "neu",
    Template$Subclone_CN$CN >  2 ~ "amp")
  
  superimpose <- NULL
  for (groups in 1:length(Groups)) {
    intersection <- NULL
    # check defined CNVs in each group
    
    CNVs <- Template$Subclone_CN %>%
            dplyr::filter(.data$Subclone %in% c(Groups[groups]),
                          .data$chr %in% c(unique(DefinedCNVs$chr)),
                          .data$CN %in% c(DefinedCNVs$CN))
    intersection <- merge(DefinedCNVs, CNVs, by = c("chr", "CN"))
    intersection <- intersection %>%
      dplyr::mutate(F_start = dplyr::case_when(.data$start < .data$CNV_start ~ 0,
                                               .data$start >= .data$CNV_start & 
                                                  .data$start <= .data$CNV_end ~ 1,
                                               .data$start > .data$CNV_end ~ 2),
                    F_end = dplyr::case_when(.data$end < .data$CNV_start ~ 0,
                                             .data$end >= .data$CNV_start & 
                                               .data$end <= .data$CNV_end ~ 1,
                                             .data$end > .data$CNV_end ~ 2),
                    seg = paste0(.data$F_start, .data$F_end),
                    final_start = dplyr::case_when(.data$seg %in% c("00", "22") 
                                                      ~ NA,
                                                   .data$seg %in% c("01", "02") 
                                                      ~ .data$CNV_start,
                                                   .data$seg %in% c("11", "12") 
                                                      ~ .data$start),
                    final_end = dplyr::case_when(.data$seg %in% c("00", "22") 
                                                    ~ NA,
                                                 .data$seg %in% c("01", "11") 
                                                    ~ .data$end,
                                                 .data$seg %in% c("02", "12") 
                                                    ~ .data$CNV_end)) %>%
      dplyr::filter(!.data$seg %in% c("00", "22")) %>%
      dplyr::mutate(cnv_range = .data$final_end - .data$final_start + 1) %>%
      dplyr::group_by(.data$CNV_region) %>%
      dplyr::summarise(cnv_range = sum(.data$cnv_range)) %>%
      as.data.frame()
    
    superimpose <- dplyr::left_join(DefinedCNVs, intersection, by = "CNV_region") %>%
                   dplyr::mutate(Subclone = Groups[groups],
                                 CNV_range = .data$CNV_end - .data$CNV_start + 1,
                                 cnv_range = ifelse(is.na(.data$cnv_range) == TRUE, 
                                                    0, .data$cnv_range),
                                 cnv_ratio = .data$cnv_range /.data$CNV_range,
                                 final_cnv = ifelse(.data$cnv_ratio >= 0.5, 1, 0)) %>%
                                 rbind(.data, superimpose)
  }
  Template$superimpose <- superimpose

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Template)
}


#' scDNA.clustering
#' 
#' Function for receiving DNA clustering output in superimpose range. Generate 
#' superimposition results between scDNA clustering and defined CNVs; this function 
#' calculates the overlap between single-cell DNA sequencing (scDNA-seq) 
#' clustering results and high-confidence defined copy number variations (CNVs), 
#' indicating whether each cluster contains the corresponding CNVs.
#'
#' @param Template A data frame where `Template$Subclone_CN` records each subclone's 
#'  unique chromosome segment template and its copy number.
#'  
#'     It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'
#' @return The function updates `Template` by adding a new matrix, `Template$DNA_cluster`, 
#'  which includes:
#'   - `DNA_cluster`: The index for final clustering result.
#'   - `DNA_Cellnum`: The number of cells in each cluster.
#'   - `DefinedCNVs`: A binary matrix (1 or 0) indicating whether each cluster 
#'                    contains the corresponding CNV.
#'
scDNA.clustering <- function(Template)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  cnv_region <- unique(Template$superimpose$CNV_region)
  cnv_region <- paste0("CNV", cnv_region)
  Subclone_ss <- unique(Template$superimpose$Subclone)
  
  cnv_matrix <- matrix(nrow = length(Subclone_ss), ncol = length(cnv_region) + 2)
  for (subclone in 1:length(Subclone_ss)) {
    Cell_num <- Template$final_cluster_output %>%
                dplyr::filter(.data$Subclone_cluster %in% Subclone_ss[subclone]) %>%
                nrow()
    CNVs <- Template$superimpose %>%
            dplyr::filter(.data$Subclone %in% Subclone_ss[subclone]) %>%
            dplyr::select(.data$final_cnv) %>%
            dplyr::pull()
    
    cnv_matrix[subclone, ] <-c(CNVs, Subclone_ss[subclone], Cell_num[1])
    
  }
  colnames(cnv_matrix) <- c(cnv_region, "DNA_cluster", "DNA_Cellnum")
  rownames(cnv_matrix) <- seq_len(nrow(cnv_matrix))
  Template$DNA_cluster <- cnv_matrix
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Template)
}


#' Totalcluster_pdf
#' 
#' For creating pdf in total clusters by CN matrix. Generate copy number profiles 
#' for scDNA-seq cell clustering results in PDF file; this function outputs copy 
#' number profiles from single-cell DNA sequencing (scDNA-seq) clustering results, 
#' saving the visualization as a PDF file.
#'
#' @param Input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Template A table recorded the clustering, pqArm clustering, re-clustering, 
#'  and subclone clustering step result for each cell. This table recorded the 
#'  clustering history in each step.
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
#' @param cellcutoff A numeric value defining the minimum number of cells required 
#'  for a cluster to be included.
#' @param step A character string specifying the name of the output clustering 
#'  step to select. Valid options are "pqArm," "Recluster," and "Subclone." 
#'  Please ensure the name matches one of these options (default: "Subclone").
#' @param FILEname A character string specifying the name of the output PDF file.
#' @param FILEpath A character string specifying the file path where the output 
#'  PDF will be saved.
#' @param sexchromosome A logical value. If `TRUE`, the output plot includes 
#'  copy number information for sex chromosomes. Defaults to `FALSE`.
#'
#' @return A copy number profile visualization saved as a PDF file.
#'
Totalcluster_pdf <- function(Input, Template, pqArm_file, cellcutoff, step = "Subclone", 
                             FILEname, FILEpath, sexchromosome)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.3_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  cluster_name = paste0(step, "_cluster")
  cellnum_name = paste0(step, "_cellnum")
  
  Cluster_No <- Totalcluster_Cluster_No(Template = Template, 
                                        cellnum_name = cellnum_name, 
                                        cellcutoff = cellcutoff, 
                                        cluster_name = cluster_name)
  Fig_seq <- list()
  height_ratio <- NULL
  cellnum_list <- NULL
  count = 0
  for (k in Cluster_No) {
    cat("Making ", FILEname, ":", "heatmap of", cluster_name, k, "\n")
    SS <- Totalcluster_SS(Template = Template, cluster_name = cluster_name, k = k)
    Cell_num <- length(SS)
    cellnum_list <- c(cellnum_list, Cell_num)
    
    # height_size <- Cell_num*10
    if (Cell_num <= 50) {
      height_size <- 300
      height_ratio <- c(height_ratio, 3)
    } else if (Cell_num > 50 && Cell_num <= 100) {
      height_size <- 500
      height_ratio <- c(height_ratio, 5)
    } else {
      height_size <- 1000
      height_ratio <- c(height_ratio, 8)
    }
    count = count + 1
    Fig_seq[[count]] <- GenomeHeatmap(Input = Input, 
                                      cellID = SS, 
                                      pqArm_file = pqArm_file, 
                                      sexchromosome = sexchromosome)
  }
  Label <- paste0(Cluster_No, " (n = ", cellnum_list, ")")
  
  combine_plot <- ggpubr::ggarrange(plotlist = Fig_seq,
                                    ncol = 1,
                                    labels = Label,
                                    common.legend = FALSE,
                                    legend = "right",
                                    hjust  = 0.8,
                                    align = "v",
                                    font.label = list(size = 45, face = "bold", 
                                                      color ="black"),
                                    heights = c(height_ratio)) +
  ggplot2::theme(plot.margin = ggplot2::margin(2,2,2,8, "cm"))
  
  ggplot2::ggsave(combine_plot,
                  filename = paste0(FILEpath, FILEname),
                  height = (300*sum(height_ratio) + 800*length(height_ratio)),
                  width = 13000,
                  units = "px",
                  limitsize = FALSE)
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}


#' Totalcluster_Cluster_No
#' 
#' Transfer the number of clusters in the data. Identify clusters meeting the 
#' thershold by number of cell; his function filters clusters based on the number 
#' of cells, returning a list of clusters that meet the specified threshold.
#'
#' @param Template  A data frame recording the pqArm clustering, re-clustering, 
#'  nd subclone clustering results for each cell. This table tracks the clustering 
#'  history at each step.
#' @param cellnum_name A character string specifying the column in `Template` that 
#'  records the number of cells in each cluster.
#' @param cellcutoff A numeric value defining the minimum number of cells required 
#'  for a cluster to be included.
#' @param cluster_name A character string specifying the column in `Template` that 
#'  records the clustering result to be filtered.
#'
#' @return A character vector containing cluster names that meet the specified 
#'  cell number threshold.
#'
Totalcluster_Cluster_No <- function(Template, cellnum_name, cellcutoff, cluster_name)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.3.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  Cluster_No <- Template %>%
                dplyr::filter(.data[[cellnum_name]] >= cellcutoff) %>%
                dplyr::pull(.data[[cluster_name]])
  Cluster_No <- sort(unique(Cluster_No))
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Cluster_No)
}


#' Totalcluster_SS
#' 
#' Transfer each cluster of cells in the data. Retrieve cell IDs from a specific 
#' cluster; this function returns a list of cell IDs belonging to a designated 
#' cluster based on the clustering results recorded in the `Template`.
#'
#' @param Template A data frame recording the pqArm clustering, re-clustering, 
#'  and subclone clustering results for each cell. This table tracks the clustering 
#'  history at each step.
#' @param cluster_name A character string specifying the column in `Template` that 
#'  records the clustering result to be queried.
#' @param k An integer specifying the designated cluster for which cell IDs should 
#'  be retrieved.
#'
#' @return A character vector containing the cell IDs that belong to the designated cluster.
#'
Totalcluster_SS <- function(Template, cluster_name, k)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.3.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  SS <- Template %>%
        dplyr::filter(.data[[cluster_name]] == k) %>%
        dplyr::pull(.data$cellID)
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(SS)
}


#' GenomeHeatmap
#' 
#' Plot copy number pattern for each cell in heatmap; this function visualizes 
#' the copy number variation (CNV) pattern for each cell, using cytoband information 
#' to annotate chromosomal regions. It supports optional inclusion of sex chromosome 
#' CNV data in the output plot.
#'
#' @param Input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param cellID A character vector specifying the `cellID`s of the cells for 
#'  heatmap plotting.
#' @param pqArm_file A table for cytoband information seen on Giemsa-stained chromosomes.
#' It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param sexchromosome Logical. If `TRUE`, the plot includes copy number information 
#'  for sex chromosomes. Defaults to `FALSE`.
#'
#' @return A `ggplot` object displaying the copy number pattern for each cell.
#'
GenomeHeatmap <- function(Input, cellID, pqArm_file, sexchromosome)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.3.3_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)

  # import CN template and total cell copy number matrix
  CN_bins_template <- NEW_CN_template(input = Input, pqArm_file = pqArm_file)
  CN_chr_template <- CN_bins_template %>%
                     dplyr::group_by(.data$chr) %>%
                     dplyr::summarise(length = max(.data$end) - min(.data$start)) %>%
                     dplyr::mutate(X_cum = c(0, cumsum(as.numeric(length[-length(length)]))))

  # Import the HC clustering cell order
  if (length(cellID) == 1) {
    cellOrder <- cellID
    Input <- Input[cellOrder]
  } else {
    cellOrder <- clusterbyHMM(input = Input, selected = cellID)
    cellOrder <- cellOrder$IDorder %>%
                 as.data.frame() %>%
                 tibble::rownames_to_column(var = "cellID")
    Input <- Input[cellOrder$cellID]
  }

  if (config_hid$msg == TRUE) {
    print(paste0("LH: start looping...6.3.3.1_cnvTree_", config_hid$v_num)) # LH: added on 02232026
  }

  # Import copy number data
  numofcell <- length(cellID)
  mat_input <- lapply(1:numofcell, function(i) segment_transform(data = Input[[i]], 
                                                                 index = i, 
                                                                 CN_chr_template = CN_chr_template))
  mat_input <- dplyr::bind_rows(mat_input)

  # chrX, chrY in plotting output is optional
  #if (sexchromosome == FALSE) {
  #  CN_bins_template <- CN_bins_template %>% 
  #                      dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
  #  CN_chr_template <- CN_chr_template %>% 
  #                     dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
  #  mat_input <- mat_input %>% 
  #               dplyr::filter(!.data$chr %in% c("chrX", "chrY"))
  #}
  # mat_input$CN[which(mat_input$CN>5)] <- 5

  # plot chr site
  Chr_tmp <- CN_chr_template %>%
             dplyr::mutate(chr = sub("^chr", "", .data$chr),
                           X_cum_end = .data$length +.data$ X_cum,
                           text_pos = (.data$length + .data$X_cum*2)/2,
                           Y_end = numofcell-0.8)

  # Fixed color template: Copy number more than 5 as same color
  color_tem <- generate_dynamic_colormap(data_matrix = mat_input$CN)

  # figure legend
  lgd_labels <- sort(unique(mat_input$CN))
  #lgd_labels[which(lgd_labels==5)] <- ">= 5"

  height = dplyr::case_when(numofcell > 60 ~ numofcell,
                            numofcell <= 60 & numofcell > 10 ~ 50,
                            numofcell <= 10 ~ 30)

  segment_length = dplyr::case_when(numofcell > 201 ~ 2,
                                    numofcell >= 101 & numofcell <= 200 ~ 4,
                                    numofcell >= 41 & numofcell <= 100 ~ 10,
                                    numofcell >= 21 & numofcell <= 40 ~ 15,
                                    numofcell >= 1 & numofcell <= 20 ~ 20)

  PlotCN_heatmap <-
    ggplot2::ggplot(mat_input) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$X_cum_start, 
                                       xend = .data$X_cum_end, 
                                       y = .data$Y_cum_start, 
                                       yend = .data$Y_cum_end,
                                       color = factor(.data$CN)), 
                          linewidth = segment_length, show.legend = TRUE) +
    ggplot2::scale_color_manual(name = "Copy Number", 
                                labels = lgd_labels, 
                                values = color_tem,
                                guide = ggplot2::guide_legend(override.aes = list(linewidth = 6))) +
    ggplot2::geom_segment(data = Chr_tmp, 
                          ggplot2::aes(x = .data$X_cum_end, 
                                       y = -0.2, 
                                       xend = .data$X_cum_end, 
                                       yend = .data$Y_end)) +
    ggplot2::geom_text(data = Chr_tmp, 
                       ggplot2::aes(x = .data$text_pos, 
                                    y = -0.05 * height, 
                                    label = .data$chr), 
                       size = 10) +
    ggplot2::labs(x = "Chromosome", y = "Cell ID", color = "Copy Number") +
    ggplot2::xlim(c(0, max(Chr_tmp$X_cum_end)+1)) +
    ggplot2::ylim(c(-0.05 *height, numofcell)) +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 2, 3, 1, "cm"),
                   plot.background = ggplot2::element_rect(fill = "white"),
                   panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::theme_void()

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(PlotCN_heatmap)
}


#' segment_transform
#' 
#' Transform copy number segment information for each cell; this function transforms 
#' copy number segment information for a single cell, organizing continuous regions 
#' with the same copy number into unified segments. It uses chromosome length and 
#' cumulative length data to standardize genomic positions.
#'
#' @param data A `GRanges` object representing copy number segments for a single 
#'  cell.
#' @param index A numeric value indicating the order or ID of the cell. Used for 
#'  labeling or indexing the output.
#' @param CN_chr_template A data frame containing chromosome information, including 
#'  chromosome name, chromosome length, and cumulative length.
#'
#' @return A data frame recording copy number segments for the cell, with continuous 
#'  regions of the same copy number merged.
#'
segment_transform <- function(data, index, CN_chr_template) # Function 6.3.3.1 #
{
  D1 <- as.data.frame(data$bins) %>% 
        dplyr::mutate(segment = rep(seq_along(rle(.data$copy.number)$values), 
                                              rle(.data$copy.number)$lengths)) %>%
        dplyr::group_by(.data$seqnames, .data$segment, .data$copy.number)
  
  D2 <- D1 %>% dplyr::summarize(start = gdata::first(.data$start),
                                end =  gdata::last(.data$end),
                                width = (.data$end) - (.data$start) + 1,
                                strand = gdata::first(.data$strand),
                                .groups = "drop") %>%
               dplyr::ungroup() %>%
               dplyr::select(-c(.data$segment, .data$strand)) %>%
               dplyr::rename("chr" = "seqnames", "CN" = "copy.number")
  
  D3 <- D2 %>% dplyr::left_join(CN_chr_template, by = "chr") %>%
               dplyr::mutate(X_cum_start = .data$start + .data$X_cum,
                             X_cum_end = .data$end + .data$X_cum,
                             Y_cum_start = index - 1,
                             Y_cum_end = index - 1)
  
  return(D3)
}


#' scDNA_CNVpattern
#' 
#' Cluster w/ or w/o CNV pattern plot. Generate a heatmap of high-confidence 
#' CNV patterns with Dendrogram; this function creates a heatmap displaying 
#' high-confidence copy number variation (CNV) patterns across clusters, with 
#' hierarchical clustering represented by a dendrogram.
#'
#' @param Input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param final_cluster A data frame containing high-confidence CNV regions, with 
#'  the following columns:
#'  
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or 
#'   "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#'   
#' @param cellcutoff A numeric value defining the minimum number of cells required 
#'  for a cluster to be included.
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
#' @param FILEpath A character string specifying the directory where the output 
#'  file will be saved.
#' @param FILEname A character string specifying the name of the output `.png` file.
#' @param sexchromosome A logical value. If `TRUE`, the heatmap includes copy number 
#'  information for sex chromosomes.
#' @param smoothheatmap A logical value. If `TRUE`, the heatmap applies smoothing over 
#'  a 10⁶ bp range in chromosome copy number data.
#'
#' @return A PNG file containing a heatmap of defined CNVs across clusters.
#'
scDNA_CNVpattern <- function(input, final_cluster, cellcutoff, pqArm_file, 
                             FILEpath, FILEname, sexchromosome, smoothheatmap)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.4_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  smoothheatmap <- config_hid$smoothheatmap
  
  if (smoothheatmap == TRUE) {
    # bins template
    bins_window = input[[1]]$bins
    bins_window$copy.number <- NULL
    
    # chromosome seperated in fixed-bin size 
    chr_df = circlize::read.chromInfo(species = pqArm_file)$df 
    
    if (sexchromosome == FALSE) {
      chr_df <- chr_df[!grepl("chrX|chrY", chr_df$chr), ]
    } else {
      chr_df = chr_df[chr_df$chr %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
    }
    chr_gr = GenomicRanges::GRanges(seqnames = chr_df[, 1], 
                                    ranges = IRanges::IRanges(chr_df[, 2] + 1, 
                                                              chr_df[, 3] + 1))
    chr_window = EnrichedHeatmap::makeWindows(chr_gr, w = 1e6, short.keep = T)
    mtch = as.data.frame(IRanges::findOverlaps(chr_window, bins_window, type = "any"))
    
    # 先平均出每個bin 的copy number ，得到average sequence
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)
    final_cluster <- final_cluster %>%
                     tidyr::unite("Subclone_no", Subclone_cluster, Subclone_cellnum, 
                                  sep = " (n = ", remove = FALSE) %>%
                     mutate(Subclone_no = paste0(Subclone_no, ")"))
    subclone_no <- unique(final_cluster$Subclone_cluster)
    
    num_mat <- NULL
    for (i in 1:length(subclone_no)) {
      cellID <- final_cluster %>% 
                dplyr::filter(.data$Subclone_cluster %in% subclone_no[i]) %>% 
                dplyr::pull(.data$cellID)
      Subclone_CN <- NEW_CN_seq(input = input, Template = cellID) #LH_02102025: modified
      averageCN <- round(apply(Subclone_CN, 1, mean), digits = 0)  # mean() for cluster of cells
      smooth_CN = rep(2, length(chr_window))
      recalculateCN = mtch %>%
                      dplyr::mutate(copy_number = averageCN[.data$subjectHits]) %>%
                      dplyr::group_by(.data$queryHits) %>%
                      dplyr::summarise(mean_copy_number = round(mean(.data$copy_number, 
                                                                     na.rm = TRUE), 
                                                                digits = 0))
      smooth_CN[as.numeric(recalculateCN$queryHits)] <- recalculateCN$mean_copy_number
      
      if(is.null(num_mat)) {
        num_mat <- smooth_CN
        num_mat <- as.matrix(num_mat) ## LH added 122325
      } else {
        num_mat <- cbind(num_mat, smooth_CN)
        num_mat <- as.matrix(num_mat) ## LH added 122325
      }
    }
    colname_num_mat <- final_cluster %>% ## LH added 051226
                       filter(Subclone_cluster %in% subclone_no) %>%
                       select(Subclone_cluster, Subclone_no) %>%
                       distinct() %>% 
                       arrange(as.numeric(Subclone_cluster))
    colnames(num_mat) <- colname_num_mat$Subclone_no
    Subclone_no <- colnames(num_mat) # LH added 051326
    
    chr_sort <- as.data.frame(chr_window) %>% 
                dplyr::select(.data$seqnames, .data$start)
    
    if (sexchromosome == TRUE) {
      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    } else if(sexchromosome == FALSE) {
      num_mat <- cbind(chr_sort, num_mat) %>%
                 dplyr::filter(!.data$seqnames %in% c("chrX", "chrY", "chrM")) %>%
                 as.matrix()
     
      num_mat <- num_mat[, -c(1, 2)]
      num_mat <- as.matrix(num_mat) ## LH added 122325
      num_mat <- apply(num_mat, c(1, 2), as.numeric)
      
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- chr[!(chr %in% c("chrX", "chrY"))]
      chr <- factor(chr, 
                    levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), 
                                            c("chrX", "chrY", "chrM")))
    }
  } else {
    ### without smoothing step ###
    # 先平均出每個bin 的copy number ，得到average sequence
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)

    final_cluster <- final_cluster %>%
                     tidyr::unite("Subclone_no", Subclone_cluster, 
                                  Subclone_cellnum, sep = " (", remove = FALSE) %>%
                     mutate(Subclone_no = paste0(Subclone_no, ")"))
    subclone_no <- unique(final_cluster$Subclone_cluster)

    num_mat <- NULL
    for(i in 1:length(subclone_no)) {
      cellID <- final_cluster %>% 
                dplyr::filter(.data$Subclone_cluster %in% subclone_no[i]) %>% 
                dplyr::pull(.data$cellID)
      Subclone_CN <- NEW_CN_seq(input = Input, Template = cellID)
      num_mat <- cbind(num_mat, round(apply(Subclone_CN, 1, mean), digits = 0))
    }
    
    chr_window = Input[[1]]$bins
    chr_window$copy.number <- NULL
    
    colname_num_mat <- final_cluster %>% ## LH added 051226
                       filter(Subclone_cluster %in% subclone_no) %>%
                       select(Subclone_cluster, Subclone_no) %>%
                       distinct() %>% 
                       arrange(as.numeric(Subclone_cluster))
    colnames(num_mat) <- colname_num_mat$Subclone_no
    Subclone_no <- colnames(num_mat) # LH added 051326
    
    chr_sort <- chr_window %>%
                as.data.frame() %>%
                dplyr::select(.data$seqnames, .data$start) %>%
                dplyr::arrange(.data$seqnames, .data$start)
    
    if (sexchromosome == TRUE) {
      num_mat <- cbind(chr_sort, num_mat) %>%
                 dplyr::arrange(.data$seqnames, .data$start) %>%
                 as.matrix()
      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    } else if (sexchromosome == FALSE) {
      num_mat <- cbind(chr_sort, num_mat) %>%
                 dplyr::arrange(.data$seqnames, .data$start) %>%
                 dplyr::filter(!.data$seqnames %in% c("chrX", "chrY", "chrM")) %>%
                 as.matrix()
      
      # heatmap annotation labels
      chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
      chr <- chr[!(chr %in% c("chrX", "chrY", "chrM"))]
      chr <- factor(chr, 
                    levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), 
                                            c("chrX", "chrY", "chrM")))
    }
    num_mat <- num_mat[, -c(1, 2)]
    num_mat <- apply(num_mat, c(1, 2), as.numeric)
  }
  # heatmap annotation labels
  chr_level <- unique(sub("^chr", "", chr))
  dynamic_colors <- generate_dynamic_colormap(data_matrix = num_mat)
  
  grDevices::png(filename = paste0(FILEpath, FILEname),
                 width = 2800,
                 height = (length(subclone_no)) * 150)
  
  legend_nrow <- min(length(dynamic_colors), 20)
  heatmap_legend <- ComplexHeatmap::Legend(
                      labels = names(dynamic_colors),  # 根據您的顏色名稱替換
                      legend_gp = grid::gpar(fill = dynamic_colors),
                      title = "Copy number",
                      nrow = legend_nrow)  # 每排顯示10個
  
  Oncoscan <- ComplexHeatmap::Heatmap(t(num_mat), 
                                      name = "Copy number", 
                                      col = dynamic_colors,
                                      column_split = chr,
                                      cluster_columns = FALSE,
                                      cluster_rows = TRUE,
                                      show_row_dend = TRUE,
                                      show_row_names = FALSE,
                                      use_raster = TRUE,
                                      show_heatmap_legend = FALSE,
                                      row_split = subclone_no,
                                      cluster_row_slices = TRUE,
                                      row_dend_width = ggplot2::unit(5, "cm"),
                                      row_dend_gp = grid::gpar(lwd = 2, col = "black"),
                                      row_title = "scDNA clusters",
                                      row_title_gp = grid::gpar(fontsize = 40, fontface = "bold"),
                                      left_annotation = ComplexHeatmap::rowAnnotation(
                                        subgroup = ComplexHeatmap::anno_text(
                                          Subclone_no, rot = 0, gp = grid::gpar(fontsize = 30))), # LH changed 051326
                                      column_title = chr_level,
                                      column_title_side = "bottom",
                                      column_title_gp = grid::gpar(fontsize = 25),
                                      border = TRUE,
                                      column_gap = ggplot2::unit(0, "points")
                                      # row_gap = unit(0, "points")
  )
  ComplexHeatmap::draw(Oncoscan,
                       annotation_legend_side = "right",
                       annotation_legend_list = list(heatmap_legend),
                       padding = ggplot2::unit(c(5, 2, 2, 2), "cm"))
  
  grDevices::dev.off()
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}


#' generate_dynamic_colormap
#' 
#' Color template for heatmaps (Customize color + RColorBrewer template). Assign 
#' colors to copy number states; this function maps numerical copy number states 
#' to a corresponding color scheme for visualization purposes.
#'
#' @param data_matrix A numeric matrix representing copy number states,
#'  where rows correspond to genomic regions and columns correspond to samples or 
#'  clusters.
#'
#' @return A matrix of the same dimensions as `data_matrix`, with each numeric 
#'  copy number state replaced by a corresponding color code.
#'
generate_dynamic_colormap <- function(data_matrix)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.4.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  unique_values <- sort(unique(as.vector(data_matrix)))

  # 固定的顏色對應表，對應數值 0~4 #orange: "#ED7D31"
  predefined_colors <- c("#D0CECE", "#8165A3", "#9BBB59", "#FFC000", "#C0504D") 
  names(predefined_colors) <- 0:4  # 為 0~4 數值建立顏色對應表

  # 定義漸層調色板，用於處理超過 4 的數值
  gradient_palette <- grDevices::colorRampPalette(c("#c51b7d", "#130303"))

  # 切分數值：0 到 4 和 大於 4
  values_below_6 <- unique_values[unique_values <= 4]
  values_above_5 <- unique_values[unique_values > 4]

  # 分配顏色：0~5 的數值使用對應的顏色
  below_6_colors <- predefined_colors[as.character(values_below_6)]

  # 大於 5 的數值使用漸層顏色
  above_5_colors <- gradient_palette(length(values_above_5))

  # 合併顏色映射
  all_colors <- c(below_6_colors, above_5_colors)
  color_mapping <- stats::setNames(all_colors, c(values_below_6, values_above_5))

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(color_mapping)
}


#' NEW_scDNA_CNVpattern
#' 
#' Cluster w/ or w/o CNV pattern plot. Generate a heatmap of high-confidence 
#' CNV patterns with Dendrogram; this function creates a heatmap displaying 
#' high-confidence copy number variation (CNV) patterns across clusters, with 
#' hierarchical clustering represented by a dendrogram.
#'
#' @param Input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param final_cluster A data frame containing high-confidence CNV regions, with 
#'  the following columns:
#'  
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `CNV_region`: The index of CNV regions.
#'   - `CN`: Copy number state, categorized as either "amp" (Amplification) or "del" (Deletion).
#'   - `CNV_start`: Start position of the CNV region.
#'   - `CNV_end`: End position of the CNV region.
#'   - `first_band`: Cytoband label of the first affected band.
#'   - `last_band`: Cytoband label of the last affected band.
#'   
#' @param cellcutoff A numeric value defining the minimum number of cells required 
#'  for a cluster to be included.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, `mm10`, `mm39`.
#'  Or a filepath of a table for cytoband information seen on Giemsa-stained chromosomes.
#'  It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param FILEpath A character string specifying the directory where the output 
#'  file will be saved.
#' @param FILEname A character string specifying the name of the output `.png` file.
#' @param sexchromosome A logical value. If `TRUE`, the heatmap includes copy number 
#'  information for sex chromosomes.
#' @param smoothheatmap A logical value. If `TRUE`, the heatmap applies smoothing over 
#'  a 10⁶ bp range in chromosome copy number data.
#'
#' @return A PNG file containing a heatmap of defined CNVs across clusters.
#'
NEW_scDNA_CNVpattern <- function(input, final_cluster, cellcutoff, pqArm_file, 
                             FILEpath, FILEname, sexchromosome, smoothheatmap)
{
  config_hid_path <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_hid_path)
  fucStep <- paste0(" 6.4_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  smoothheatmap <- config_hid$smoothheatmap
  
  if (smoothheatmap == TRUE) {
    # bins template
    bins_window = input[[1]]$bins
    bins_window$copy.number <- NULL
    
    # chromosome separated in fixed-bin size 
    chr_df = circlize::read.chromInfo(species = pqArm_file)$df 
    
    # Dynamically track targets based on input configuration setting
    if (sexchromosome == FALSE) {
      chr_df <- chr_df[!grepl("chrX|chrY|chrM", chr_df$chr), ]
    } else {
      chr_df = chr_df[chr_df$chr %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]
    }
    
    chr_gr = GenomicRanges::GRanges(seqnames = chr_df[, 1], 
                                    ranges = IRanges::IRanges(chr_df[, 2] + 1, 
                                                              chr_df[, 3] + 1))
    chr_window = EnrichedHeatmap::makeWindows(chr_gr, w = 1e6, short.keep = T)
    mtch = as.data.frame(IRanges::findOverlaps(chr_window, bins_window, type = "any"))
    
    # 先平均出每個bin 的copy number ，得到average sequence
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)
    final_cluster <- final_cluster %>%
      tidyr::unite("Subclone_no", Subclone_cluster, Subclone_cellnum, 
                   sep = " (n = ", remove = FALSE) %>%
      mutate(Subclone_no = paste0(Subclone_no, ")"))
    subclone_no <- unique(final_cluster$Subclone_cluster)
    
    num_mat <- NULL
    for (i in 1:length(subclone_no)) {
      cellID <- final_cluster %>% 
        dplyr::filter(.data$Subclone_cluster %in% subclone_no[i]) %>% 
        dplyr::pull(.data$cellID)
      Subclone_CN <- NEW_CN_seq(input = input, Template = cellID) 
      averageCN <- round(apply(Subclone_CN, 1, mean), digits = 0)  
      smooth_CN = rep(2, length(chr_window))
      recalculateCN = mtch %>%
        dplyr::mutate(copy_number = averageCN[.data$subjectHits]) %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::summarise(mean_copy_number = round(mean(.data$copy_number, 
                                                       na.rm = TRUE), 
                                                  digits = 0))
      smooth_CN[as.numeric(recalculateCN$queryHits)] <- recalculateCN$mean_copy_number
      
      if(is.null(num_mat)) {
        num_mat <- smooth_CN
        num_mat <- as.matrix(num_mat) 
      } else {
        num_mat <- cbind(num_mat, smooth_CN)
        num_mat <- as.matrix(num_mat) 
      }
    }
    colname_num_mat <- final_cluster %>% 
      filter(Subclone_cluster %in% subclone_no) %>%
      select(Subclone_cluster, Subclone_no) %>%
      distinct() %>% 
      arrange(as.numeric(Subclone_cluster))
    colnames(num_mat) <- colname_num_mat$Subclone_no
    Subclone_no <- colnames(num_mat) 
    
    chr_sort <- as.data.frame(chr_window) %>% 
      dplyr::select(.data$seqnames, .data$start)
    
    # Setup heatmap annotations based on user parameters and data structure safety thresholds
    chr <- as.character(sort(GenomicRanges::seqnames(chr_window)))
    if (sexchromosome == FALSE) {
      num_mat <- cbind(chr_sort, num_mat) %>%
        dplyr::filter(!.data$seqnames %in% c("chrX", "chrY", "chrM")) %>%
        as.matrix()
      
      num_mat <- num_mat[, -c(1, 2), drop = FALSE]
      num_mat <- apply(num_mat, c(1, 2), as.numeric)
      
      chr <- chr[!(chr %in% c("chrX", "chrY", "chrM"))]
      chr <- factor(chr, levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), 
                                                 c("chrX", "chrY", "chrM")))
    } else {
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    }
    
  } else {
    ### without smoothing step ###
    final_cluster <- final_cluster %>% dplyr::filter(.data$Subclone_cellnum >= cellcutoff)
    final_cluster <- final_cluster %>%
      tidyr::unite("Subclone_no", Subclone_cluster, 
                   Subclone_cellnum, sep = " (", remove = FALSE) %>%
      mutate(Subclone_no = paste0(Subclone_no, ")"))
    subclone_no <- unique(final_cluster$Subclone_cluster)
    
    num_mat <- NULL
    for(i in 1:length(subclone_no)) {
      cellID <- final_cluster %>% 
        dplyr::filter(.data$Subclone_cluster %in% subclone_no[i]) %>% 
        dplyr::pull(.data$cellID)
      Subclone_CN <- NEW_CN_seq(input = input, Template = cellID) # Fixed standard variable name capitalization
      num_mat <- cbind(num_mat, round(apply(Subclone_CN, 1, mean), digits = 0))
    }
    
    chr_window = input[[1]]$bins
    chr_window$copy.number <- NULL
    
    colname_num_mat <- final_cluster %>% 
      filter(Subclone_cluster %in% subclone_no) %>%
      select(Subclone_cluster, Subclone_no) %>%
      distinct() %>% 
      arrange(as.numeric(Subclone_cluster))
    colnames(num_mat) <- colname_num_mat$Subclone_no
    Subclone_no <- colnames(num_mat) 
    
    chr_sort <- chr_window %>%
      as.data.frame() %>%
      dplyr::select(.data$seqnames, .data$start) %>%
      dplyr::arrange(.data$seqnames, .data$start)
    
    # Superimpose chromosomal definitions conditionally based on incoming structures
    num_mat <- cbind(chr_sort, num_mat) %>%
      dplyr::arrange(.data$seqnames, .data$start)
    
    chr <- as.character(num_mat$seqnames)
    
    if (sexchromosome == FALSE) {
      num_mat <- num_mat %>% 
        dplyr::filter(!.data$seqnames %in% c("chrX", "chrY", "chrM"))
      chr <- as.character(num_mat$seqnames)
      chr <- factor(chr, levels = dplyr::setdiff(levels(GenomicRanges::seqnames(chr_window)), 
                                                 c("chrX", "chrY", "chrM")))
    } else {
      chr <- factor(chr, levels = levels(GenomicRanges::seqnames(chr_window)))
    }
    
    num_mat <- as.matrix(num_mat)
    num_mat <- num_mat[, -c(1, 2), drop = FALSE]
    num_mat <- apply(num_mat, c(1, 2), as.numeric)
  }
  # heatmap annotation labels
  chr_level <- unique(sub("^chr", "", chr))
  dynamic_colors <- generate_dynamic_colormap(data_matrix = num_mat)
  
  grDevices::png(filename = paste0(FILEpath, FILEname),
                 width = 2800,
                 height = (length(subclone_no)) * 150)
  
  legend_nrow <- min(length(dynamic_colors), 20)
  heatmap_legend <- ComplexHeatmap::Legend(
    labels = names(dynamic_colors),  
    legend_gp = grid::gpar(fill = dynamic_colors),
    title = "Copy number",
    nrow = legend_nrow)  
  
  Oncoscan <- ComplexHeatmap::Heatmap(t(num_mat), 
                                      name = "Copy number", 
                                      col = dynamic_colors,
                                      column_split = chr,
                                      cluster_columns = FALSE,
                                      cluster_rows = TRUE,
                                      show_row_dend = TRUE,
                                      show_row_names = FALSE,
                                      use_raster = TRUE,
                                      show_heatmap_legend = FALSE,
                                      row_split = subclone_no,
                                      cluster_row_slices = TRUE,
                                      row_dend_width = ggplot2::unit(5, "cm"),
                                      row_dend_gp = grid::gpar(lwd = 2, col = "black"),
                                      row_title = "scDNA clusters",
                                      row_title_gp = grid::gpar(fontsize = 40, fontface = "bold"),
                                      left_annotation = ComplexHeatmap::rowAnnotation(
                                        subgroup = ComplexHeatmap::anno_text(
                                          Subclone_no, rot = 0, gp = grid::gpar(fontsize = 30))), 
                                      column_title = chr_level,
                                      column_title_side = "bottom",
                                      column_title_gp = grid::gpar(fontsize = 25),
                                      border = TRUE,
                                      column_gap = ggplot2::unit(0, "points")
  )
  ComplexHeatmap::draw(Oncoscan,
                       annotation_legend_side = "right",
                       annotation_legend_list = list(heatmap_legend),
                       padding = ggplot2::unit(c(5, 2, 2, 2), "cm"))
  
  grDevices::dev.off()
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}
