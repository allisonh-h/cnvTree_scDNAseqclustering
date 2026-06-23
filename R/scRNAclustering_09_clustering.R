#' Superimpose scRNA-seq CNVs onto scDNA-seq CNVs
#'
#' This function integrates high-confidence "Determined CNVs" (typically from 
#' scDNA-seq) with experimental CNV regions (typically from scRNA-seq). 
#' It calculates the genomic overlap ratio for each cell group to determine if a 
#' specific CNV event is functionally present.
#'
#' @param inputFILE A list object containing \code{cnv_region} (data.frame) and 
#'   \code{cnv_grouping} (data.frame) from the infercnv pipeline.
#' @param DeterminedCNVs A data.frame containing the scDNA-seq CNV regions. 
#'   Must include columns: \code{chr}, \code{CN}, \code{CNV_start}, \code{CNV_end}, 
#'   and \code{CNV_region}.
#' @param cnv_ratio A numeric value between 0 and 1 serving as the overlap threshold. 
#' If the ratio of intersection length to total target CNV length meets or exceeds 
#' this value, the variation is flagged as present (\code{final_cnv = 1})
#'
#' @details 
#' The function iterates through each cell group and performs the following:
#' \enumerate{
#'   \item \strong{Intersection}: Matches regions based on Chromosome (\code{chr}) 
#'   and Copy Number state (\code{CN}).
#'   \item \strong{Segment Logic}: Uses a coordinate-based flag system 
#'   (\code{F_start}, \code{F_end}) to identify how the scRNA-seq data overlaps 
#'   with scDNA-seq Data.
#'   \item \strong{Ratio Calculation}: Determines the \code{cnv_ratio} 
#'   (percentage of the scDNA-seq Data covered by scRNA-seq Data).
#'   \item \strong{Final Call}: Assigns \code{final_cnv = 1} if the overlap ratio 
#'   is at least 50% (\eqn{\ge 0.5}), and \code{0} otherwise.
#' }
#'
#' @return The original \code{inputFILE} list with an additional \code{superimpose} 
#'   data.frame containing the overlap statistics and final CNV calls.
#'
superimpose.data <- function(inputFILE, DeterminedCNVs, cnv_ratio)
{
  fucStep <- paste0(" 9.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  Groups <- unique(inputFILE$cnv_grouping$cell_group_name)
  superimpose <- NULL
  
  for (groups in 1:length(Groups)) {
    intersection <- NULL
    # check defined CNVs in each group
    CNVs <- inputFILE$cnv_region %>%
            filter(cell_group_name %in% c(Groups[groups]),
                   chr %in% c(unique(DeterminedCNVs$chr)),
                   CN %in% c(unique(DeterminedCNVs$CN)))

    intersection <- merge(DeterminedCNVs, CNVs, by = c("chr", "CN"))
    intersection <- intersection %>%
                    mutate(F_start = case_when(start < CNV_start ~ 0,
                                               start >= CNV_start & start <= CNV_end ~ 1,
                                               start > CNV_end ~ 2),
                           F_end = case_when(end < CNV_start ~ 0,
                                             end >= CNV_start & end <= CNV_end ~ 1,
                                             end > CNV_end ~ 2),
                           seg = paste0(F_start, F_end),
                           final_start = case_when(seg %in% c("00", "22") ~ NA,
                                                   seg %in% c("01", "02") ~ CNV_start,
                                                   seg %in% c("11", "12") ~ start),
                           final_end = case_when(seg %in% c("00", "22") ~ NA,
                                                 seg %in% c("01", "11") ~ end,
                                                 seg %in% c("02", "12") ~ CNV_end)) %>%
                    filter(!seg %in% c("00", "22")) %>%
                    mutate(cnv_range = final_end - final_start + 1) %>%
                    group_by(CNV_region) %>%
                    summarise(cnv_range = sum(cnv_range)) %>%
                    as.data.frame(.)
    
    superimpose <- left_join(DeterminedCNVs, intersection, by = "CNV_region") %>%
                   mutate(cell_group_name = Groups[groups],
                          CNV_range = CNV_end - CNV_start + 1,
                          cnv_range = ifelse(is.na(cnv_range) == TRUE, 0, cnv_range),
                          Cnv_ratio = cnv_range / CNV_range,
                          final_cnv = ifelse(Cnv_ratio >= cnv_ratio, 1, 0)) %>%
                   rbind(., superimpose)
  }
  inputFILE$superimpose <- superimpose
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(inputFILE)
}


#' Expand Group-Level CNV Calls to Individual Cell Level
#'
#' This function transforms summarized group-level CNV results into a cell-by-CNV
#' binary matrix. It maps the \code{final_cnv} status of each cell group back to
#' every individual cell (\code{cellID}) belonging to that group.
#'
#' @param inputFILE A list object containing the \code{superimpose} data.frame 
#'   (output from \code{superimpose.data}) and the \code{cnv_grouping} mapping.
#' @param DeterminedCNVs A data.frame of the scDNA-seq regions used to 
#'   standardize the column names (CNV regions).
#'
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies unique cell groups present in the superimposed data.
#'   \item For each group, retrieves the list of individual \code{cellID}s.
#'   \item Extracts the binary CNV pattern (\code{final_cnv} vector) for that group.
#'   \item Replicates the group pattern for every cell in that group using \code{sapply}.
#'   \item Combines all groups into a single master matrix using \code{rbind}.
#' }
#'
#' @return A matrix where rows are \code{cellID}s and columns are \code{CNV_regions}. 
#'   Values are binary (1 for presence, 0 for absence of CNV).
#'
superimpose.FileLevel <- function(inputFILE, DeterminedCNVs) 
{
  fucStep <- paste0(" 9.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  Groups <- unique(inputFILE$superimpose$cell_group_name)
  
  output <- NULL
  for (groups in 1:length(Groups)) {
    cellID_list <- inputFILE$cnv_grouping %>% 
                   filter(cell_group_name %in% c(Groups[groups])) %>%
                   pull(cellID)
    cnv_pattern <- inputFILE$superimpose %>%
                   filter(cell_group_name %in% c(Groups[groups])) %>%
                   arrange(CNV_region) %>%
                   pull(final_cnv)
    result_table <- data.frame(matrix(ncol = length(cnv_pattern), 
                                      nrow = length(cellID_list)))
    
    result_table <- t(sapply(1:length(cellID_list), function(i) cnv_pattern))
    rownames(result_table) <- cellID_list
    
    output <- rbind(output, result_table)
  }
  
  colnames(output) <- paste0("CNV", DeterminedCNVs$CNV_region)
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(output)
}