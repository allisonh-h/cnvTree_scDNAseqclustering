#' infercnv_cnvregion
#' 
#' Filter and Format InferCNV HMM Predictions; this function reads InferCNV 
#' HMM prediction files and extracts CNV regions for a specified set of cell groups. 
#' It filters out neutral states (state 3) and calculates the size and type 
#' (amplification or deletion) for each region.
#'
#' @param input_dir_RNA A character string specifying the directory path containing 
#'  the InferCNV output files (specifically matching "cnv_regions" and 
#'  "HMM_CNV_predictions") or a tabular file.
#' @param selected_groups A character vector of group names to filter for. 
#' @param RNAdataSource Character. The format type of the RNA data source:
#'  \itemize{
#'   \item \code{"1"}: Tabular data file (read via \code{readRDS}).
#'   \item \code{"2"}: InferCNV output folder (searches for "HMM_CNV_predictions" 
#'   and "cnv_regions" files).
#'  }
#'   
#' @details 
#' The function performs the following transformations:
#' \itemize{
#'   \item Filters out genomic regions where `state == 3` (neutral/diploid).
#'   \item Calculates `CNV_size` as the absolute difference between start and end 
#'    positions.
#'   \item Categorizes regions into `CN` (Copy Number) types: 
#'         "amp" (amplification, state > 3) or 
#'         "del" (deletion, state < 3).
#'   \item Cleans `cell_group_name` by removing the selection prefix.
#' }
#'
#' @return A `data.frame` (or `tibble`) containing filtered CNV regions with 
#' additional columns: `CNV_size` and `CN`.
#' 
infercnv_cnvregion <- function(input_dir_RNA, selected_groups, RNAdataSource)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 8.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  if (RNAdataSource == "1") {
    cnv_regions <- readRDS(input_dir_RNA)
    
  } else if (RNAdataSource == "2") {
    #======= new =======
    matching_files <- list.files(path = input_dir_RNA, 
                                 pattern = "HMM_CNV_predictions.*cnv_regions", 
                                 recursive = TRUE, full.names = TRUE)
    if (length(matching_files) == 0) {
      stop("No matching CNV prediction files found in: ", input_dir_RNA)
    }
    
    if (length(matching_files) > 1) {
      # Count path separators (handles both '/' and '\')
      file_depths <- lengths(gregexpr("[/\\\\]", matching_files))
      RNA_file <- matching_files[which.max(file_depths)]
    } else {
      RNA_file <- matching_files[1]
    }
    cnv_regions <- read.table(RNA_file, sep = "\t", quote = "", comment.char = "", 
                              fill = TRUE, header = TRUE)
    
    cat("Successfully loaded file from:\n", RNA_file, "\n")
    
    #======= new =======
    #======= old =======
    #all_dirs <- list.dirs(path = input_dir_RNA, 
    #                      recursive = TRUE, 
    #                      full.names = TRUE)
    #dir_depths <- lengths(gregexpr("/", all_dirs))
    #max_depth <- max(dir_depths)
    #folders_2nd <- all_dirs[dir_depths == (max_depth - 1)][1] # InferCNV files
    #File <- list.files(folders_2nd, full.names = TRUE)
    #File <- File[grepl("cnv_regions", File) & grepl("HMM_CNV_predictions", File)]

    #cnv_regions <- read.table(File, sep = "\t", 
    #                          quote = "", 
    #                          comment.char = "", 
    #                          fill = TRUE, 
    #                          header = TRUE)
    #======= old =======
  } else {
    stop("Please define RNAdataSource for:
          '1' Tabular data file
          '2' InferCNV output folder")
  }
  if (length(selected_groups) == 1) {
    cnv_region <- cnv_regions %>% 
                  filter(str_detect(cell_group_name, paste0("^(", paste(selected_groups, collapse = "|"), ")")), 
                         state != 3) %>% 
                  mutate(cell_group_name = str_replace(
                      cell_group_name, 
                      paste0("^(", stringr::str_flatten(selected_groups, collapse = "|"), ")\\."), 
                             ""),
                  CNV_size = abs(start - end),
                  CN = ifelse(state > 3, "amp", "del"))
  } else { # LH added 012025 --->>
    cnv_regionlist <- vector("list", length(selected_groups))
    cnv_regionlist <- vector("list", length(selected_groups))
    
    for (i in seq_along(selected_groups)) {
      selected_groups_index <- selected_groups[i]
      cnv_regionlist[[i]] <- cnv_regions %>% 
        filter(
          str_detect(cell_group_name, paste0("^", selected_groups_index)), 
          state != 3
        ) %>% 
        mutate(
          cell_group_name = str_replace(
            cell_group_name, 
            paste0("^(", stringr::str_flatten(selected_groups, collapse = "|"), ")\\."), 
            ""
          ),
          CNV_size = abs(start - end),
          CN_RNA = ifelse(state > 3, "amp", "del"))
    }
    cnv_region <- bind_rows(cnv_regionlist)
  } # <<--- LH added 012025 
    check_dims(x = cnv_region) ## LH added 012026
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(cnv_region)
}


#' check_dims
#' 
#' checks Dimensions of an Object; this function checks whether an object has 
#' valid dimensions. It stops execution if the object has no dimensions or has 
#' zero rows
#'
#' @param x An R Object to check
#' @param name
#'
check_dims <- function(x, name = "selected_groups") 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 8.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  d <- dim(x)
  if (is.null(d)) {
    stop(name, " has no dimensions")
  }
  if (d[1] == 0) {
    stop(name, " has 0 rows (", d[1], " x ", d[2], ")", 
         "; selectied_groups not presented in the files")
  }
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}


#' infercnv_cnvgrouping
#' 
#' Maps individual cells to Selected CNV Groups; this function parses the 
#' \code{infercnv.observation_groupings.txt} file to create a mapping between 
#' individual cell IDs and their assigned CNV groups. It can return 'all 
#' observations' or 'filter for specific user-selected groups'.
#'
#' @param input_dir_RNA A character string specifying the directory path containing 
#'  the InferCNV output files (specifically matching "cnv_regions" and 
#'  "HMM_CNV_predictions") or a tabular file.
#' @param selected_groups A character vector of group names to filter for.
#' @param RNAdataSource Character. The format type of the RNA data source:
#'   \itemize{
#'     \item \code{"1"}: Tabular data file (read via \code{readRDS}).
#'     \item \code{"2"}: InferCNV output folder (searches for "HMM_CNV_predictions" 
#'     and "cnv_regions" files).
#'   }
#'
#' @details 
#' The function operates in two modes:
#' \enumerate{
#'   \item \strong{Global Mode}: If "all_observations" is present in \code{selected_groups}, 
#'   it processes the entire groupings file.
#'   \item \strong{Filtered Mode}: It uses regex to retain only cells whose 
#'   \code{Dendrogram.Group} matches the \code{selected_groups} prefix.
#' }
#' It converts the original row names (cell bar codes) into a dedicated \code{cellID} column.
#'
#' @return A \code{tibble} or \code{data.frame} with columns:
#' \itemize{
#'   \item \code{cellID}: The unique identifier/barcode for each cell.
#'   \item \code{Dendrogram.Group}: The original group assignment from InferCNV.
#'   \item \code{cell_group_name}: A copy of the dendrogram group for downstream compatibility.
#' }
#'
infercnv_cnvgrouping <- function(input_dir_RNA, selected_groups, RNAdataSource) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 8.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  if (RNAdataSource == "1") {
    cnv_groupings <- utils::read.delim2(input_dir_RNA, sep = " ")
  } else if (RNAdataSource == "2") {
    #======= new =======
    matching_files <- list.files(path = input_dir_RNA, 
                                 pattern = "infercnv.observation_groupings.txt", 
                                 recursive = TRUE, full.names = TRUE)
    if (length(matching_files) == 0) {
      stop("No /infercnv.observation_groupings.txt found in: ", input_dir_RNA)
    }
    
    if (length(matching_files) > 1) {
      # Count path separators (handles both '/' and '\')
      file_depths <- lengths(gregexpr("[/\\\\]", matching_files))
      RNA_file2 <- matching_files[which.max(file_depths)]
    } else {
      RNA_file2 <- matching_files[1]
    }
    cnv_groupings <- read.table(RNA_file2, header = TRUE)
    #======= new =======
    #======= old =======
    #all_dirs <- list.dirs(path = input_dir_RNA, 
    #                      recursive = TRUE, 
    #                      full.names = TRUE)
    #dir_depths <- lengths(gregexpr("/", all_dirs))
    #max_depth <- max(dir_depths)
    #folders_2nd <- all_dirs[dir_depths == (max_depth - 1)][1] # InferCNV files
    #File <- paste0(folders_2nd, "/infercnv.observation_groupings.txt")
    #cnv_groupings <- read.table(File, header = TRUE)
    #======= old =======
  } else {
    stop("Please define RNAdataSource for:
          '1' Tabular data file
          '2' InferCNV output folder")
  }
  if (c("all_observations") %in% selected_groups) { # LH 012025 added --->>
    cnv_grouping <- cnv_groupings %>% 
                    tibble::rownames_to_column(., "cellID") %>%
                    mutate(cell_group_name = Dendrogram.Group)
  } else { # <<--- LH 012025 added
    
    cnv_grouping <- cnv_groupings %>%
                    #Filter if EITHER column matches the selected_groups pattern
                    filter(if_any(any_of(c("Dendrogram.Group", "cell_group.name")), 
                    ~ str_detect(.x, paste0("^", selected_groups, collapse = "|")))) %>%
                    tibble::rownames_to_column("cellID") %>%
                    mutate(cell_group_name = coalesce(
                      if ("Dendrogram.Group" %in% names(.)) Dendrogram.Group else NULL,
                      if ("cell_group.name" %in% names(.)) cell_group.name else NULL))
  }
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(cnv_grouping)
}


#' infercnv_groups_summary
#' 
#' Summarize InferCNV Observation Groupings; this function reads the 
#' `infercnv.observation_groupings.txt` file from a  specified InferCNV output 
#' directory. It calculates the total number of cells assigned to each dendrogram 
#' group and filters the final summary to include only the specified groups.
#'
#' @param input_dir_RNA A character string specifying the directory path containing 
#'  the InferCNV output files (specifically matching "cnv_regions" and 
#'  "HMM_CNV_predictions") or a tabular file.
#' @param selected_groups A character vector of group names to filter the summary by. 
#'  
#' @return A data frame containing two columns:
#'   \describe{
#'     \item{group}{The name of the dendrogram group.}
#'     \item{cell_num}{The total number of cells belonging to that group.}
#'   }
#'
infercnv_groups_summary <- function(input_dir_RNA, selected_groups) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" infercnv_groups_summary_", cnvTree_v_num) # function 8.3
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  File <- paste0(input_dir_RNA, "/infercnv.observation_groupings.txt")
  cnv_groupings <- read.table(File, header = TRUE)
  cnv_grouping_sum <- table(cnv_groupings$Dendrogram.Group) %>%
                      as.data.frame(.) %>%
                      setNames(c("group", "cell_num")) %>%
                      filter(grepl(paste(selected_groups, collapse = "|"), group))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(cnv_grouping_sum)
}
