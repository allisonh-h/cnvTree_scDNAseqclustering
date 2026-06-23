#' Convert copy number data to GRanges Format
#'
#' This function reads a `.rds` file or `.txt` file containing copy number 
#' variation (CNV) data and converts it into a list of `GRanges` objects, where 
#' each element corresponds to a single cell.
#'
#' @param input_dir_DNA A `.rds` file or `.txt` file containing a data frame with 
#'  the following required columns:
#'  
#'    - `cellID`: Unique identifier for each cell.
#'    - `seqnames`: Chromosome or sequence name.
#'    - `start`: Start position of the segment.
#'    - `end`: End position of the segment.
#'    - `copy.number`: Copy number value for the segment.
#'    
#' @param cores An integer specifying the number of CPU cores to use for parallel 
#'  processing, default = 1.
#' @return A named list where each element represents a cell, containing its 
#'  corresponding genomic segments as a `GRanges` object.
#'
#' @export
#'
changeFormat <- function(input_dir_DNA, cores, sexchromosome) 
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  fucStep <- paste0(" 1.0_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  # Locked variable
  cores = config_hid$cores
  sexchromosome = config_hid$sexchromosome
  
  if(endsWith(input_dir_DNA, ".rds") == TRUE) {
    Bin_CN <- readRDS(input_dir_DNA)
  } else if(endsWith(input_dir_DNA, ".txt") == TRUE) {
    Bin_CN <- utils::read.delim2(input_dir_DNA, sep = " ")
  } else {
    "Please input either .rds or .txt format as input."
    stop(changeFormat)
  }

  if (sexchromosome == FALSE) {
    rows_to_remove <- grepl("chrX|chrY|chrM", Bin_CN$seqnames)
    Bin_CN <- Bin_CN[!rows_to_remove, ]
  }
  
  data.table::setDT(Bin_CN)
  
  clean_chroms <- unique(str_remove(Bin_CN$seqnames, "(?i)chr"))
  is_numeric <- !is.na(suppressWarnings(as.numeric(clean_chroms)))
  sorted_nums <- sort(as.numeric(clean_chroms[is_numeric]))
  sorted_sex  <- sort(clean_chroms[!is_numeric])
  
  correct_level_order <- c(as.character(sorted_nums), sorted_sex)
  
  Bin_CN <- Bin_CN %>%
            mutate(seqnames_sort = factor(str_remove(seqnames, "(?i)chr"), 
                   levels = correct_level_order)) %>%
            arrange(seqnames_sort)
  # ----------------------------------------------------------------------------
  
  # Set up parallel environment
  future::plan(future::multisession, workers = cores)

  # Split data by "cellID"
  Bin_CN_list <- split(Bin_CN, by = "cellID", keep.by = FALSE)
  template_data <- Bin_CN_list[[1]][,setdiff(names(Bin_CN_list[[1]]), 
                                                   "copy.number"), 
                                                   with = FALSE]
  template_gr <- GenomicRanges::makeGRangesFromDataFrame(template_data, 
                                                         keep.extra.columns = FALSE)
  NewFormat <- list()
  NewFormat <- lapply(seq_along(Bin_CN_list), function(i) {
  #NewFormat <- future.apply::future_lapply(seq_along(Bin_CN_list), function(i) {##LH
    # Extract current Bin_CN_list element
    cell_data <- Bin_CN_list[[i]]

    # Copy the template GRanges object and add "copy.number"
    Bins <- template_gr
    S4Vectors::mcols(Bins)$copy.number <- cell_data$copy.number

    # Compute breakpoints
    shifted_cn <- c(NA, cell_data$copy.number[-nrow(cell_data)])  # 向前平移
    breakpoint_rows <- cell_data$copy.number != shifted_cn & !is.na(shifted_cn)

    breakpoints <- cell_data[breakpoint_rows, 
                             c("seqnames", "start", "end", "copy.number"), 
                             with = FALSE]

    # Convert breakpoints to GRanges if not empty
    Breakpoints <- if (nrow(breakpoints) == 0) {
      NULL
    } else {
      GenomicRanges::GRanges(seqnames = breakpoints$seqnames,
                             ranges = IRanges::IRanges(start = breakpoints$start, 
                                                       end = breakpoints$end),
                             copy.number = breakpoints$copy.number)
    }
    # Return a list with the new format
    list(ID = names(Bin_CN_list)[i], 
         bins = Bins, 
         breakpoints = Breakpoints)
  })

  names(NewFormat) <- names(Bin_CN_list)

  on.exit(future::plan(future::sequential), add = TRUE)
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(NewFormat)
}


#' Generate copy number matrix for selected cells
#'
#' This function extracts copy number variations from a list of `GRanges` objects
#' and organizes them into an integer matrix. The matrix contains selected cells 
#' as columns, with genomic regions (fixed bins) as rows.
#'
#' @param input A named list where each element represents a single cell as a 
#'  `GRanges` object.
#' @param Template A character vector containing the `cellID`s of selected cells 
#'  to be included in the matrix.
#'
#' @return An integer matrix:
#'   - Columns represent the selected `cellID`s.
#'   - Rows represent genomic regions, separated into fixed bins.
#'
NEW_CN_seq <- function(input, Template)
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  fucStep <- paste0(" 1.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  c <- lapply(Template, function(i) {
    obj <- input[[i]]
    if (!is.null(obj) && "bins" %in% names(obj) && 
        "copy.number" %in% names(S4Vectors::mcols(obj[["bins"]]))) {
      S4Vectors::mcols(obj[["bins"]])$copy.number
    } else {
      warning(sprintf("Missing data or 'copy.number' column for template '%s'", i))
      NA_real_  # Use real if copy numbers are decimals
    }
  })
  
  c <- as.data.frame(c, stringsAsFactors = FALSE) #避免字符串變為因子
  names(c) <- Template
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(c)
}


#' Run the cnvTree Setup Pipeline
#' 
#' This initializes the output environment, sets up necessary folder structures 
#' for RNA data, and extracts the available cell groups/types for selection.
#' 
#' @param config List. Internal configuration parameters containing versioning 
#'   and message flags (e.g., \code{config$input_dir_RNA}).
#' @param output_dir Character. Path to the directory where results and output 
#'   folders will be created and saved.
#' @param input_dir_RNA Character. Path to the input single-cell RNA-seq 
#'   data directory (e.g., inferCNV output folder).
#' @param RNAdataSource Character. The format or source type of the RNA data 
#'   (e.g., "dataframe" or "inferCNV").
#'
#' @return A character vector or list of identified cell groups/types extracted 
#'   by \code{select_groups()} for the next steps of the pipeline.
#' 
#' @export 
#' 
run_cnvTree_Pipeline <- function(output_dir,
                                 input_dir_RNA, RNAdataSource) 
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  print(config_path_hid)
  config_hid <- read_yaml(config_path_hid)
  fucStep <- paste0(" 1.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  #Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    message("Created output folder: ", output_dir)
  }
  selected_groups <- select_groups(input_dir_RNA = input_dir_RNA, 
                                   RNAdataSource = RNAdataSource)

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(selected_groups)
}


#' Setup Directory Paths Based on RNA Data Source Type
#'
#' Scans the single-cell RNA input directory to dynamically identify and select 
#' the target data folder based on the structure of the specified data source.
#'
#' @param input_dir_RNA Character. The root input directory path to recursively scan.
#' @param output_dir Character. Path to the output folder.
#' @param RNAdataSource Integer or Character. The type of RNA data source being used:
#'   \itemize{
#'     \item \code{1}: Tabular data file / matrix structure.
#'     \item \code{2}: InferCNV output folder structure.
#'   }
#'
#' @return Character. The absolute file path of the selected target folder.
#' 
setup_pipeline_folders <- function(input_dir_RNA, output_dir, RNAdataSource)
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" 1.2.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)

  all_dirs <- list.dirs(path = input_dir_RNA, 
                        recursive = TRUE, 
                        full.names = TRUE)
  
  dir_depths <- lengths(gregexpr("/", all_dirs))

  #Ensure the folder actually contains sub-directories to prevent errors
  if (length(dir_depths) == 0) {
    stop("Error: No directories found in the specified input_dir_RNA.")
  }
  max_depth <- max(dir_depths)

  #Identify folders at specific depths
  if (RNAdataSource == 1) {
    folders_2nd <- all_dirs[dir_depths == (max_depth)][1]
  } else if (RNAdataSource == 2) {
    folders_2nd <- all_dirs[dir_depths == (max_depth - 1)][1]
  } else {
    stop("Please define RNAdataSource for:
          '1' Tabular data file
          '2' InferCNV output folder")
  }

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(target_folder)
}


#' Selects Cell Groups Interactively from CNV region File
#'
#' This function scans a specified directory for CNV prediction files, extracts 
#' unique cell group names, and provides an interactive menu for the user to 
#' select one or more groups for downstream analysis.
#'
#' @param input_dir_RNA A character string specifying the path to the directory containing 
#'   the CNV results (specifically looking for files containing "cnv_regions" 
#'   and "HMM_CNV_predictions").
#' @param RNAdataSource Integer or Character. The type of RNA data source being used:
#'   \itemize{
#'     \item \code{1}: Tabular data file / matrix structure.
#'     \item \code{2}: InferCNV output folder structure.
#'   }
#'
#' @details 
#' The function identifies the relevant file in the directory, reads the 
#' `cell_group_name` column, and strips suffixes (anything after the first dot) 
#' to identify unique groups. It then opens an interactive menu. 
#' 
#' User Flow:
#' 1. Select an initial group. If '0' is pressed, the function execution stops.
#' 2. A repeat loop allows for multiple additional selections.
#'
#' @return A character vector containing the unique names of the selected cell groups.
#' 
select_groups <- function(input_dir_RNA, RNAdataSource)
{
  config_path_hid <- system.file("cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" 1.2.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)

  if (RNAdataSource == "1") {
    print("datatype: dataframe")
    ext <- tools::file_ext(input_dir_RNA)
    if (ext == "rds") {
      cnv_regions <- readRDS(input_dir_RNA)
    } else if (ext == "txt") {
      cnv_regions <- read.table(input_dir_RNA, sep = " ", header = TRUE) 
    } else {
      stop("Unsupported format. Please provide a .rds, .txt, or .tsv file.")
    }
    
  } else if (RNAdataSource == "2") { # InferCNV data
    print("datatype: InferCNV files")
    all_dirs <- list.dirs(path = input_dir_RNA, recursive = TRUE, full.names = TRUE)
    dir_depths <- lengths(gregexpr("/", all_dirs))
    max_depth <- max(dir_depths)
    folders_pth <- all_dirs[dir_depths == (max_depth - 1)][1]
    File <- list.files(folders_pth, full.names = TRUE)
    File <- File[grepl("cnv_regions", File) & 
                   grepl("HMM_CNV_predictions", File)]
    cnv_regions <- read.table(File, sep = "\t", quote = "", comment.char = "", 
                              fill = TRUE, header = TRUE)
    
  } else {
    stop("Please set `config$RNAdataSource = 1 or 2`")
  }
  cell_group_names <- cnv_regions$cell_group_name
  extracted_cell_group_names <- unique(sub("\\..*", "", cell_group_names))
  
  cat (
    "\ncnvTree Selection:",
    "  1. Enter a group number to select it.",
    "  2. Enter 0 to finish.",
    "  [Constraint]: Do not select the group used as the InferCNV reference.\n", 
    sep = "\n"
  )
  
  menu_choices <- c(extracted_cell_group_names)
  exit_index <- length(menu_choices) # This is the index number R assigns to our exit string
  
  choice_index_0 <- menu(
    choices = menu_choices, 
    title = "Select a group index:"
  )
  
  if (choice_index_0 == 0) {
    stop("No group selected, stopping.") 
  } else {
    print(paste0("First Group Selected: ", extracted_cell_group_names[choice_index_0])) 
  }
  
  repeat {
    n <- length(extracted_cell_group_names)
    choice_index <- menu(
      choices = menu_choices, 
      title = "cnvTree: Choose more groups or enter 0 to finish."
    )
    
    # Exit condition triggered by native 0 or selecting our appended "88" option
    if (choice_index == 0) {
      DebugMsg(fucStep, "end", msg = config_hid$msg)
      break 
      
    } else if (choice_index %in% choice_index_0) {
      # CRITICAL FIX: Wrapped unique() in sort() because identical(c(2,1), c(1,2)) is FALSE.
      # Sorting ensures it successfully catches when all groups are selected.
      if (!identical(sort(unique(choice_index_0)), seq_len(n))) {
        print("cnvTree: Group already selected. Select another group, or enter 0 to finish.") 
      } else {
        print("cnvTree: All groups selected. Select 0 to finish.") 
      }
      
    } else {
      choice_index_0 <- c(choice_index_0, choice_index)
      print(paste0("Group selected: ", extracted_cell_group_names[choice_index])) 
    }
  }
  
  # Return unique selected groups
  selected_groups <- unique(extracted_cell_group_names[choice_index_0])
  message("Final selected groups: ")
  message(paste0(selected_groups, collapse = ", "))

  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(selected_groups)
}


#' Load and Standardize AneuFinder HMM Objects
#'
#' This function imports AneuFinder Hidden Markov Model (HMM) files from a specified 
#' directory. It ensures that the list names and internal object IDs match the file 
#' basenames rather than absolute system paths.
#'
#' @param input_dir_DNA A string specifying the directory containing the \code{.RData} 
#'   or \code{.rds} HMM files (typically the 'method-edivisive' output folder).
#'
#' @return A named list of AneuFinder HMM objects with standardized names and IDs.
#'
ProcessHmmList <- function(input_dir_DNA) 
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" 1.3_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  files <- list.files(input_dir_DNA, full.names = TRUE)
  if (length(files) == 0) stop ("No hmm files found at: ", input_dir_DNA)
  
  hmms <- loadFromFiles(files)
  names(hmms) <- basename(names(hmms)) # Strip Path from Names
  # Sync Internal IDs
  for (i in seq_along(hmms)) {
    hmms[[i]]$ID <- names(hmms)[i]
  }
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(hmms)
}

