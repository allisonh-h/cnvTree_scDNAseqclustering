#' scRNA_input.infercnv
#' 
#' Import and Parse InferCNV Outputs; this function parses InferCNV output 
#' directories to extract CNV regions and cell groupings. It is capable of handling 
#' both single, consolidated output folders and outputs that have been split across 
#' multiple directories. It also performs a quality control check to ensure there 
#' are no duplicated `cellID`s across the merged datasets.
#'
#' @param input_dir_RNA A character string specifying the main directory path 
#'   containing the InferCNV output.
#' @param selected_groups A character vector of group names to retain during 
#'   parsing.
#' @param RNAdataSource 
#'
#' @return A nested list object containing the parsed InferCNV data. The list 
#'   contains a `Round_1` element, which holds sub-lists (one per processed 
#'   file/folder) consisting of:
#'   \describe{
#'     \item{cnv_region}{The parsed CNV regions.}
#'     \item{cnv_grouping}{The parsed cell groupings.}
#'   }
#'
scRNA_input.infercnv <- function(input_dir_RNA, selected_groups, RNAdataSource) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 11.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  cnv_Output <- NULL
  cnv_region <- infercnv_cnvregion(input_dir_RNA = input_dir_RNA, 
                                   selected_groups = selected_groups,
                                   RNAdataSource = RNAdataSource)
  cnv_grouping <- infercnv_cnvgrouping(input_dir_RNA = input_dir_RNA, 
                                       selected_groups = selected_groups,
                                       RNAdataSource = RNAdataSource)
  result_list <- list(cnv_region = cnv_region, 
                      cnv_grouping = cnv_grouping)
  cnv_Output$Round_1 <- result_list 
  
  merged_cnv_groupings <- lapply(cnv_Output, function(Round) {
    round_data <- lapply(Round, function(File) File$cnv_grouping)
    do.call(rbind, round_data)
  })
  
  for(i in 1: length(merged_cnv_groupings)) {
    # Check for duplicated cellID in the i-th data frame of merged_cnv_groupings
    dup <- duplicated(merged_cnv_groupings[[i]]$cellID)
    if (any(dup)) {
      print(paste("LH: duplicated cellID found in index", i)) 
    }
    if (any(duplicated(merged_cnv_groupings[[i]]$cellID))) {
      duplicated_value <- merged_cnv_groupings[[i]]$cellID[duplicated(merged_cnv_groupings[[i]]$cellID)]
      stop(paste("Duplicated cellID in Round", i, ":", duplicated_value))
    }
  }
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(cnv_Output)
}


#' scRNA_superimpose
#' 
#' Batch Superimpose scDNA-seq CNVs onto scRNA-seq Datasets; this wrapper function 
#' processes multiple scRNA-seq output files, superimposing scDNA-seq CNV regions 
#' onto InferCNV observational data. It supports both flat list structures and 
#' nested list structures (split outputs).
#'
#' @param RNA_output A list of objects, where each object contains \code{cnv_region} 
#'  and \code{cnv_grouping}. If \code{splitOutput} is TRUE, this should be a 
#'  nested list of lists.
#' @param DeterminedCNVs A data.frame containing the scDNA-seq CNV regions 
#'  standardized by chromosome and copy number state.
#' @param cnv_ratio A numeric value defining the parameter filtering scale or  
#'  abundance ratio required to validate a superimposed CNV region call.
#'
#' @details 
#'  The function operates in two modes based on the global \code{splitOutput} flag:
#' \itemize{
#'   \item \strong{Single File Mode}: Processes a flat list of files and aggregates 
#'   them into a single output labeled "Round_singleFile".
#'   \item \strong{Split Mode}: Processes a nested list (e.g., from multiple 
#'   experimental runs or time points), aggregating results by "Round_N".
#' }
#' It sequentially calls \code{\link{superimpose.data}} to calculate overlap ratios 
#' and \code{\link{superimpose.FileLevel}} to expand those calls to a cell-by-CNV matrix.
#'
#' @return A named list of data.frames. Each data.frame is a binary matrix where 
#'   rows are cells and columns are CNV regions.
#'
scRNA_superimpose <- function(RNA_output, DeterminedCNVs, cnv_ratio) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 11.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  if (length(RNA_output) == 0) {
    message("Please input RNA_output.")
  } else {
    superimpose_output <- list()
    for (Files in 1:length(RNA_output)) {
      S_output <- NULL
      RNA_output[[Files]] <- superimpose.data(inputFILE = RNA_output[[Files]],                      
                                              DeterminedCNVs = DeterminedCNVs,
                                              cnv_ratio = cnv_ratio)
      Output <- superimpose.FileLevel(inputFILE = RNA_output[[Files]], 
                                      DeterminedCNVs = DeterminedCNVs)
      S_output <- rbind(S_output, Output)
    }
    Folder_name <- paste0("Round_", "noVoting")
    S_output <- as.data.frame(S_output)
    superimpose_output[[Folder_name]] = S_output
  }

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(superimpose_output)
}


#' scRNA_output
#' 
#' Final output coordinator for the scRNA-seq CNV analysis pipeline; this function
#' extracts the clustered cell groupings, formats the summarized RNA cluster data, 
#' and triggers the generation of the final CNV pattern heatmap. All results 
#' are exported directly to a specified output directory.
#'
#' @param Summary A list object containing the final clustering results. It must 
#'  contain a `$Voting_result` element (typically generated by the `scRNA_clustering()` 
#'  function).
#' @param output.dir A character string specifying the directory path where output 
#'  files will be saved.
#' @param DeterminedCNVs A data frame detailing the high-confidence CNV regions. 
#'  This is passed directly to the `scRNA_CNVpattern()` plotting function.
#' @param cellcutoff A numeric value specifying the minimum number of cells required 
#'  for a cluster to be retained.
#' @param cellcutoffRNA A numeric threshold specifying the minimum number of cells 
#'  required to retain an RNA cluster during the formatting step.
#' @param filterZero A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},  
#'  removes columns where the sum of all cell clusters is zero.
#'   
#' @return 
#' This function does not return an R object. Instead, it generates and saves three  
#' files to the `output.dir`:
#' \itemize{
#'   \item \code{cnvTree.scRNAseq_grouping}: A tabular file of cell IDs and their 
#'    CNV patterns.
#'   \item \code{cnvTree.scRNAseq_RNAcluster}: A tabular file of the formatted 
#'    superimpose output.
#'   \item \code{cnvTree.scRNAseq_Fig_CNVpattern.png}: A heatmap visualization 
#'    of the CNV patterns.
#' }
#'
scRNA_output <- function(Summary, output_dir, DeterminedCNVs, cellcutoff, 
                         cellcutoffRNA, filterZero) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 11.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  message("=== Step 05: Output cnvTree scRNA results ===")
  timestamp <- format(Sys.time(), "%m%d_%H")
  
  # 1. cellID clustering output
  scRNA_grouping <- as.data.frame(Summary$Round_noVoting) %>%
                    tibble::rownames_to_column(var = "cellID")
  filename = paste0("/", timestamp, "_", "cnvTree.scRNAseq_grouping")
  writeOutput(data = scRNA_grouping, filename = filename, path = output_dir)
  
  # 2. superimpose output
  scRNA_output <- scRNA_output.format(inputFILE = Summary, 
                                      cellcutoffRNA = cellcutoffRNA,
                                      DeterminedCNVs = DeterminedCNVs,
                                      filterZero = filterZero)
  filename = paste0("/", timestamp, "_", "cnvTree.scRNAseq_RNAcluster")
  writeOutput(data = scRNA_output, filename = filename, path = output_dir)
  
  # 3. RNA cnv plot: final votes
  filename = paste0("/", timestamp, "_", "cnvTree.scRNAseq_Fig_CNVpattern.png")
  CNVpattern(Input = scRNA_output, 
             FILEname = filename,
             FILEpath = output_dir,
             patternType = "scRNA")
  
  # 4. DNA clusters output
  Data_final_DNA <- scDNA_output.format(inputFILE = output_dir, 
                                      cellcutoff = cellcutoff,
                                      filterZero = filterZero)
  filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_Fig_CNVpattern.png")
  CNVpattern(Input = Data_final_DNA, 
             FILEname = filename,
             FILEpath = output_dir,
             patternType = "scDNA")
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}


#' cnvTree_scRNAclustering
#' 
#' Execute the Full cnvTree scRNA-seq Clustering Pipeline; this master wrapper 
#' function sequentially executes the entire scRNA-seq CNV analysis pipeline. 
#' It parses the InferCNV outputs, superimposes high-confidence CNV regions, 
#' performs consensus clustering across multiple iterations (if applicable), and 
#' generates the final tabular files and heatmap visualizations in the specified 
#' directory.
#'
#' @param input_dir_RNA A character string specifying the main directory path containing the 
#'  InferCNV outputs.
#' @param selected_groups A character vector of group names to retain during the 
#'  parsing step.
#' @param output.dir A character string specifying the directory path where the final 
#'  results (tables and plots) will be saved.
#' @param DeterminedCNVs A data frame detailing the high-confidence CNV regions (chromosome, 
#'   copy number state, and cytoband boundaries).
#' @param cellcutoff A numeric value specifying the minimum number of cells required 
#'  for a cluster to be retained.
#' @param cellcutoffRNA A numeric threshold specifying the minimum number of cells 
#'   required to retain an RNA cluster during the formatting step.
#' @param filterZero A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},  
#'   removes columns where the sum of all cell clusters is zero.
#' @param cnv_ratio A numeric value defining the parameter filtering scale or  
#'  abundance ratio required to validate a superimposed CNV region call.
#' @param RNAdataSource Character. The format type of the RNA data source:
#'  \itemize{
#'   \item \code{"1"}: Tabular data file (read via \code{readRDS}).
#'   \item \code{"2"}: InferCNV output folder (searches for "HMM_CNV_predictions" 
#'   and "cnv_regions" files).
#'  }
#'
#' @return 
#' This function does not return an R object to the environment. It acts as an 
#' orchestrator, executing side effects that save three distinct files (grouping  
#' data, superimpose data, and a CNV heatmap) directly to the specified `output.dir`.
#'
#' @export
#' 
cnvTree_scRNAclustering <- function(input_dir_RNA, selected_groups, output_dir, 
                                    cellcutoff, cellcutoffRNA, RNAdataSource, 
                                    cnv_ratio, filterZero) 
{ 
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 11.4_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------

  ptm <- startTimed("Start cnvTree_scRNAclustering...")
  
  scDNA_file <- list.files(path = output_dir, 
                           pattern = "DefinedCNVregion.*\\.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)
  if (length(scDNA_file) == 0) {
    stop("Error: No DefinedCNVregion file found")
  } else if (length(scDNA_file) >= 2) {
    stop("Error: Please provide only one DefinedCNVregion file")
  }
  
  Determine_CNVs <- read.table(scDNA_file, header = TRUE)

  RNA_output <- scRNA_input.infercnv(input_dir_RNA = input_dir_RNA, 
                                     selected_groups = selected_groups,
                                     RNAdataSource = RNAdataSource)

  Superimpose_output <- scRNA_superimpose(RNA_output = RNA_output, 
                                          DeterminedCNVs = Determine_CNVs,
                                          cnv_ratio = cnv_ratio)

  Superimpose_output$Round_noVoting <- Superimpose_output$Round_noVoting %>% 
      mutate(Pattern = apply(., 1, function(row) paste0(row, collapse = "_")))
    
  scRNA_output(Summary = Superimpose_output, 
               output_dir = output_dir, 
               cellcutoff = cellcutoff,
               cellcutoffRNA = cellcutoffRNA,
               DeterminedCNVs = Determine_CNVs,
               filterZero = filterZero)
  endTimed(ptm)
  print(paste("Output directory is set at: ", config$output_dir))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Superimpose_output) #LH changed 042026
}