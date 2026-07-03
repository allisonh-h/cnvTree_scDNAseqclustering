#' scRNA_output.format
#' 
#' Format single-cell RNA copy number variation (CNV) data; this function processes, 
#' filters, and formats single-cell RNA-seq CNV data based on cell cutoffs and 
#' frequency filters. It clusters CNV patterns, renames genomic columns with 
#' readable cytoband notation, adjusts values based on amplification/deletion  
#' status, and optionally filters out columns containing only zeroes.
#'
#' @param inputFILE A list containing a data frame named \code{Round_noVoting}.  
#'   This data frame must contain a \code{Pattern} column (underscore-separated 
#'   binary CNV strings) along with individual CNV state columns.
#' @param cellcutoffRNA An integer specifying the minimum number of single cells  
#'   required to keep a specific CNV pattern.
#' @param filterZero A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},  
#'   removes columns where the sum of all cell clusters is zero.
#' @param DeterminedCNVs A data frame containing metadata for the CNVs. It must  
#'   include the columns: \code{chr} (chromosome name), \code{CN} (type of variant,  
#'   e.g., "amp"), \code{first_band}, and \code{last_band} (cytoband regions).
#'
#' @return A data frame where rows represent identified RNA cell clusters (named  
#'   with cluster IDs and cell counts) and columns represent distinct CNV events.  
#'   Values are adjusted by variant direction (positive for amplifications, 
#'   negative for deletions).
#' 
scRNA_output.format <- function(inputFILE, cellcutoffRNA, filterZero, DeterminedCNVs)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 10.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  inputFILE <- inputFILE
  summary <- inputFILE$Round_noVoting %>%
             select(Pattern) %>% #select(all_of(Pattern))
             table(.) %>%
             as.data.frame(.) %>%
             setNames(c("Pattern", "RNA_Cellnum"))
  output <- tidyr::separate(summary, Pattern, 
                            into = paste0("CNV", 
                                          1:(ncol(inputFILE$Round_noVoting) - 1)), 
                            sep = "_") %>%
            filter(RNA_Cellnum >= cellcutoffRNA)
  output_matrix <- as.matrix(output)
  output_matrix <- ifelse(output_matrix == 1, 1, 0)
  output$sum_cnv <- apply(output_matrix, 1, sum)
  output <- output %>%
            dplyr::arrange(sum_cnv, RNA_Cellnum ) %>%
            dplyr::mutate(RNA_cluster = 1:n())
  
  # --- Data Preparation ---
  Data_clean <- as.data.frame(output)
  rownames(Data_clean) <- paste0(Data_clean$RNA_cluster, " (n = ", Data_clean$RNA_Cellnum, ")")
  Data_clean <- Data_clean %>% 
                select(!c(RNA_Cellnum, sum_cnv, RNA_cluster)) %>% 
                mutate(across(everything(), as.numeric)) 
  custom_colnames <- DeterminedCNVs %>%
                     mutate(chr = sub("^chr", "", chr),
                            CN_type = ifelse(CN == "amp", "+", "-"),
                            chr_cytoband = paste0(CN_type, " ", chr, " (", 
                                                  first_band, "-", last_band, ")")) %>%
                     pull(chr_cytoband)
  # Subset your custom names to match the number of columns
  target_colnames <- custom_colnames[1:ncol(Data_clean)]
  # Find and print the duplicates
  duplicates <- target_colnames[duplicated(target_colnames)]
  print("Duplicated column names; required to fix the logic on fuction 10.1:")
  print(unique(duplicates))
  # Make the subset names unique
  unique_colnames <- make.unique(target_colnames)
  colnames(Data_clean) <- unique_colnames
  
  if (filterZero == TRUE) {
    # 計算總和、過濾 0 的欄位、並移除總和列（保留原本設定好的 rownames）
    Total_row <- Data_clean %>% 
                 summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
    Data <- Data_clean %>% bind_rows(Total_row) %>% as.data.frame() 
    rownames(Data) <- c(rownames(Data_clean), "Total_Sum")
    Data_final_RNA <- Data %>% 
                      select(where(~ last(.x) != 0)) %>% 
                      dplyr::slice(-n()) %>% 
                      as.data.frame()
  } else {
    Data_final_RNA <- Data_clean
  }
  for (i in seq_len(ncol(Data_final_RNA))) {
    sign_char <- substr(colnames(Data_final_RNA)[i], 1, 1)
    if (sign_char == "-") {
      Data_final_RNA[, i] <- as.integer(Data_final_RNA[, i]) * -1
    }
  }
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Data_final_RNA)
}


#' CNVpattern
#' 
#' Generate a Heatmap of CNV Patterns; this function visualizes the presence of 
#' Copy Number Variations (CNVs) across different single-cell DNA sequencing 
#' (scDNA-seq) clusters. It takes cluster data and high-confidence CNV regions, 
#' formats the labels (e.g., adding cell counts to clusters and formatting 
#' cytobands), and generates a highly customized, publication-ready heatmap 
#' using `ComplexHeatmap`. The resulting plot is saved directly to a specified 
#' directory as a PNG file.
#'
#' @param Input A data frame containing scRNA-seq cluster data. It must include 
#'  the columns `RNA_cluster`, `RNA_Cellnum`, and `sum_cnv`. The remaining columns 
#'  should represent the CNV matrix data to be plotted.
#' @param patternType A character string defining the data profiling schema. 
#'  Must be either \code{"scDNA"} or \code{"scRNA"}. This alters the heatmap 
#'  title and neutral background color definitions.
#' @param FILEpath A character string specifying the output directory path where 
#'  the heatmap image should be saved. Must include a trailing slash if not handled 
#'  globally.
#' @param FILEname A character string specifying the desired name of the output 
#'  PNG file (e.g., `"CNV_pattern_heatmap.png"`).
#'
#' @return This function does not return an R object. Instead, it generates and 
#'  saves a `.png` image file containing the CNV heatmap to the specified path.
#'
CNVpattern <- function(Input, FILEpath, FILEname, patternType) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 10.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  # 將資料轉換為文字矩陣; 保持2D行列結構&把 -1, 0, 1 安全地轉成文字給col_fun
  mat_data <- as.matrix(Input)
  mode(mat_data) <- "character" 
  # 抓出純染色體號碼，用來比對何時要劃分界線
  col_names_clean <- colnames(mat_data)
  chr_labels <- gsub("^[-+]\\s+([0-9XY]+)\\s+\\(.*$", "\\1", col_names_clean)
  
  # 確保名稱是文字，對應矩陣裡的 "-1", "0", "1"
  if (patternType == "scDNA") {
    CNVpattern_title <- "Defined CNV regions \n(scDNA clusters)"
    col <- c("-1" = "#8165A3", "0" = "#9BBB59", "1" = "#FFC000")
  } else if (patternType == "scRNA") {
    CNVpattern_title <- "High-confidence CNV regions \n(DefinedCNVs * scRNA clusters)"
    #col <- c("-1" = "#2166ac", "0" = "#EFEDF6", "1" = "#b2182b")
    col <- c("-1" = "#8165A3", "0" = "#EFEDF6", "1" = "#FFC000")
  }
  # 固定網格與 PNG 尺寸動態計算
  cell_size_w <- grid::unit(3, "cm")  
  cell_size_h <- grid::unit(1.0, "cm")
  png_w <- (ncol(mat_data) * 150) + 2500 
  png_h <- (nrow(mat_data) * 80) + 1800
  # 輸出 PNG
  png(filename = paste0(FILEpath, FILEname),
      width = png_w,
      height = png_h,
      res = 100)
  # ---熱圖建構 (使用轉換後的 mat_data) ---
  Oncoscan <- ComplexHeatmap::Heatmap(mat_data,
                                      name = "CNV exist",
                                      width = ncol(mat_data) * cell_size_w,
                                      height = nrow(mat_data) * cell_size_h,
                                      cluster_rows = FALSE,
                                      cluster_columns = FALSE,
                                      show_column_dend = FALSE,
                                      show_row_dend = FALSE,
                                      col = col,
                                      rect_gp = grid::gpar(col = "white", lwd = 2), 
                                      
                                      # 左側列標題
                                      #row_title = "scDNA \nclusters",
                                      #row_title_side = "left",
                                      #row_title_rot = 0, 
                                      #row_title_gp = grid::gpar(fontsize = 32, 
                                      #                          fontface = "bold", 
                                      #                          lineheight = 0.9),
                                      
                                      # 左側列名稱 (細胞群與數量)
                                      row_names_side = "left",
                                      row_names_gp = grid::gpar(fontsize = 28),
                                      show_row_names = TRUE,
                                      
                                      # 上方主要標題
                                      column_title = CNVpattern_title,
                                      column_title_side = "top",
                                      column_title_gp = grid::gpar(fontsize = 40, 
                                                                   fontface = "bold", 
                                                                   gridvjust = 2), 
                                      
                                      # 下方欄位名稱 (染色體與 Cytoband 區間)
                                      column_names_side = "bottom",
                                      column_names_gp = grid::gpar(fontsize = 28),
                                      column_names_rot = 90,
                                      column_names_centered = FALSE,
                                      border = TRUE,
                                      border_gp = grid::gpar(col = "black", lwd = 4), # 加粗外框
                                      
                                      layer_fun = function(j, i, x, y, w, h, fill) {
                                        NC <- ncol(mat_data)
                                        NR <- nrow(mat_data)
                                        if(any(i == NR)) {
                                          for(col_idx in 1:(NC - 1)) {
                                            if(chr_labels[col_idx] != chr_labels[col_idx + 1]) {
                                              this_x <- x[which(j == col_idx & i == 1)]
                                              next_x <- x[which(j == col_idx + 1 & i == 1)]
                                              mid_x <- (this_x + next_x) / 2
                                              grid.lines(
                                                x = grid::unit.c(mid_x, mid_x),
                                                y = grid::unit.c(y[which(j == col_idx & i == NR)] - h[1]/2, 
                                                                 y[which(j == col_idx & i == 1)] + h[1]/2),
                                                gp = grid::gpar(col = "yellow", lwd = 10) # 加粗染色體交界線
                                              )
                                            }
                                          }
                                        }
                                      },
                                      
                                      # 圖例 (Legend) 調整
                                      heatmap_legend_param = list(
                                        title_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                        labels_gp = grid::gpar(fontsize = 20),
                                        grid_width = grid::unit(1.5, "cm"),
                                        grid_height = grid::unit(1.5, "cm"),
                                        at = c("-1", "0", "1"),
                                        labels = c("Del", "Neu", "Amp")
                                      )
  )
  # 加大 padding 確保四週巨大的字體不會超出圖片邊界
  draw(Oncoscan, padding = grid::unit(c(3, 3, 3, 3), "cm"))
  
  invisible(dev.off())
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}


#' Format Single-Cell DNA Copy Number Variant (CNV) Cluster Data
#'
#' @description
#' Reads CNV regions and single-cell DNA subclone cluster data from an output_dir,
#' merges the datasets, formats genomic cytobands, adjusts values based on copy 
#' number types (amplifications vs. deletions), and optionally filters out zero-sum 
#' features.
#'
#' @param inputFILE Character. The directory path where the source text files are located.
#' @param cellcutoff Numeric. The minimum cell count threshold required to retain a 
#'   clone.
#' @param filterZero Logical. If \code{TRUE}, columns (CNV regions) with a total 
#'   sum of 0 across all subclones will be removed.
#'
#' @return A matrix (or data frame structure) where rows represent subclones and 
#'   columns represent formatted genomic regions (cytobands) containing copy number 
#'   values.
#' 
scDNA_output.format <- function(inputFILE, cellcutoff, filterZero) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 10.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------- start below ---------------------
  
  inputFILE <- inputFILE
  # 1. Load Defined CNV Regions
  scDNA_file <- list.files(path = inputFILE, 
                           pattern = "DefinedCNVregion.*\\.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)
  if (length(scDNA_file) == 0) stop("No DefinedCNVregion file found.")
  Determine_CNVs <- read.table(scDNA_file[1], header = TRUE) %>% 
    mutate(CNV_region = as.character(CNV_region),
           chr_index  = sub("chr", "", chr))
  
  # 2. Load and filter DNA Cluster Data
  scDNAcluster_file <- list.files(path = inputFILE, 
                                  pattern = "DNAcluster.*\\.txt$", 
                                  recursive = TRUE, 
                                  full.names = TRUE)
  if (length(scDNAcluster_file) == 0) stop("No DNAcluster file found.")
  scDNAcluster <- read.table(scDNAcluster_file[1], header = TRUE) %>% 
    filter(DNA_Cellnum >= cellcutoff) %>%
    mutate(Subclone_no = paste0(DNA_cluster, " n = (", DNA_Cellnum, ")"))
  
  # 3. Clean and reshape cluster data
  scDNAcluster_clean <- scDNAcluster %>% 
    select(-DNA_cluster, -DNA_Cellnum) %>% 
    pivot_longer(cols = -Subclone_no, 
                 names_to = "CNV_region", 
                 values_to = "Value") %>% 
    mutate(CNV_region = sub("CNV", "", CNV_region),
           Value = as.numeric(Value)) %>% 
    pivot_wider(names_from = Subclone_no, values_from = Value)
  
  # 4. Merge datasets and construct cytoband headers
  combined_df <- full_join(Determine_CNVs, scDNAcluster_clean, by = "CNV_region") %>% 
    mutate(CN_type = if_else(CN == "amp", "+", "-"),
           chr_cytoband = paste0(CN_type, " ", chr_index, 
                                 " (", first_band, "-", last_band, ") "))
  new_rownames <- combined_df$chr_cytoband
  is_negative_row <- combined_df$CN_type == "-"
  meta_cols <- c("CNV_region", "first_band", "last_band", "CNV_start", 
                 "CNV_end", "chr", "CN", "chr_index", "CN_type", "chr_cytoband")
  matrix_data <- combined_df %>% select(-any_of(meta_cols)) %>% as.matrix()
  rownames(matrix_data) <- new_rownames
  # Invert sign for deletion regions (-)
  matrix_data[is_negative_row, ] <- matrix_data[is_negative_row, ] * -1
  final_matrix <- t(matrix_data)
  
  # 5. Filter out columns where the sum of values is 0
  if (filterZero) {
    non_zero_cols <- colSums(final_matrix, na.rm = TRUE) != 0
    Data_final_DNA <- final_matrix[, non_zero_cols, drop = FALSE]
  } else {
    Data_final_DNA <- final_matrix
  }
  
  return(Data_final_DNA)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}